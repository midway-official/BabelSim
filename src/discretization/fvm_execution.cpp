#include "internal/compute_backend.h"
#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "internal/boundary_evaluation.h"
#include "internal/fvm_execution.h"

#include "babelsim/operators.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {
namespace detail {
namespace {

template <typename T>
void requireCellField(const Field<T>& field, const Mesh& mesh, const char* name) {
    field.validateStorage();
    if (&field.mesh() != &mesh || field.location() != FieldLocation::Cell) {
        throw std::invalid_argument(std::string(name) + " must be a cell field on the run mesh");
    }
}

void requireFaceField(const ScalarField& field, const Mesh& mesh, const char* name) {
    field.validateStorage();
    if (&field.mesh() != &mesh || field.location() != FieldLocation::Face) {
        throw std::invalid_argument(std::string(name) + " must be a face field on the run mesh");
    }
}

void requireDistinct(const void* input, const void* result) {
    if (input == result)
        throw std::invalid_argument("explicit operator input and result must not alias");
}

}  // 匿名命名空间

struct FvmExecution::Implementation {
    Implementation(const Mesh& mesh_value, std::unique_ptr<ComputeBackend> backend_value)
        : mesh(&mesh_value),
          backend(std::move(backend_value)),
          gradient_workspace(mesh_value, FieldLocation::Cell, "grad"),
          face_coefficient_workspace(mesh_value, FieldLocation::Face, "faceCoefficient"),
          face_flux_workspace(mesh_value, FieldLocation::Face, "faceFlux")
    {
        mesh->validate();
        if (meshData(*mesh).ghost_layers != 0 && meshData(*mesh).ghost_layers < 3)
            throw std::invalid_argument("FVM composite operators require ghostLayers >= 3");
        if (!backend) throw std::invalid_argument("FVM execution requires a compute backend");
    }

    template <typename T>
    void synchronize(Field<T>& field) {
        backend->synchronize(field);
    }

    // 通量增量作用于选定几何区域；物理边界不等于并行分区边界。
    void addFaceIncrement(ScalarField& target, double factor, math::FaceRegion region) {
        requireFaceField(target, *mesh, "flux increment target");
        if (region != math::FaceRegion::All && region != math::FaceRegion::Interior)
            throw std::invalid_argument("unknown face region");
        for (Index face : detail::meshData(*mesh).owned_faces) {
            if (region == math::FaceRegion::Interior &&
                detail::meshData(*mesh).face_neighbour[face] == invalid_index) continue;
            detail::fieldData(target)[face] += factor * detail::fieldData(face_flux_workspace)[face];
        }
        synchronize(target);
    }

    const Mesh* mesh;
    std::unique_ptr<ComputeBackend> backend;
    VectorField gradient_workspace;
    ScalarField face_coefficient_workspace;
    ScalarField face_flux_workspace;
};


FvmExecution::FvmExecution(const Mesh& mesh,
                           std::unique_ptr<ComputeBackend> backend)
    : m_implementation(std::make_unique<Implementation>(
          mesh, std::move(backend))) {}
FvmExecution::~FvmExecution() = default;
ComputeBackend& FvmExecution::backend() { return *m_implementation->backend; }
const Mesh& FvmExecution::mesh() const { return *m_implementation->mesh; }

double FvmExecution::relativeChange(
    const VectorField& current,
    const VectorField& previous) const
{
    const Implementation& state = *m_implementation;
    requireCellField(current, *state.mesh, "relative-change current");
    requireCellField(previous, *state.mesh, "relative-change previous");
    double difference_squared = 0.0;
    double current_squared = 0.0;
    for (Index cell : detail::meshData(*state.mesh).owned_cells) {
        difference_squared += squaredNorm(detail::fieldData(current)[cell] - detail::fieldData(previous)[cell]);
        current_squared += squaredNorm(detail::fieldData(current)[cell]);
    }
    const double local[2] = {difference_squared, current_squared};
    double global[2]{};
    state.backend->sum(local, global, 2);
    return std::sqrt(global[0]) / std::max(std::sqrt(global[1]), 1e-30);
}

double FvmExecution::relativeChange(
    const ScalarField& current,
    const ScalarField& previous) const
{
    const Implementation& state = *m_implementation;
    requireCellField(current, *state.mesh, "relative-change current");
    requireCellField(previous, *state.mesh, "relative-change previous");
    double difference_squared = 0.0;
    double current_squared = 0.0;
    for (Index cell : detail::meshData(*state.mesh).owned_cells) {
        const double difference = detail::fieldData(current)[cell] - detail::fieldData(previous)[cell];
        difference_squared += difference * difference;
        current_squared += detail::fieldData(current)[cell] * detail::fieldData(current)[cell];
    }
    const double local[2] = {difference_squared, current_squared};
    double global[2]{};
    state.backend->sum(local, global, 2);
    return std::sqrt(global[0]) / std::max(std::sqrt(global[1]), 1e-30);
}

double FvmExecution::relativeMagnitude(
    const ScalarField& value,
    const ScalarField& reference) const
{
    const Implementation& state = *m_implementation;
    requireCellField(value, *state.mesh, "relative-magnitude value");
    requireCellField(reference, *state.mesh, "relative-magnitude reference");
    double value_squared = 0.0;
    double reference_squared = 0.0;
    for (Index cell : detail::meshData(*state.mesh).owned_cells) {
        value_squared += detail::fieldData(value)[cell] * detail::fieldData(value)[cell];
        reference_squared += detail::fieldData(reference)[cell] * detail::fieldData(reference)[cell];
    }
    const double local[2] = {value_squared, reference_squared};
    double global[2]{};
    state.backend->sum(local, global, 2);
    return std::sqrt(global[0]) / std::max(std::sqrt(global[1]), 1e-30);
}

FluxBalance FvmExecution::fluxBalance(const ScalarField& face_flux) const {
    const Implementation& state = *m_implementation;
    requireFaceField(face_flux, *state.mesh, "face flux");
    FluxBalance result;
    double squared = 0.0;
    double scale = 0.0;
    for (Index cell : detail::meshData(*state.mesh).owned_cells) {
        double imbalance = 0.0;
        for (Index face : state.mesh->cellFaces(cell)) {
            const std::size_t index = static_cast<std::size_t>(face);
            const double outward = detail::meshData(*state.mesh).face_owner[index] == cell
                ? detail::fieldData(face_flux)[face] : -detail::fieldData(face_flux)[face];
            imbalance += outward;
            scale += std::abs(outward);
        }
        result.l1 += std::abs(imbalance);
        squared += imbalance * imbalance;
        result.maximum = std::max(result.maximum, std::abs(imbalance));
    }
    const double local[3] = {result.l1, squared, scale};
    double global[3]{};
    state.backend->sum(local, global, 3);
    const double local_maximum = result.maximum;
    state.backend->maximum(&local_maximum, &result.maximum, 1);
    result.l1 = global[0];
    result.l2 = std::sqrt(global[1]);
    result.relative = global[0] / std::max(global[2], 1e-30);
    return result;
}

bool FvmExecution::all(bool local_condition) const {
    return m_implementation->backend->all(local_condition);
}

PerformanceCounters FvmExecution::performance() const {
    return m_implementation->backend->performance();
}

void FvmExecution::evaluate(math::ScalarGradient operation, VectorField& result) {
    Implementation& state = *m_implementation;
    operation.options.requireComplete("scalar gradient");
    requireCellField(operation.field, *state.mesh, "gradient input");
    requireCellField(result, *state.mesh, "gradient result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    gradient(operation.field, result, operation.options.gradient.value());
    state.synchronize(result);
    if (result.calculatedBoundary()) {
        for (Index face = 0; face < state.mesh->faceCount(); ++face) {
            if (!state.mesh->boundaryFace(face)) continue;
            const Index owner = state.mesh->owner(face);
            const Vec3 n = state.mesh->faceNormal(face);
            const Vec3 d = state.mesh->faceCentre(face) - state.mesh->cellCentre(owner);
            Vec3 g = fieldData(static_cast<const VectorField&>(result))[owner];
            const auto bc = FieldAccess::condition(operation.field, face);
            const double sn = bc.type == BoundaryType::FixedValue
                ? (bc.value - fieldData(operation.field)[owner] - dot(g, d - dot(d,n)*n)) / dot(d,n)
                : bc.type == BoundaryType::FixedGradient ? bc.value : 0.0;
            g += (sn - dot(g,n)) * n;
            FieldAccess::setTrace(result, face, g);
        }
    }

}

void FvmExecution::evaluate(math::VectorGradient operation, TensorField& result) {
    Implementation& state = *m_implementation;
    operation.options.requireComplete("vector gradient");
    requireCellField(operation.field, *state.mesh, "gradient input");
    requireCellField(result, *state.mesh, "gradient result");
    state.synchronize(const_cast<VectorField&>(operation.field));
    gradient(operation.field, result, operation.options.gradient.value());
    state.synchronize(result);
    if (result.calculatedBoundary()) {
        for (Index face = 0; face < state.mesh->faceCount(); ++face) {
            if (!state.mesh->boundaryFace(face)) continue;
            const Index owner = state.mesh->owner(face);
            const Vec3 n = state.mesh->faceNormal(face);
            const Vec3 d = state.mesh->faceCentre(face) - state.mesh->cellCentre(owner);
            Tensor3 g = fieldData(static_cast<const TensorField&>(result))[owner];
            const auto bc = FieldAccess::condition(operation.field, face);
            if (bc.type == BoundaryType::Symmetry) g = symmetricBoundaryValue(g, n);
            else for (int component = 0; component < 3; ++component) {
                const double sn = bc.type == BoundaryType::FixedValue
                    ? (bc.value[component] - fieldData(operation.field)[owner][component] -
                       dot(g[component], d - dot(d,n)*n)) / dot(d,n)
                    : bc.type == BoundaryType::FixedGradient ? bc.value[component] : 0.0;
                g[component] += (sn - dot(g[component],n)) * n;
            }
            FieldAccess::setTrace(result, face, g);
        }
    }

}

void FvmExecution::evaluate(math::FaceFlux operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    operation.options.requireComplete("face flux");
    requireFaceField(result, *state.mesh, "flux result");
    state.synchronize(const_cast<VectorField&>(operation.velocity));
    flux(operation.velocity, result,
         operation.options.interpolation.value(),
         operation.options.gradient.value());
    state.synchronize(result);
    if (operation.velocity.location() == FieldLocation::Cell)
        const_cast<VectorField&>(operation.velocity).setBoundaryFlux(result);
}

void FvmExecution::evaluate(math::FaceDivergence operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireFaceField(operation.flux, *state.mesh, "divergence input");
    requireCellField(result, *state.mesh, "divergence result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    divergence(operation.flux, result);
    state.synchronize(result);
    FieldAccess::extrapolateTrace(result);
}

void FvmExecution::evaluate(math::VectorDivergence operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    operation.options.requireComplete("vector divergence");
    requireCellField(operation.field, *state.mesh, "divergence input");
    requireCellField(result, *state.mesh, "divergence result");
    state.synchronize(const_cast<VectorField&>(operation.field));
    divergence(
        operation.field, result, operation.options.interpolation.value(),
        operation.options.gradient.value());
    state.synchronize(result);
    FieldAccess::extrapolateTrace(result);
}

void FvmExecution::evaluate(math::TensorDivergence operation, VectorField& result) {
    Implementation& state = *m_implementation;
    operation.options.requireComplete("tensor divergence");
    requireCellField(operation.field, *state.mesh, "tensor divergence input");
    requireCellField(result, *state.mesh, "tensor divergence result");
    state.synchronize(const_cast<TensorField&>(operation.field));
    divergence(operation.field, result, operation.options.interpolation.value(),
               operation.options.gradient.value());
    state.synchronize(result);
    FieldAccess::extrapolateTrace(result);
}

void FvmExecution::evaluate(math::ScalarConvection operation, ScalarField& result) {
    requireDistinct(&operation.field, &result);
    Implementation& state = *m_implementation;
    operation.options.requireComplete("scalar convection");
    requireFaceField(operation.flux, *state.mesh, "convection flux");
    requireCellField(operation.field, *state.mesh, "convection field");
    requireCellField(result, *state.mesh, "convection result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    const_cast<ScalarField&>(operation.field).setBoundaryFlux(operation.flux);
    state.synchronize(const_cast<ScalarField&>(operation.field));
    convection(
        operation.flux, operation.field, result,
        operation.options.convection.value(),
        operation.options.interpolation.value(),
        operation.options.gradient.value());
    state.synchronize(result);
    FieldAccess::extrapolateTrace(result);
}

void FvmExecution::evaluate(math::VectorConvection operation, VectorField& result) {
    requireDistinct(&operation.field, &result);
    Implementation& state = *m_implementation;
    operation.options.requireComplete("vector convection");
    requireFaceField(operation.flux, *state.mesh, "convection flux");
    requireCellField(operation.field, *state.mesh, "convection field");
    requireCellField(result, *state.mesh, "convection result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    const_cast<VectorField&>(operation.field).setBoundaryFlux(operation.flux);
    state.synchronize(const_cast<VectorField&>(operation.field));
    convection(
        operation.flux, operation.field, result,
        operation.options.convection.value(),
        operation.options.interpolation.value(),
        operation.options.gradient.value());
    state.synchronize(result);
    FieldAccess::extrapolateTrace(result);
}

void FvmExecution::evaluate(math::ScalarInterpolation operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    operation.options.requireComplete("scalar interpolation");
    requireCellField(operation.field, *state.mesh, "interpolation input");
    requireFaceField(result, *state.mesh, "interpolation result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    interpolate(
        operation.field, result, operation.options.interpolation.value(),
        operation.options.gradient.value());
    state.synchronize(result);
}

void FvmExecution::evaluate(math::VectorInterpolation operation, VectorField& result) {
    Implementation& state = *m_implementation;
    operation.options.requireComplete("vector interpolation");
    requireCellField(operation.field, *state.mesh, "interpolation input");
    if (&result.mesh() != state.mesh || result.location() != FieldLocation::Face) {
        throw std::invalid_argument("vector interpolation result must be a face field");
    }
    state.synchronize(const_cast<VectorField&>(operation.field));
    interpolate(
        operation.field, result, operation.options.interpolation.value(),
        operation.options.gradient.value());
    state.synchronize(result);
}

void FvmExecution::evaluate(math::ScalarReconstruction operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "reconstruction input");
    requireCellField(operation.gradient, *state.mesh, "reconstruction gradient");
    requireFaceField(result, *state.mesh, "reconstruction result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    state.synchronize(const_cast<VectorField&>(operation.gradient));
    reconstruct(operation.field, operation.gradient, result);
    state.synchronize(result);
}

void FvmExecution::evaluate(math::VectorReconstruction operation, VectorField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "reconstruction input");
    requireCellField(operation.gradient, *state.mesh, "reconstruction gradient");
    if (&result.mesh() != state.mesh || result.location() != FieldLocation::Face) {
        throw std::invalid_argument("vector reconstruction result must be a face field");
    }
    state.synchronize(const_cast<VectorField&>(operation.field));
    state.synchronize(const_cast<TensorField&>(operation.gradient));
    reconstruct(operation.field, operation.gradient, result);
    state.synchronize(result);
}

void FvmExecution::evaluate(math::ScalarLaplacian operation, ScalarField& result) {
    requireDistinct(&operation.field, &result);
    requireDistinct(operation.coefficient_field, &result);
    Implementation& state = *m_implementation;
    operation.options.requireComplete("scalar laplacian");
    requireCellField(operation.field, *state.mesh, "laplacian input");
    requireCellField(result, *state.mesh, "laplacian result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    if (operation.coefficient_field == nullptr) {
        laplacian(
            operation.coefficient, operation.field, result,
            operation.options.gradient.value(),
            operation.options.diffusion.value());
    } else {
        const ScalarField& coefficient = *operation.coefficient_field;
        if (&coefficient.mesh() != state.mesh ||
            (coefficient.location() != FieldLocation::Cell &&
             coefficient.location() != FieldLocation::Face)) {
            throw std::invalid_argument("laplacian diffusivity must be a cell or face field");
        }
        state.synchronize(const_cast<ScalarField&>(coefficient));
        if (coefficient.location() == FieldLocation::Cell) {
            interpolate(
                coefficient, state.face_coefficient_workspace,
                operation.options.coefficientInterpolation.value_or(operation.options.interpolation.value()),
                operation.options.coefficientGradient.value_or(operation.options.gradient.value()));
            state.synchronize(state.face_coefficient_workspace);
            laplacian(
                state.face_coefficient_workspace, operation.field, result,
            operation.options.gradient.value(),
            operation.options.diffusion.value());
        } else {
            laplacian(
                coefficient, operation.field, result,
            operation.options.gradient.value(),
            operation.options.diffusion.value());
        }
    }
    state.synchronize(result);
    FieldAccess::extrapolateTrace(result);
}

void FvmExecution::subtract(
    const ScalarField& coefficient,
    math::ScalarGradient operation,
    VectorField& target)
{
    Implementation& state = *m_implementation;
    requireCellField(coefficient, *state.mesh, "gradient multiplier");
    requireCellField(operation.field, *state.mesh, "gradient input");
    requireCellField(target, *state.mesh, "gradient correction target");
    state.synchronize(const_cast<ScalarField&>(coefficient));
    state.synchronize(const_cast<ScalarField&>(operation.field));
    state.gradient_workspace.useCalculatedBoundary();
    evaluate(operation, state.gradient_workspace);
    target.addProduct(-1.0, coefficient, state.gradient_workspace);
    state.synchronize(target);
}

void FvmExecution::evaluate(math::ScalarDiffusionFlux operation, ScalarField& target) {
    requireDistinct(&operation.coefficient, &target);
    Implementation& state = *m_implementation;
    operation.options.requireComplete("scalar diffusion flux");
    requireCellField(operation.field, *state.mesh, "diffusion-flux field");
    requireFaceField(target, *state.mesh, "diffusion-flux target");
    if (&operation.coefficient.mesh() != state.mesh ||
        (operation.coefficient.location() != FieldLocation::Cell &&
         operation.coefficient.location() != FieldLocation::Face)) {
        throw std::invalid_argument(
            "diffusion-flux coefficient must be a cell or face field on the run mesh");
    }
    if (operation.gradient != nullptr)
        requireCellField(*operation.gradient, *state.mesh, "diffusion-flux gradient");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    state.synchronize(const_cast<ScalarField&>(operation.coefficient));
    const VectorField* reconstructed_gradient = operation.gradient;
    if (reconstructed_gradient == nullptr) {
        gradient(operation.field, state.gradient_workspace,
                 operation.options.gradient.value());
        reconstructed_gradient = &state.gradient_workspace;
    }
    state.synchronize(const_cast<VectorField&>(*reconstructed_gradient));

    const ScalarField* face_coefficient = &operation.coefficient;
    if (operation.coefficient.location() == FieldLocation::Cell) {
        interpolate(
            operation.coefficient, state.face_coefficient_workspace,
            operation.options.coefficientInterpolation.value_or(operation.options.interpolation.value()),
            operation.options.coefficientGradient.value_or(operation.options.gradient.value()));
        state.synchronize(state.face_coefficient_workspace);
        face_coefficient = &state.face_coefficient_workspace;
    }
    diffusionFlux(
        *face_coefficient, operation.field, *reconstructed_gradient,
        target, operation.options.diffusion.value());
    state.synchronize(target);
}

void FvmExecution::add(math::FaceFlux operation, ScalarField& target, math::FaceRegion region) {
    requireFaceField(target, *m_implementation->mesh, "flux increment target");
    // Cell 输入的入口/出口判定使用目标中的当前通量，而不是上一次工作区残值。
    if (operation.velocity.location() == FieldLocation::Cell) {
        m_implementation->synchronize(target);
        m_implementation->face_flux_workspace.assign(target);
    }
    evaluate(operation, m_implementation->face_flux_workspace);
    m_implementation->addFaceIncrement(target, 1.0, region);
}

void FvmExecution::subtract(math::ScalarDiffusionFlux operation, ScalarField& target,
                            math::FaceRegion region) {
    requireFaceField(target, *m_implementation->mesh, "flux increment target");
    requireDistinct(&operation.coefficient, &target);
    evaluate(operation, m_implementation->face_flux_workspace);
    m_implementation->addFaceIncrement(target, -1.0, region);
}

void FvmExecution::evaluate(math::NormalGradient operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    operation.options.requireComplete("normal gradient");
    requireCellField(operation.field, *state.mesh, "normal-gradient input");
    requireFaceField(result, *state.mesh, "normal-gradient result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    const VectorField* scalar_gradient = operation.gradient;
    if (scalar_gradient) {
        requireCellField(*scalar_gradient, *state.mesh, "normal-gradient supplied gradient");
        state.synchronize(const_cast<VectorField&>(*scalar_gradient));
    } else {
        gradient(operation.field, state.gradient_workspace,
                 operation.options.gradient.value());
        state.synchronize(state.gradient_workspace);
        scalar_gradient = &state.gradient_workspace;
    }
    for (Index face : detail::meshData(*state.mesh).owned_faces) {
        detail::fieldData(result)[face] = integratedNormalGradient(operation.field, *scalar_gradient, face,
            operation.options.diffusion.value()) / state.mesh->faceArea(face);
    }
    state.synchronize(result);
}

}  // detail 命名空间

}  // babelsim 命名空间
