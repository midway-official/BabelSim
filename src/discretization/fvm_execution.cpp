#include "internal/compute_backend.h"
#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "internal/fvm_execution.h"

#include "babelsim/discrete_equation.h"
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

int residualSign(bool lhs, int term_sign) {
    return (lhs ? 1 : -1) * term_sign;
}

void requireDistinct(const void* input, const void* result) {
    if (input == result)
        throw std::invalid_argument("explicit operator input and result must not alias");
}

}  // 匿名命名空间

struct FvmExecution::Implementation {
    template <typename T>
    struct History {
        explicit History(const Field<T>& field)
            : source(&field), previous(field), older(field)
        {}

        const Field<T>* source;
        Field<T> previous;
        Field<T> older;
        bool has_older = false;

        void advance(const Field<T>& current) {
            older.assign(previous);
            previous.assign(current);
            has_older = true;
        }
    };

    using ScalarHistory = History<double>;
    using VectorHistory = History<Vec3>;

    Implementation(const Mesh& mesh_value, const Methods& methods_value,
                   std::unique_ptr<ComputeBackend> backend_value,
                   double initial_delta_t)
        : mesh(&mesh_value),
          methods(methods_value),
          backend(std::move(backend_value)),
          gradient_workspace(mesh_value, FieldLocation::Cell, "grad"),
          face_coefficient_workspace(mesh_value, FieldLocation::Face, "faceCoefficient"),
          face_flux_workspace(mesh_value, FieldLocation::Face, "faceFlux"),
          delta_t(initial_delta_t)
    {
        mesh->validate();
        if (!backend) throw std::invalid_argument("FVM execution requires a compute backend");
    }

    ScalarHistory& history(const ScalarField& field) {
        for (ScalarHistory& value : scalar_histories) {
            if (value.source == &field) return value;
        }
        scalar_histories.emplace_back(field);
        return scalar_histories.back();
    }

    VectorHistory& history(const VectorField& field) {
        for (VectorHistory& value : vector_histories) {
            if (value.source == &field) return value;
        }
        vector_histories.emplace_back(field);
        return vector_histories.back();
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
    Methods methods;
    std::unique_ptr<ComputeBackend> backend;
    std::vector<ScalarHistory> scalar_histories;
    std::vector<VectorHistory> vector_histories;
    VectorField gradient_workspace;
    ScalarField face_coefficient_workspace;
    ScalarField face_flux_workspace;
    double delta_t;
    bool has_time_step = false;
};


FvmExecution::FvmExecution(const Mesh& mesh, const Methods& methods,
                           std::unique_ptr<ComputeBackend> backend, double delta_t)
    : m_implementation(std::make_unique<Implementation>(
          mesh, methods, std::move(backend), delta_t)) {}
FvmExecution::~FvmExecution() = default;

void FvmExecution::beginStep(double delta_t) {
    Implementation& state = *m_implementation;
    if (state.has_time_step) {
        for (auto& history : state.scalar_histories) history.advance(*history.source);
        for (auto& history : state.vector_histories) history.advance(*history.source);
    }
    state.delta_t = delta_t;
    state.has_time_step = true;
}

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
        for (Index face : detail::meshData(*state.mesh).cell_faces[static_cast<std::size_t>(cell)]) {
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

namespace {

void addScalarSource(
    ScalarDiscreteEquation& equation, const Mesh& mesh, const ScalarEquationTerm& term, int canonical)
{
    const double scale = -static_cast<double>(canonical) * term.coefficient;
    if (term.field == nullptr) {
        for (Index cell : detail::meshData(mesh).owned_cells) {
            equation.source[static_cast<std::size_t>(cell)] +=
                scale * detail::meshData(mesh).cell_volumes[static_cast<std::size_t>(cell)];
        }
        return;
    }
    requireCellField(*term.field, mesh, "scalar source");
    for (Index cell : detail::meshData(mesh).owned_cells) {
        equation.source[static_cast<std::size_t>(cell)] +=
            scale * detail::fieldData((*term.field))[cell] * detail::meshData(mesh).cell_volumes[static_cast<std::size_t>(cell)];
    }
}

void addVectorSource(
    VectorDiscreteEquation& equation, const Mesh& mesh, const VectorEquationTerm& term, int canonical)
{
    const double scale = -static_cast<double>(canonical) * term.coefficient;
    if (term.vector_field != nullptr) requireCellField(*term.vector_field, mesh, "vector source");
    for (Index cell : detail::meshData(mesh).owned_cells) {
        equation.source[static_cast<std::size_t>(cell)] +=
            scale * (term.vector_field != nullptr ? detail::fieldData((*term.vector_field))[cell] : term.vector_source) *
                detail::meshData(mesh).cell_volumes[static_cast<std::size_t>(cell)];
    }
}

}  // 匿名命名空间

SolveResult FvmExecution::solve(
    const ScalarEquationDefinition& expression,
    EquationControl equation_control)
{
    Implementation& state = *m_implementation;
    if (!(equation_control.relaxation > 0.0 && equation_control.relaxation <= 1.0) ||
        !std::isfinite(equation_control.reference_value))
        throw std::invalid_argument("invalid equation relaxation or reference value");
    const ScalarField* unknown_pointer = nullptr;
    const auto inspect_unknown = [&unknown_pointer](const std::vector<ScalarEquationTerm>& terms) {
        for (const ScalarEquationTerm& term : terms) {
            if (term.kind == EquationTermKind::Source || term.field == nullptr) continue;
            if (unknown_pointer != nullptr && unknown_pointer != term.field) {
                throw std::invalid_argument("a scalar equation must have one transported field");
            }
            unknown_pointer = term.field;
        }
    };
    inspect_unknown(expression.lhs.m_terms);
    inspect_unknown(expression.rhs.m_terms);
    if (unknown_pointer == nullptr) {
        throw std::invalid_argument("scalar equation has no unknown field");
    }
    ScalarField& unknown = const_cast<ScalarField&>(*unknown_pointer);
    requireCellField(unknown, *state.mesh, "scalar unknown");
    ScalarDiscreteEquation equation(*state.mesh);

    const auto add = [&](const std::vector<ScalarEquationTerm>& terms, bool lhs) {
        for (const ScalarEquationTerm& term : terms) {
            const int canonical = residualSign(lhs, term.sign);
            switch (term.kind) {
                case EquationTermKind::TimeDerivative: {
                    if (canonical != 1 || term.field != &unknown) {
                        throw std::invalid_argument("ddt must be a positive left-hand-side term");
                    }
                    Implementation::ScalarHistory& history = state.history(unknown);
                    // BDF2 首步还没有两个历史层，用 Euler 启动。
                    const TimeMethod method = state.methods.time == TimeMethod::BDF2 &&
                        !history.has_older ? TimeMethod::Euler : state.methods.time;
                    if (term.coefficient_field != nullptr) {
                        requireCellField(
                            *term.coefficient_field, *state.mesh, "time coefficient");
                        state.synchronize(const_cast<ScalarField&>(*term.coefficient_field));
                        addTimeDerivative(
                            equation, history.previous, state.delta_t, *term.coefficient_field,
                            method,
                            history.has_older ? &history.older : nullptr);
                    } else {
                        addTimeDerivative(
                            equation, history.previous, state.delta_t, term.coefficient,
                            method,
                            history.has_older ? &history.older : nullptr);
                    }
                    break;
                }
                case EquationTermKind::Convection:
                    if (canonical != 1 || term.field != &unknown || term.flux == nullptr) {
                        throw std::invalid_argument("implicit convection must be a positive left-hand-side term");
                    }
                    requireFaceField(*term.flux, *state.mesh, "convection flux");
                    state.synchronize(const_cast<ScalarField&>(*term.flux));
                    state.synchronize(unknown);
                    addConvection(
                        equation, *term.flux, unknown,
                        state.methods.convectionFor(unknown.name()),
                        state.methods.interpolationFor(unknown.name()),
                        state.methods.gradientFor(unknown.name()),
                        term.coefficient);
                    break;
                case EquationTermKind::Laplacian:
                    if (canonical != -1 || term.field != &unknown) {
                        throw std::invalid_argument("laplacian must appear on the right-hand side or negated on the left");
                    }
                    state.synchronize(unknown);
                    if (term.coefficient_field != nullptr) {
                        const ScalarField& coefficient = *term.coefficient_field;
                        if (&coefficient.mesh() != state.mesh ||
                            (coefficient.location() != FieldLocation::Cell &&
                             coefficient.location() != FieldLocation::Face)) {
                            throw std::invalid_argument("diffusivity must be a cell or face field");
                        }
                        state.synchronize(const_cast<ScalarField&>(coefficient));
                        const ScalarField* face_coefficient = &coefficient;
                        if (coefficient.location() == FieldLocation::Cell) {
                            interpolate(
                                coefficient, state.face_coefficient_workspace,
                                state.methods.interpolationFor(unknown.name()),
                                state.methods.gradientFor(unknown.name()));
                            state.synchronize(state.face_coefficient_workspace);
                            face_coefficient = &state.face_coefficient_workspace;
                        }
                        addDiffusion(
                            equation, *face_coefficient, unknown,
                            state.methods.gradientFor(unknown.name()),
                            state.methods.diffusionFor(unknown.name()));
                    } else {
                        addDiffusion(
                            equation, term.coefficient, unknown,
                            state.methods.gradientFor(unknown.name()),
                            state.methods.diffusionFor(unknown.name()));
                    }
                    break;
                case EquationTermKind::Source:
                    addScalarSource(equation, *state.mesh, term, canonical);
                    break;
                case EquationTermKind::Gradient:
                    throw std::invalid_argument("gradient is not a scalar implicit term");
            }
        }
    };
    add(expression.lhs.m_terms, true);
    add(expression.rhs.m_terms, false);

    if (equation_control.relaxation != 1.0) {
        for (Index cell : detail::meshData(*state.mesh).owned_cells) {
            const std::size_t index = static_cast<std::size_t>(cell);
            const double added = equation.diagonal[index] * (1.0 / equation_control.relaxation - 1.0);
            equation.diagonal[index] += added;
            equation.source[index] += added * detail::fieldData(unknown)[cell];
        }
    }
    if (equation_control.fix_reference) {
        for (Index cell : detail::meshData(*state.mesh).owned_cells) {
            if (detail::globalCellId(*state.mesh, cell) != 0) continue;
            const std::size_t index = static_cast<std::size_t>(cell);
            if (!(equation.diagonal[index] > 0.0)) {
                throw std::runtime_error("reference equation diagonal is invalid");
            }
            equation.source[index] += equation.diagonal[index] * equation_control.reference_value;
            equation.diagonal[index] += equation.diagonal[index];
            break;
        }
    }

    return state.backend->solve(equation, unknown);
}

std::array<SolveResult, 3> FvmExecution::solve(
    const VectorEquationDefinition& expression,
    VectorEquationControl equation_control)
{
    Implementation& state = *m_implementation;
    if (!(equation_control.relaxation > 0.0 &&
          equation_control.relaxation <= 1.0)) {
        throw std::invalid_argument("vector equation relaxation must be in (0, 1]");
    }
    if (equation_control.mobility != nullptr) {
        requireCellField(*equation_control.mobility, *state.mesh, "equation mobility");
    }
    const VectorField* unknown_pointer = nullptr;
    const auto inspect_unknown = [&unknown_pointer, &equation_control](const std::vector<VectorEquationTerm>& terms) {
        for (const VectorEquationTerm& term : terms) {
            if (equation_control.mobility != nullptr &&
                (term.scalar_field == equation_control.mobility ||
                 term.coefficient_field == equation_control.mobility ||
                 term.flux == equation_control.mobility))
                throw std::invalid_argument("equation response must not alias an input field");
            if (term.kind == EquationTermKind::Source || term.vector_field == nullptr) continue;
            if (unknown_pointer != nullptr && unknown_pointer != term.vector_field) {
                throw std::invalid_argument("a vector equation must have one transported field");
            }
            unknown_pointer = term.vector_field;
        }
    };
    inspect_unknown(expression.lhs.m_terms);
    inspect_unknown(expression.rhs.m_terms);
    if (unknown_pointer == nullptr) {
        throw std::invalid_argument("vector equation has no unknown field");
    }
    VectorField& unknown = const_cast<VectorField&>(*unknown_pointer);
    requireCellField(unknown, *state.mesh, "vector unknown");
    VectorDiscreteEquation equation(*state.mesh);

    const auto add = [&](const std::vector<VectorEquationTerm>& terms, bool lhs) {
        for (const VectorEquationTerm& term : terms) {
            const int canonical = residualSign(lhs, term.sign);
            switch (term.kind) {
                case EquationTermKind::TimeDerivative: {
                    if (canonical != 1 || term.vector_field != &unknown) {
                        throw std::invalid_argument("ddt must be a positive left-hand-side term");
                    }
                    Implementation::VectorHistory& history = state.history(unknown);
                    // BDF2 首步还没有两个历史层，用 Euler 启动。
                    const TimeMethod method = state.methods.time == TimeMethod::BDF2 &&
                        !history.has_older ? TimeMethod::Euler : state.methods.time;
                    if (term.coefficient_field != nullptr) {
                        requireCellField(
                            *term.coefficient_field, *state.mesh, "time coefficient");
                        state.synchronize(const_cast<ScalarField&>(*term.coefficient_field));
                        addTimeDerivative(
                            equation, history.previous, state.delta_t, *term.coefficient_field,
                            method,
                            history.has_older ? &history.older : nullptr);
                    } else {
                        addTimeDerivative(
                            equation, history.previous, state.delta_t, term.coefficient,
                            method,
                            history.has_older ? &history.older : nullptr);
                    }
                    break;
                }
                case EquationTermKind::Convection:
                    if (canonical != 1 || term.vector_field != &unknown || term.flux == nullptr) {
                        throw std::invalid_argument("implicit convection must be a positive left-hand-side term");
                    }
                    requireFaceField(*term.flux, *state.mesh, "convection flux");
                    state.synchronize(const_cast<ScalarField&>(*term.flux));
                    state.synchronize(unknown);
                    addConvection(
                        equation, *term.flux, unknown,
                        state.methods.convectionFor(unknown.name()),
                        state.methods.interpolationFor(unknown.name()),
                        state.methods.gradientFor(unknown.name()),
                        term.coefficient);
                    break;
                case EquationTermKind::Laplacian:
                    if (canonical != -1 || term.vector_field != &unknown) {
                        throw std::invalid_argument("laplacian must appear on the right-hand side or negated on the left");
                    }
                    state.synchronize(unknown);
                    if (term.coefficient_field != nullptr) {
                        const ScalarField& coefficient = *term.coefficient_field;
                        if (&coefficient.mesh() != state.mesh ||
                            (coefficient.location() != FieldLocation::Cell &&
                             coefficient.location() != FieldLocation::Face)) {
                            throw std::invalid_argument("diffusivity must be a cell or face field");
                        }
                        state.synchronize(const_cast<ScalarField&>(coefficient));
                        const ScalarField* face_coefficient = &coefficient;
                        if (coefficient.location() == FieldLocation::Cell) {
                            interpolate(
                                coefficient, state.face_coefficient_workspace,
                                state.methods.interpolationFor(unknown.name()),
                                state.methods.gradientFor(unknown.name()));
                            state.synchronize(state.face_coefficient_workspace);
                            face_coefficient = &state.face_coefficient_workspace;
                        }
                        addDiffusion(
                            equation, *face_coefficient, unknown,
                            state.methods.gradientFor(unknown.name()),
                            state.methods.diffusionFor(unknown.name()));
                    } else {
                        addDiffusion(
                            equation, term.coefficient, unknown,
                            state.methods.gradientFor(unknown.name()),
                            state.methods.diffusionFor(unknown.name()));
                    }
                    break;
                case EquationTermKind::Gradient:
                    if (term.scalar_field == nullptr) {
                        throw std::invalid_argument("gradient term has no scalar field");
                    }
                    requireCellField(*term.scalar_field, *state.mesh, "gradient field");
                    state.synchronize(const_cast<ScalarField&>(*term.scalar_field));
                    gradient(
                        *term.scalar_field, state.gradient_workspace,
                        state.methods.gradientFor(term.scalar_field->name()));
                    state.synchronize(state.gradient_workspace);
                    for (Index cell : detail::meshData(*state.mesh).owned_cells) {
                        equation.source[static_cast<std::size_t>(cell)] -=
                            static_cast<double>(canonical) *
                            detail::meshData(*state.mesh).cell_volumes[static_cast<std::size_t>(cell)] *
                            detail::fieldData(state.gradient_workspace)[cell];
                    }
                    break;
                case EquationTermKind::Source:
                    addVectorSource(equation, *state.mesh, term, canonical);
                    break;
            }
        }
    };
    add(expression.lhs.m_terms, true);
    add(expression.rhs.m_terms, false);

    if (equation_control.relaxation != 1.0) {
        for (Index face : detail::meshData(*state.mesh).owned_faces) {
            equation.upper[static_cast<std::size_t>(face)] *= equation_control.relaxation;
            equation.lower[static_cast<std::size_t>(face)] *= equation_control.relaxation;
        }
    }
    for (Index cell : detail::meshData(*state.mesh).owned_cells) {
        const std::size_t index = static_cast<std::size_t>(cell);
        if (!(equation.diagonal[index] > 0.0) || !std::isfinite(equation.diagonal[index])) {
            throw std::runtime_error("vector equation diagonal is not positive and finite");
        }
        if (equation_control.relaxation != 1.0) {
            equation.source[index] = equation_control.relaxation * equation.source[index] +
                (1.0 - equation_control.relaxation) *
                equation.diagonal[index] * detail::fieldData(unknown)[cell];
        }
        if (equation_control.mobility != nullptr) {
            detail::fieldData((*equation_control.mobility))[cell] =
                detail::meshData(*state.mesh).cell_volumes[index] / equation.diagonal[index];
        }
    }
    if (equation_control.mobility != nullptr) {
        state.synchronize(*equation_control.mobility);
    }

    return state.backend->solve(equation, unknown);
}

void FvmExecution::evaluate(math::ScalarGradient operation, VectorField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "gradient input");
    requireCellField(result, *state.mesh, "gradient result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    gradient(operation.field, result, state.methods.gradientFor(operation.field.name()));
    state.synchronize(result);
}

void FvmExecution::evaluate(math::VectorGradient operation, TensorField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "gradient input");
    requireCellField(result, *state.mesh, "gradient result");
    state.synchronize(const_cast<VectorField&>(operation.field));
    gradient(operation.field, result, state.methods.gradientFor(operation.field.name()));
    state.synchronize(result);
}

void FvmExecution::evaluate(math::FaceFlux operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireFaceField(result, *state.mesh, "flux result");
    state.synchronize(const_cast<VectorField&>(operation.velocity));
    flux(operation.velocity, result,
         state.methods.interpolationFor(operation.velocity.name()),
         state.methods.gradientFor(operation.velocity.name()));
    state.synchronize(result);
}

void FvmExecution::evaluate(math::FaceDivergence operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireFaceField(operation.flux, *state.mesh, "divergence input");
    requireCellField(result, *state.mesh, "divergence result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    divergence(operation.flux, result);
    state.synchronize(result);
}

void FvmExecution::evaluate(math::VectorDivergence operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "divergence input");
    requireCellField(result, *state.mesh, "divergence result");
    state.synchronize(const_cast<VectorField&>(operation.field));
    divergence(
        operation.field, result, state.methods.interpolationFor(operation.field.name()),
        state.methods.gradientFor(operation.field.name()));
    state.synchronize(result);
}

void FvmExecution::evaluate(math::ScalarConvection operation, ScalarField& result) {
    requireDistinct(&operation.field, &result);
    Implementation& state = *m_implementation;
    requireFaceField(operation.flux, *state.mesh, "convection flux");
    requireCellField(operation.field, *state.mesh, "convection field");
    requireCellField(result, *state.mesh, "convection result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    state.synchronize(const_cast<ScalarField&>(operation.field));
    convection(
        operation.flux, operation.field, result,
        state.methods.convectionFor(operation.field.name()),
        state.methods.interpolationFor(operation.field.name()),
        state.methods.gradientFor(operation.field.name()));
    state.synchronize(result);
}

void FvmExecution::evaluate(math::VectorConvection operation, VectorField& result) {
    requireDistinct(&operation.field, &result);
    Implementation& state = *m_implementation;
    requireFaceField(operation.flux, *state.mesh, "convection flux");
    requireCellField(operation.field, *state.mesh, "convection field");
    requireCellField(result, *state.mesh, "convection result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    state.synchronize(const_cast<VectorField&>(operation.field));
    convection(
        operation.flux, operation.field, result,
        state.methods.convectionFor(operation.field.name()),
        state.methods.interpolationFor(operation.field.name()),
        state.methods.gradientFor(operation.field.name()));
    state.synchronize(result);
}

void FvmExecution::evaluate(math::ScalarInterpolation operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "interpolation input");
    requireFaceField(result, *state.mesh, "interpolation result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    interpolate(
        operation.field, result, state.methods.interpolationFor(operation.field.name()),
        state.methods.gradientFor(operation.field.name()));
    state.synchronize(result);
}

void FvmExecution::evaluate(math::VectorInterpolation operation, VectorField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "interpolation input");
    if (&result.mesh() != state.mesh || result.location() != FieldLocation::Face) {
        throw std::invalid_argument("vector interpolation result must be a face field");
    }
    state.synchronize(const_cast<VectorField&>(operation.field));
    interpolate(
        operation.field, result, state.methods.interpolationFor(operation.field.name()),
        state.methods.gradientFor(operation.field.name()));
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
    requireCellField(operation.field, *state.mesh, "laplacian input");
    requireCellField(result, *state.mesh, "laplacian result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    if (operation.coefficient_field == nullptr) {
        laplacian(
            operation.coefficient, operation.field, result,
            state.methods.gradientFor(operation.field.name()),
            state.methods.diffusionFor(operation.field.name()));
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
                state.methods.interpolationFor(operation.field.name()),
                state.methods.gradientFor(operation.field.name()));
            state.synchronize(state.face_coefficient_workspace);
            laplacian(
                state.face_coefficient_workspace, operation.field, result,
                state.methods.gradientFor(operation.field.name()),
                state.methods.diffusionFor(operation.field.name()));
        } else {
            laplacian(
                coefficient, operation.field, result,
                state.methods.gradientFor(operation.field.name()),
                state.methods.diffusionFor(operation.field.name()));
        }
    }
    state.synchronize(result);
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
    gradient(
        operation.field, state.gradient_workspace,
        state.methods.gradientFor(operation.field.name()));
    state.synchronize(state.gradient_workspace);
    target.addProduct(-1.0, coefficient, state.gradient_workspace);
    state.synchronize(target);
}

void FvmExecution::evaluate(math::ScalarDiffusionFlux operation, ScalarField& target) {
    requireDistinct(&operation.coefficient, &target);
    Implementation& state = *m_implementation;
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
                 state.methods.gradientFor(operation.field.name()));
        reconstructed_gradient = &state.gradient_workspace;
    }
    state.synchronize(const_cast<VectorField&>(*reconstructed_gradient));

    const ScalarField* face_coefficient = &operation.coefficient;
    if (operation.coefficient.location() == FieldLocation::Cell) {
        interpolate(
            operation.coefficient, state.face_coefficient_workspace,
            state.methods.interpolationFor(operation.field.name()),
            state.methods.gradientFor(operation.field.name()));
        state.synchronize(state.face_coefficient_workspace);
        face_coefficient = &state.face_coefficient_workspace;
    }
    diffusionFlux(
        *face_coefficient, operation.field, *reconstructed_gradient,
        target, state.methods.diffusionFor(operation.field.name()));
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
    requireCellField(operation.field, *state.mesh, "normal-gradient input");
    requireFaceField(result, *state.mesh, "normal-gradient result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    gradient(operation.field, state.gradient_workspace,
             state.methods.gradientFor(operation.field.name()));
    state.synchronize(state.gradient_workspace);
    for (Index face : detail::meshData(*state.mesh).owned_faces) {
        detail::fieldData(result)[face] = integratedNormalGradient(operation.field, state.gradient_workspace, face,
            state.methods.diffusionFor(operation.field.name())) / state.mesh->faceArea(face);
    }
    state.synchronize(result);
}

}  // detail 命名空间

}  // babelsim 命名空间
