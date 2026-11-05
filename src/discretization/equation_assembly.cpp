#include "babelsim/operators.h"
#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "operator_support.h"

#include <cmath>
#include <optional>
#include <stdexcept>

namespace babelsim {
using namespace detail::operators;
namespace {

template <typename T>
void requireEquation(const DiscreteEquation<T>& equation, const Mesh& mesh) {
    if (equation.mesh != &mesh) {
        throw std::invalid_argument("equation and field use different meshes");
    }
}

template <typename T>
void addConvectionImpl(DiscreteEquation<T>& equation, const ScalarField& face_flux,
                       const Field<T>& transported, ConvectionMethod method,
                       InterpolationMethod interpolation_method, GradientMethod gradient_method,
                       double flux_scale) {
    const Mesh& mesh = transported.mesh();
    requireEquation(equation, mesh);
    requireField(face_flux, mesh, FieldLocation::Face, "face flux");
    requireField(transported, mesh, FieldLocation::Cell, "transported");
    const_cast<Field<T>&>(transported).setBoundaryFlux(face_flux);

    if (interpolation_method != InterpolationMethod::Linear &&
        interpolation_method != InterpolationMethod::Corrected) {
        throw std::invalid_argument("unsupported convection interpolation method");
    }
    if (!std::isfinite(flux_scale)) {
        throw std::invalid_argument("convection flux scale must be finite");
    }
    using GradientField = typename GradientFieldType<T>::Type;
    std::optional<GradientField> transported_gradient;
    if (method == ConvectionMethod::LinearUpwind ||
        (method == ConvectionMethod::Central &&
         interpolation_method == InterpolationMethod::Corrected)) {
        transported_gradient.emplace(mesh, FieldLocation::Cell, "grad(" + transported.name() + ')');
        gradient(transported, *transported_gradient, gradient_method);
    }

    for (Index face : detail::meshData(mesh).owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double integrated_flux = flux_scale * detail::fieldData(face_flux)[face];
        if (neighbour != invalid_index) {
            if (method == ConvectionMethod::Upwind || method == ConvectionMethod::LinearUpwind) {
                equation.diagonal[static_cast<std::size_t>(owner)] +=
                    std::max(integrated_flux, 0.0);
                equation.upper[f] += std::min(integrated_flux, 0.0);
                equation.diagonal[static_cast<std::size_t>(neighbour)] +=
                    std::max(-integrated_flux, 0.0);
                equation.lower[f] += std::min(-integrated_flux, 0.0);
                if (method == ConvectionMethod::LinearUpwind) {
                    const Index upwind = integrated_flux >= 0.0 ? owner : neighbour;
                    const T correction =
                        linearUpwindCorrection(transported, *transported_gradient, upwind, face);
                    equation.source[static_cast<std::size_t>(owner)] -=
                        integrated_flux * correction;
                    equation.source[static_cast<std::size_t>(neighbour)] +=
                        integrated_flux * correction;
                }
            } else if (method == ConvectionMethod::Central) {
                const double weight = ownerWeight(mesh, face);
                equation.diagonal[static_cast<std::size_t>(owner)] += integrated_flux * weight;
                equation.upper[f] += integrated_flux * (1.0 - weight);
                equation.lower[f] -= integrated_flux * weight;
                equation.diagonal[static_cast<std::size_t>(neighbour)] -=
                    integrated_flux * (1.0 - weight);
                if (transported_gradient) {
                    const T correction =
                        correctedInternalFaceValue(transported, *transported_gradient, face) -
                        internalFaceValue(transported, face);
                    equation.source[static_cast<std::size_t>(owner)] -=
                        integrated_flux * correction;
                    equation.source[static_cast<std::size_t>(neighbour)] +=
                        integrated_flux * correction;
                }
            } else {
                throw std::invalid_argument("unsupported convection method");
            }
            continue;
        }

        const auto& condition = detail::FieldAccess::condition(transported, face, integrated_flux);
        if (method == ConvectionMethod::Upwind || method == ConvectionMethod::LinearUpwind) {
            if (integrated_flux >= 0.0) {
                equation.diagonal[static_cast<std::size_t>(owner)] += integrated_flux;
                if (method == ConvectionMethod::LinearUpwind) {
                    equation.source[static_cast<std::size_t>(owner)] -=
                        integrated_flux *
                        linearUpwindCorrection(transported, *transported_gradient, owner, face);
                }
            } else {
                equation.source[static_cast<std::size_t>(owner)] -=
                    integrated_flux * boundaryFaceValue(transported, face, integrated_flux);
            }
            continue;
        }
        if (method != ConvectionMethod::Central) {
            throw std::invalid_argument("unsupported convection method");
        }
        switch (condition.type) {
        case BoundaryType::FixedValue:
            equation.source[static_cast<std::size_t>(owner)] -= integrated_flux * condition.value;
            break;
        case BoundaryType::FixedGradient: {
            const double distance = boundaryNormalDistance(mesh, face);
            equation.diagonal[static_cast<std::size_t>(owner)] += integrated_flux;
            equation.source[static_cast<std::size_t>(owner)] -=
                integrated_flux * distance * condition.value;
            break;
        }
        case BoundaryType::ZeroGradient:
        case BoundaryType::Symmetry:
            equation.diagonal[static_cast<std::size_t>(owner)] += integrated_flux;
            break;
        case BoundaryType::InletOutlet:
            if (integrated_flux >= 0.0) {
                equation.diagonal[static_cast<std::size_t>(owner)] += integrated_flux;
            } else {
                equation.source[static_cast<std::size_t>(owner)] -=
                    integrated_flux * condition.value;
            }
            break;
        }
        if (transported_gradient) {
            const T correction = reconstructedBoundaryValue(transported, *transported_gradient,
                                                            face, integrated_flux) -
                                 boundaryFaceValue(transported, face, integrated_flux);
            equation.source[static_cast<std::size_t>(owner)] -= integrated_flux * correction;
        }
    }
}

template <typename T>
void requireTimeHistory(const Mesh& mesh, TimeMethod method, const Field<T>* older) {
    if (method == TimeMethod::BDF2) {
        if (older == nullptr) {
            throw std::invalid_argument("BDF2 requires two previous fields");
        }
        requireField(*older, mesh, FieldLocation::Cell, "older");
    } else if (method != TimeMethod::Euler) {
        throw std::invalid_argument("unsupported time method");
    }
}

template <typename T>
void addCellTimeDerivative(DiscreteEquation<T>& equation, Index cell, double coefficient,
                           const Field<T>& previous, TimeMethod method, const Field<T>* older) {
    const auto c = static_cast<std::size_t>(cell);
    // coefficient is capacity * volume / dt. Keep the same arithmetic order
    // for a constant capacity and a cell field, including the BDF2 RHS.
    if (method == TimeMethod::Euler) {
        equation.diagonal[c] += coefficient;
        equation.source[c] += coefficient * detail::fieldData(previous)[cell];
    } else {
        equation.diagonal[c] += 1.5 * coefficient;
        equation.source[c] += coefficient * (2.0 * detail::fieldData(previous)[cell] -
                                             0.5 * detail::fieldData(*older)[cell]);
    }
}

template <typename T>
void addTimeDerivativeImpl(DiscreteEquation<T>& equation, const Field<T>& previous, double dt,
                           double density, TimeMethod method, const Field<T>* older) {
    const Mesh& mesh = previous.mesh();
    requireEquation(equation, mesh);
    requireField(previous, mesh, FieldLocation::Cell, "previous");
    if (method == TimeMethod::Steady) {
        return;
    }
    if (!(dt > 0.0) || !(density > 0.0) || !std::isfinite(dt) || !std::isfinite(density)) {
        throw std::invalid_argument("time step and density must be positive and finite");
    }
    requireTimeHistory(mesh, method, older);
    for (Index cell : detail::meshData(mesh).owned_cells) {
        const auto c = static_cast<std::size_t>(cell);
        const double coefficient = density * detail::meshData(mesh).cell_volumes[c] / dt;
        addCellTimeDerivative(equation, cell, coefficient, previous, method, older);
    }
}

template <typename T>
void addFieldTimeDerivativeImpl(DiscreteEquation<T>& equation, const Field<T>& previous, double dt,
                                const ScalarField& capacity, TimeMethod method,
                                const Field<T>* older) {
    const Mesh& mesh = previous.mesh();
    requireEquation(equation, mesh);
    requireField(previous, mesh, FieldLocation::Cell, "previous");
    requireField(capacity, mesh, FieldLocation::Cell, "time coefficient");
    if (method == TimeMethod::Steady)
        return;
    if (!(dt > 0.0) || !std::isfinite(dt)) {
        throw std::invalid_argument("time step must be positive and finite");
    }
    requireTimeHistory(mesh, method, older);
    for (Index cell : detail::meshData(mesh).owned_cells) {
        const std::size_t index = static_cast<std::size_t>(cell);
        if (!(detail::fieldData(capacity)[cell] > 0.0) ||
            !std::isfinite(detail::fieldData(capacity)[cell])) {
            throw std::invalid_argument("time coefficient must be positive and finite");
        }
        const double coefficient =
            detail::fieldData(capacity)[cell] * detail::meshData(mesh).cell_volumes[index] / dt;
        addCellTimeDerivative(equation, cell, coefficient, previous, method, older);
    }
}

} // namespace

void addConvection(ScalarDiscreteEquation& equation, const ScalarField& face_flux,
                   const ScalarField& transported, ConvectionMethod method,
                   InterpolationMethod interpolation_method, GradientMethod gradient_method,
                   double flux_scale) {
    addConvectionImpl(equation, face_flux, transported, method, interpolation_method,
                      gradient_method, flux_scale);
}

void addConvection(VectorDiscreteEquation& equation, const ScalarField& face_flux,
                   const VectorField& transported, ConvectionMethod method,
                   InterpolationMethod interpolation_method, GradientMethod gradient_method,
                   double flux_scale) {
    addConvectionImpl(equation, face_flux, transported, method, interpolation_method,
                      gradient_method, flux_scale);
}

namespace {

struct FaceDiffusivity {
    double constant = 0.0;
    const ScalarField* field = nullptr;

    double operator[](Index face) const {
        return field == nullptr ? constant : detail::fieldData((*field))[face];
    }
};

void addScalarDiffusion(ScalarDiscreteEquation& equation, FaceDiffusivity diffusivity,
                        const ScalarField& scalar, GradientMethod gradient_method,
                        DiffusionMethod diffusion_method,
                        const VectorField* reconstructed_gradient = nullptr,
                        bool assemble_matrix = true) {
    const Mesh& mesh = scalar.mesh();
    requireEquation(equation, mesh);
    requireField(scalar, mesh, FieldLocation::Cell, "scalar");

    std::optional<VectorField> gradient_workspace;
    const VectorField* cell_gradient = reconstructed_gradient;
    if (cell_gradient)
        requireField(*cell_gradient, mesh, FieldLocation::Cell, "diffusion gradient");
    if (correctedDiffusion(diffusion_method)) {
        if (!cell_gradient) {
            gradient_workspace.emplace(mesh, FieldLocation::Cell, "grad(" + scalar.name() + ')');
            gradient(scalar, *gradient_workspace, gradient_method);
            cell_gradient = &*gradient_workspace;
        }
    } else if (diffusion_method != DiffusionMethod::Orthogonal) {
        throw std::invalid_argument("unsupported diffusion method");
    }

    for (Index face : detail::meshData(mesh).owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double gamma = diffusivity[face];
        if (gamma < 0.0 || !std::isfinite(gamma)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        if (gamma == 0.0) {
            continue;
        }
        const double coefficient = gamma * detail::meshData(mesh).face_orthogonal_coefficients[f];
        if (neighbour != invalid_index) {
            if (assemble_matrix) {
                equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
                equation.diagonal[static_cast<std::size_t>(neighbour)] += coefficient;
                equation.upper[f] -= coefficient;
                equation.lower[f] -= coefficient;
            }
            if (correctedDiffusion(diffusion_method)) {
                const double orthogonal =
                    detail::meshData(mesh).face_orthogonal_coefficients[f] *
                    (detail::fieldData(scalar)[neighbour] - detail::fieldData(scalar)[owner]);
                const double correction =
                    gamma *
                    (integratedNormalGradient(scalar, *cell_gradient, face, diffusion_method) -
                     orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
                equation.source[static_cast<std::size_t>(neighbour)] -= correction;
            }
            continue;
        }

        const auto& condition = detail::FieldAccess::condition(scalar, face);
        if (condition.type == BoundaryType::FixedValue) {
            if (assemble_matrix)
                equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.source[static_cast<std::size_t>(owner)] += coefficient * condition.value;
            if (correctedDiffusion(diffusion_method)) {
                const double orthogonal = detail::meshData(mesh).face_orthogonal_coefficients[f] *
                                          (condition.value - detail::fieldData(scalar)[owner]);
                equation.source[static_cast<std::size_t>(owner)] +=
                    gamma *
                    (integratedNormalGradient(scalar, *cell_gradient, face, diffusion_method) -
                     orthogonal);
            }
        } else if (condition.type == BoundaryType::FixedGradient) {
            equation.source[static_cast<std::size_t>(owner)] +=
                gamma * condition.value * detail::meshData(mesh).face_areas[f];
        }
    }
}

void addVectorDiffusion(VectorDiscreteEquation& equation, FaceDiffusivity diffusivity,
                        const VectorField& vector, GradientMethod gradient_method,
                        DiffusionMethod diffusion_method) {
    const Mesh& mesh = vector.mesh();
    requireEquation(equation, mesh);
    requireField(vector, mesh, FieldLocation::Cell, "vector");
    std::optional<TensorField> cell_gradient;
    if (correctedDiffusion(diffusion_method)) {
        cell_gradient.emplace(mesh, FieldLocation::Cell, "grad(" + vector.name() + ')');
        gradient(vector, *cell_gradient, gradient_method);
    } else if (diffusion_method != DiffusionMethod::Orthogonal) {
        throw std::invalid_argument("unsupported vector diffusion method");
    }

    for (Index face : detail::meshData(mesh).owned_faces) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = detail::meshData(mesh).face_owner[f];
        const Index neighbour = detail::meshData(mesh).face_neighbour[f];
        const double gamma = diffusivity[face];
        if (gamma < 0.0 || !std::isfinite(gamma)) {
            throw std::invalid_argument("face diffusivity must be non-negative and finite");
        }
        if (gamma == 0.0) {
            continue;
        }
        const double coefficient = gamma * detail::meshData(mesh).face_orthogonal_coefficients[f];
        if (neighbour != invalid_index) {
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.diagonal[static_cast<std::size_t>(neighbour)] += coefficient;
            equation.upper[f] -= coefficient;
            equation.lower[f] -= coefficient;
            if (correctedDiffusion(diffusion_method)) {
                const Vec3 orthogonal =
                    detail::meshData(mesh).face_orthogonal_coefficients[f] *
                    (detail::fieldData(vector)[neighbour] - detail::fieldData(vector)[owner]);
                const Vec3 correction = gamma * (integratedNormalGradient(vector, *cell_gradient,
                                                                          face, diffusion_method) -
                                                 orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
                equation.source[static_cast<std::size_t>(neighbour)] -= correction;
            }
            continue;
        }

        const auto& condition = detail::FieldAccess::condition(vector, face);
        if (condition.type == BoundaryType::FixedValue ||
            condition.type == BoundaryType::Symmetry) {
            const Vec3 boundary_value = boundaryFaceValue(vector, face, -1.0);
            equation.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            equation.source[static_cast<std::size_t>(owner)] += coefficient * boundary_value;
            if (correctedDiffusion(diffusion_method)) {
                const Vec3 orthogonal = detail::meshData(mesh).face_orthogonal_coefficients[f] *
                                        (boundary_value - detail::fieldData(vector)[owner]);
                const Vec3 correction = gamma * (integratedNormalGradient(vector, *cell_gradient,
                                                                          face, diffusion_method) -
                                                 orthogonal);
                equation.source[static_cast<std::size_t>(owner)] += correction;
            }
        } else if (condition.type == BoundaryType::FixedGradient) {
            equation.source[static_cast<std::size_t>(owner)] +=
                gamma * detail::meshData(mesh).face_areas[f] * condition.value;
        }
    }
}

} // 匿名命名空间

void refreshDiffusionRhs(ScalarDiscreteEquation& equation, const ScalarField& diffusivity,
                         const ScalarField& scalar, const VectorField* gradient,
                         DiffusionMethod method) {
    requireField(diffusivity, scalar.mesh(), FieldLocation::Face, "face diffusivity");
    std::fill(equation.source.begin(), equation.source.end(), 0.0);
    addScalarDiffusion(equation, {0.0, &diffusivity}, scalar, GradientMethod::GreenGauss, method,
                       gradient, false);
}

void addDiffusion(ScalarDiscreteEquation& equation, double diffusivity, const ScalarField& scalar,
                  GradientMethod gradient_method, DiffusionMethod diffusion_method) {
    if (diffusivity < 0.0 || !std::isfinite(diffusivity)) {
        throw std::invalid_argument("diffusivity must be non-negative and finite");
    }
    if (diffusivity != 0.0) {
        addScalarDiffusion(equation, {diffusivity, nullptr}, scalar, gradient_method,
                           diffusion_method);
    }
}

void addDiffusion(ScalarDiscreteEquation& equation, const ScalarField& face_diffusivity,
                  const ScalarField& scalar, GradientMethod gradient_method,
                  DiffusionMethod diffusion_method) {
    requireField(face_diffusivity, scalar.mesh(), FieldLocation::Face, "face diffusivity");
    addScalarDiffusion(equation, {0.0, &face_diffusivity}, scalar, gradient_method,
                       diffusion_method);
}

void addDiffusion(VectorDiscreteEquation& equation, double diffusivity, const VectorField& vector,
                  GradientMethod gradient_method, DiffusionMethod diffusion_method) {
    if (diffusivity < 0.0 || !std::isfinite(diffusivity)) {
        throw std::invalid_argument("diffusivity must be non-negative and finite");
    }
    if (diffusivity != 0.0) {
        addVectorDiffusion(equation, {diffusivity, nullptr}, vector, gradient_method,
                           diffusion_method);
    }
}

void addDiffusion(ScalarDiscreteEquation& equation, double diffusivity, const ScalarField& scalar,
                  const VectorField& reconstructed_gradient, DiffusionMethod method) {
    if (diffusivity < 0.0 || !std::isfinite(diffusivity))
        throw std::invalid_argument("diffusivity must be non-negative and finite");
    addScalarDiffusion(equation, {diffusivity, nullptr}, scalar, GradientMethod::GreenGauss, method,
                       &reconstructed_gradient);
}

void addDiffusion(ScalarDiscreteEquation& equation, const ScalarField& face_diffusivity,
                  const ScalarField& scalar, const VectorField& reconstructed_gradient,
                  DiffusionMethod method) {
    requireField(face_diffusivity, scalar.mesh(), FieldLocation::Face, "face diffusivity");
    addScalarDiffusion(equation, {0.0, &face_diffusivity}, scalar, GradientMethod::GreenGauss,
                       method, &reconstructed_gradient);
}

void addDiffusion(VectorDiscreteEquation& equation, const ScalarField& face_diffusivity,
                  const VectorField& vector, GradientMethod gradient_method,
                  DiffusionMethod diffusion_method) {
    requireField(face_diffusivity, vector.mesh(), FieldLocation::Face, "face diffusivity");
    addVectorDiffusion(equation, {0.0, &face_diffusivity}, vector, gradient_method,
                       diffusion_method);
}

void addTimeDerivative(ScalarDiscreteEquation& equation, const ScalarField& previous, double dt,
                       double density, TimeMethod method, const ScalarField* older) {
    addTimeDerivativeImpl(equation, previous, dt, density, method, older);
}

void addTimeDerivative(ScalarDiscreteEquation& equation, const ScalarField& previous, double dt,
                       const ScalarField& volumetric_capacity, TimeMethod method,
                       const ScalarField* older) {
    addFieldTimeDerivativeImpl(equation, previous, dt, volumetric_capacity, method, older);
}

void addTimeDerivative(VectorDiscreteEquation& equation, const VectorField& previous, double dt,
                       const ScalarField& density, TimeMethod method, const VectorField* older) {
    addFieldTimeDerivativeImpl(equation, previous, dt, density, method, older);
}

void addTimeDerivative(VectorDiscreteEquation& equation, const VectorField& previous, double dt,
                       double density, TimeMethod method, const VectorField* older) {
    addTimeDerivativeImpl(equation, previous, dt, density, method, older);
}

} // namespace babelsim
