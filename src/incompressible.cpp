#include "babelsim/incompressible.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace babelsim {
namespace {

bool finitePositive(double value) {
    return value > 0.0 && std::isfinite(value);
}

bool healthy(const SolveResult& result) {
    return result.status != SolveStatus::NumericalFailure &&
        std::isfinite(result.final_residual) &&
        std::isfinite(result.relative_residual);
}

}  // namespace

void FluidProperties::validate() const {
    if (!finitePositive(density) || !finitePositive(dynamic_viscosity)) {
        throw std::invalid_argument("density and dynamic viscosity must be positive");
    }
}

void SimpleControl::validate() const {
    if (max_iterations <= 0 ||
        !(velocity_relaxation > 0.0 && velocity_relaxation <= 1.0) ||
        !(pressure_relaxation > 0.0 && pressure_relaxation <= 1.0) ||
        !finitePositive(continuity_tolerance) ||
        !finitePositive(velocity_tolerance)) {
        throw std::invalid_argument("SIMPLE controls are invalid");
    }
    velocity_solver.validate();
    pressure_solver.validate();
}

SimpleSolver::SimpleSolver(
    IncompressibleFields& fields,
    FluidProperties fluid,
    Methods methods,
    SimpleControl control)
    : SimpleSolver(fields, fluid, methods, control, ParallelContext{})
{}

SimpleSolver::SimpleSolver(
    IncompressibleFields& fields,
    FluidProperties fluid,
    Methods methods,
    SimpleControl control,
    ParallelContext parallel)
    : fields_(fields),
      mesh_(fields.velocity.mesh()),
      fluid_(fluid),
      methods_(methods),
      control_(control),
      parallel_(parallel),
      momentum_assembly_(mesh_),
      pressure_assembly_(mesh_),
      velocity_linear_solver_(control.velocity_solver),
      pressure_linear_solver_(control.pressure_solver),
      momentum_(mesh_),
      pressure_equation_(mesh_),
      pressure_correction_(mesh_, FieldLocation::Cell, "pPrime"),
      pressure_gradient_(mesh_, FieldLocation::Cell, "gradP"),
      correction_gradient_(mesh_, FieldLocation::Cell, "gradPPrime"),
      mass_flux_(mesh_, FieldLocation::Face, "rhoPhi"),
      mobility_(static_cast<std::size_t>(mesh_.cellCount()), 0.0),
      previous_velocity_(static_cast<std::size_t>(mesh_.cellCount())),
      pressure_source_(Eigen::VectorXd::Zero(mesh_.ownedCellCount())),
      pressure_solution_(Eigen::VectorXd::Zero(mesh_.ownedCellCount()))
{
    if (&fields.pressure.mesh() != &mesh_ || &fields.face_flux.mesh() != &mesh_ ||
        fields.velocity.location() != FieldLocation::Cell ||
        fields.pressure.location() != FieldLocation::Cell ||
        fields.face_flux.location() != FieldLocation::Face) {
        throw std::invalid_argument("incompressible fields are inconsistent");
    }
    fluid_.validate();
    control_.validate();
    parallel_.validate();
    if (parallel_.distributed() !=
        (mesh_.ownedCellCount() < mesh_.cellCount())) {
        throw std::invalid_argument(
            "parallel context and mesh ownership are inconsistent");
    }
    if (parallel_.distributed()) {
        halo_ = std::make_unique<HaloExchange>(mesh_, parallel_);
        distributed_velocity_solver_ = std::make_unique<DistributedLinearSolver>(
            mesh_, parallel_, control_.velocity_solver);
        distributed_pressure_solver_ = std::make_unique<DistributedLinearSolver>(
            mesh_, parallel_, control_.pressure_solver);
    }
    for (auto& solution : velocity_solution_) {
        solution.resize(mesh_.ownedCellCount());
    }
    for (auto& source : momentum_source_) {
        source.resize(mesh_.ownedCellCount());
    }
    for (Index patch = 0; patch < static_cast<Index>(mesh_.patches.size()); ++patch) {
        const auto pressure_type = fields_.pressure.boundary(patch).type;
        if (pressure_type == BoundaryType::FixedValue) {
            has_fixed_pressure_ = true;
            pressure_correction_.setBoundary(
                patch, BoundaryCondition<double>::fixedValue(0.0));
        } else if (pressure_type == BoundaryType::Symmetry) {
            pressure_correction_.setBoundary(
                patch, BoundaryCondition<double>::symmetry());
        } else {
            pressure_correction_.setBoundary(
                patch, BoundaryCondition<double>::zeroGradient());
        }
    }
    if (halo_) {
        halo_->exchange(fields_.velocity);
        halo_->exchange(fields_.pressure);
    }
    flux(fields_.velocity, fields_.face_flux, methods_.interpolation);
}

void SimpleSolver::assembleMomentum(const TimeState& time) {
    momentum_.reset();
    for (Index face = 0; face < mesh_.faceCount(); ++face) {
        mass_flux_[face] = fluid_.density * fields_.face_flux[face];
    }
    if (methods_.time != TimeMethod::Steady) {
        if (time.previous == nullptr) {
            throw std::invalid_argument("transient SIMPLE requires a previous velocity field");
        }
        addTimeDerivative(
            momentum_, *time.previous, time.dt, fluid_.density,
            methods_.time, time.older);
    }
    addConvection(
        momentum_, mass_flux_, fields_.velocity, methods_.convection);
    addDiffusion(
        momentum_, fluid_.dynamic_viscosity, fields_.velocity,
        methods_.gradient, methods_.diffusion);
    gradient(fields_.pressure, pressure_gradient_, methods_.gradient);
    for (Index cell = 0; cell < mesh_.cellCount(); ++cell) {
        const auto c = static_cast<std::size_t>(cell);
        momentum_.source[c] -= mesh_.cell_volumes[c] * pressure_gradient_[cell];
    }

    // This is TaihoCFD's algebraically scaled form of standard equation
    // under-relaxation: keep aP, scale off-diagonals and the physical source.
    const double alpha = control_.velocity_relaxation;
    for (Index face = 0; face < mesh_.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        momentum_.upper[f] *= alpha;
        momentum_.lower[f] *= alpha;
    }
    for (Index cell = 0; cell < mesh_.cellCount(); ++cell) {
        const auto c = static_cast<std::size_t>(cell);
        if (!(momentum_.diagonal[c] > 0.0) ||
            !std::isfinite(momentum_.diagonal[c])) {
            throw std::runtime_error("momentum diagonal is not positive and finite");
        }
        momentum_.source[c] =
            alpha * momentum_.source[c] +
            (1.0 - alpha) * momentum_.diagonal[c] * fields_.velocity[cell];
        mobility_[c] = mesh_.cell_volumes[c] / momentum_.diagonal[c];
    }
    if (halo_) {
        halo_->exchange(mobility_);
    }
}

std::array<SolveResult, 3> SimpleSolver::solveMomentum() {
    momentum_assembly_.update(momentum_);
    const auto& matrix = momentum_assembly_.matrix();
    assembleSource(momentum_, momentum_source_);
    if (distributed_velocity_solver_) {
        if (!velocity_pattern_ready_) {
            distributed_velocity_solver_->compute(matrix, momentum_);
            velocity_pattern_ready_ = true;
        } else {
            distributed_velocity_solver_->factorize(matrix, momentum_);
        }
    } else {
        if (!velocity_pattern_ready_) {
            velocity_linear_solver_.compute(matrix);
            velocity_pattern_ready_ = true;
        } else {
            velocity_linear_solver_.factorize(matrix);
        }
    }
    std::array<SolveResult, 3> results;
    for (std::size_t component = 0; component < 3; ++component) {
        Eigen::VectorXd& solution = velocity_solution_[component];
        for (Index cell : mesh_.owned_cells) {
            solution[mesh_.ownedIndex(cell)] = fields_.velocity[cell][component];
        }
        results[component] = distributed_velocity_solver_
            ? distributed_velocity_solver_->solve(
                  momentum_source_[component], solution)
            : velocity_linear_solver_.solve(
                  momentum_source_[component], solution);
        for (Index cell : mesh_.owned_cells) {
            fields_.velocity[cell][component] =
                solution[mesh_.ownedIndex(cell)];
        }
    }
    if (halo_) {
        halo_->exchange(fields_.velocity);
    }
    return results;
}

void SimpleSolver::momentumInterpolation() {
    gradient(fields_.pressure, pressure_gradient_, methods_.gradient);
    for (Index face = 0; face < mesh_.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh_.face_owner[f];
        const Index neighbour = mesh_.face_neighbour[f];
        if (neighbour == invalid_index) {
            const Vec3 boundary_velocity = boundaryFaceValue(
                fields_.velocity, face, fields_.face_flux[face]);
            fields_.face_flux[face] = dot(
                boundary_velocity, mesh_.face_area_vectors[f]);
            continue;
        }
        const double weight = mesh_.face_owner_weights[f];
        const double face_mobility =
            weight * mobility_[static_cast<std::size_t>(owner)] +
            (1.0 - weight) * mobility_[static_cast<std::size_t>(neighbour)];
        const Vec3 face_velocity =
            weight * fields_.velocity[owner] +
            (1.0 - weight) * fields_.velocity[neighbour];
        const Vec3 interpolated_pressure_response =
            weight * mobility_[static_cast<std::size_t>(owner)] *
                pressure_gradient_[owner] +
            (1.0 - weight) * mobility_[static_cast<std::size_t>(neighbour)] *
                pressure_gradient_[neighbour];
        fields_.face_flux[face] =
            dot(face_velocity, mesh_.face_area_vectors[f]) +
            dot(interpolated_pressure_response, mesh_.face_area_vectors[f]) -
            face_mobility * mesh_.face_orthogonal_coefficients[f] *
                (fields_.pressure[neighbour] - fields_.pressure[owner]);
    }
}

void SimpleSolver::assemblePressureCorrection() {
    pressure_equation_.reset();
    pressure_correction_.fill(0.0);
    for (Index face = 0; face < mesh_.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh_.face_owner[f];
        const Index neighbour = mesh_.face_neighbour[f];
        pressure_equation_.source[static_cast<std::size_t>(owner)] -=
            fields_.face_flux[face];
        if (neighbour != invalid_index) {
            pressure_equation_.source[static_cast<std::size_t>(neighbour)] +=
                fields_.face_flux[face];
            const double weight = mesh_.face_owner_weights[f];
            const double face_mobility =
                weight * mobility_[static_cast<std::size_t>(owner)] +
                (1.0 - weight) * mobility_[static_cast<std::size_t>(neighbour)];
            const double coefficient =
                face_mobility * mesh_.face_orthogonal_coefficients[f];
            pressure_equation_.diagonal[static_cast<std::size_t>(owner)] += coefficient;
            pressure_equation_.diagonal[static_cast<std::size_t>(neighbour)] += coefficient;
            pressure_equation_.upper[f] -= coefficient;
            pressure_equation_.lower[f] -= coefficient;
        } else if (fields_.pressure.boundary(mesh_.face_patch[f]).type ==
                   BoundaryType::FixedValue) {
            pressure_equation_.diagonal[static_cast<std::size_t>(owner)] +=
                mobility_[static_cast<std::size_t>(owner)] *
                mesh_.face_orthogonal_coefficients[f];
        }
    }

    if (!has_fixed_pressure_ && parallel_.rank == 0) {
        if (mesh_.owned_cells.empty()) {
            throw std::runtime_error("closed domain has no pressure reference cell");
        }
        const auto reference = static_cast<std::size_t>(mesh_.owned_cells.front());
        if (!(pressure_equation_.diagonal[reference] > 0.0)) {
            throw std::runtime_error("pressure reference diagonal is invalid");
        }
        pressure_equation_.diagonal[reference] +=
            pressure_equation_.diagonal[reference];
    }
}

void SimpleSolver::correctPressureAndVelocity() {
    for (Index cell : mesh_.owned_cells) {
        fields_.pressure[cell] +=
            control_.pressure_relaxation * pressure_correction_[cell];
    }
    gradient(pressure_correction_, correction_gradient_, methods_.gradient);
    for (Index cell : mesh_.owned_cells) {
        fields_.velocity[cell] -=
            mobility_[static_cast<std::size_t>(cell)] * correction_gradient_[cell];
    }

    for (Index face = 0; face < mesh_.faceCount(); ++face) {
        const auto f = static_cast<std::size_t>(face);
        const Index owner = mesh_.face_owner[f];
        const Index neighbour = mesh_.face_neighbour[f];
        if (neighbour != invalid_index) {
            const double weight = mesh_.face_owner_weights[f];
            const double face_mobility =
                weight * mobility_[static_cast<std::size_t>(owner)] +
                (1.0 - weight) * mobility_[static_cast<std::size_t>(neighbour)];
            fields_.face_flux[face] +=
                face_mobility * mesh_.face_orthogonal_coefficients[f] *
                (pressure_correction_[owner] - pressure_correction_[neighbour]);
        } else if (fields_.pressure.boundary(mesh_.face_patch[f]).type ==
                   BoundaryType::FixedValue) {
            fields_.face_flux[face] +=
                mobility_[static_cast<std::size_t>(owner)] *
                mesh_.face_orthogonal_coefficients[f] * pressure_correction_[owner];
        }
    }
    if (halo_) {
        halo_->exchange(fields_.pressure);
        halo_->exchange(fields_.velocity);
    }
}

ContinuityMetrics SimpleSolver::continuity() const {
    ContinuityMetrics result;
    double squared = 0.0;
    double scale = 0.0;
    for (Index cell : mesh_.owned_cells) {
        double imbalance = 0.0;
        for (Index face : mesh_.cell_faces[static_cast<std::size_t>(cell)]) {
            const auto f = static_cast<std::size_t>(face);
            const double outward = mesh_.face_owner[f] == cell
                ? fields_.face_flux[face]
                : -fields_.face_flux[face];
            imbalance += outward;
            scale += std::abs(outward);
        }
        result.l1 += std::abs(imbalance);
        squared += imbalance * imbalance;
        result.maximum = std::max(result.maximum, std::abs(imbalance));
    }
    const double local_sums[3] = {result.l1, squared, scale};
    double global_sums[3]{};
    parallel_.sum(local_sums, global_sums, 3);
    const double local_maximum = result.maximum;
    parallel_.maximum(&local_maximum, &result.maximum, 1);
    result.l1 = global_sums[0];
    result.l2 = std::sqrt(global_sums[1]);
    result.relative = global_sums[0] / std::max(global_sums[2], 1e-30);
    return result;
}

SimpleIterationResult SimpleSolver::iterate(const TimeState& time) {
    previous_velocity_ = fields_.velocity.values();
    assembleMomentum(time);
    SimpleIterationResult result;
    result.velocity = solveMomentum();
    momentumInterpolation();
    assemblePressureCorrection();

    pressure_assembly_.update(pressure_equation_);
    assembleSource(pressure_equation_, pressure_source_);
    const auto& pressure_matrix = pressure_assembly_.matrix();
    pressure_solution_.setZero();
    if (distributed_pressure_solver_) {
        if (!pressure_pattern_ready_) {
            distributed_pressure_solver_->compute(
                pressure_matrix, pressure_equation_);
            pressure_pattern_ready_ = true;
        } else {
            distributed_pressure_solver_->factorize(
                pressure_matrix, pressure_equation_);
        }
    } else {
        if (!pressure_pattern_ready_) {
            pressure_linear_solver_.compute(pressure_matrix);
            pressure_pattern_ready_ = true;
        } else {
            pressure_linear_solver_.factorize(pressure_matrix);
        }
    }
    result.pressure = distributed_pressure_solver_
        ? distributed_pressure_solver_->solve(
              pressure_source_, pressure_solution_)
        : pressure_linear_solver_.solve(
              pressure_source_, pressure_solution_);
    for (Index cell : mesh_.owned_cells) {
        pressure_correction_[cell] =
            pressure_solution_[mesh_.ownedIndex(cell)];
    }
    if (halo_) {
        halo_->exchange(pressure_correction_);
    }
    correctPressureAndVelocity();

    double difference_squared = 0.0;
    double velocity_squared = 0.0;
    for (Index cell : mesh_.owned_cells) {
        const Vec3 difference = fields_.velocity[cell] -
            previous_velocity_[static_cast<std::size_t>(cell)];
        difference_squared += squaredNorm(difference);
        velocity_squared += squaredNorm(fields_.velocity[cell]);
    }
    double correction_squared = 0.0;
    double pressure_squared = 0.0;
    for (Index cell : mesh_.owned_cells) {
        correction_squared += pressure_correction_[cell] * pressure_correction_[cell];
        pressure_squared += fields_.pressure[cell] * fields_.pressure[cell];
    }
    const double local_norms[4] = {
        difference_squared,
        velocity_squared,
        correction_squared,
        pressure_squared,
    };
    double global_norms[4]{};
    parallel_.sum(local_norms, global_norms, 4);
    result.relative_velocity_change =
        std::sqrt(global_norms[0]) /
        std::max(std::sqrt(global_norms[1]), 1e-30);
    result.relative_pressure_correction =
        std::sqrt(global_norms[2]) /
        std::max(std::sqrt(global_norms[3]), 1e-30);
    result.continuity = continuity();
    result.healthy =
        std::all_of(result.velocity.begin(), result.velocity.end(), healthy) &&
        healthy(result.pressure) &&
        std::isfinite(result.relative_velocity_change) &&
        std::isfinite(result.continuity.relative);
    result.converged = result.healthy &&
        std::all_of(
            result.velocity.begin(), result.velocity.end(),
            [](const SolveResult& solve_result) { return solve_result.converged(); }) &&
        result.pressure.converged() &&
        result.continuity.relative <= control_.continuity_tolerance &&
        result.relative_velocity_change <= control_.velocity_tolerance;
    return result;
}

}  // namespace babelsim
