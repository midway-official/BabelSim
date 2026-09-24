#include "babelsim/solver_control.h"

#include <cmath>
#include <stdexcept>

namespace babelsim {

void LinearSolverConfig::validate() const {
    if (!(absolute_tolerance > 0.0) || !std::isfinite(absolute_tolerance) ||
        !(relative_tolerance > 0.0) || !std::isfinite(relative_tolerance) ||
        max_iterations <= 0) {
        throw std::invalid_argument("linear solver tolerances or iteration limit are invalid");
    }
    switch (solver) {
    case LinearSolverType::ConjugateGradient:
    case LinearSolverType::BiCGSTAB:
    case LinearSolverType::GMRES:
    case LinearSolverType::FGMRES:
        break;
    default:
        throw std::invalid_argument("unsupported Krylov solver type");
    }
    switch (preconditioner) {
    case PreconditionerType::None:
    case PreconditionerType::IncompleteCholesky:
    case PreconditionerType::Hypre:
    case PreconditionerType::GAMG:
    case PreconditionerType::BlockJacobi:
    case PreconditionerType::ASM:
    case PreconditionerType::Jacobi:
        break;
    default:
        throw std::invalid_argument("unsupported PETSc preconditioner type");
    }
    if (preconditioner == PreconditionerType::GAMG &&
        (amg_max_levels <= 0 || amg_coarse_size <= 0 || amg_smoothing_steps <= 0)) {
        throw std::invalid_argument("GAMG hierarchy controls must be positive");
    }
}

}  // namespace babelsim
