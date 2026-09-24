#include "internal/petsc_linear_system.h"

#include "internal/petsc_session.h"

#include <petscpc.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>

namespace babelsim::detail {
namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

void checkPetsc(PetscErrorCode code, const char* operation) {
    if (code == 0) return;
    const char* message = nullptr;
    PetscErrorMessage(code, &message, nullptr);
    throw std::runtime_error(std::string(operation) + " failed with PETSc error " +
                             std::to_string(code) +
                             (message ? std::string(": ") + message : std::string{}));
}

bool sameValues(const std::vector<PetscScalar>& left,
                const std::vector<PetscScalar>& right) {
    if (left.size() != right.size()) return false;
    for (std::size_t i = 0; i < left.size(); ++i) {
        // Exact equality is intentional: only an unchanged effective matrix may
        // reuse a setup PC. NaN coefficients conservatively force an update.
        if (left[i] != right[i]) return false;
    }
    return true;
}

const char* kspType(LinearSolverType type) {
    switch (type) {
    case LinearSolverType::ConjugateGradient: return KSPCG;
    case LinearSolverType::BiCGSTAB: return KSPBCGS;
    case LinearSolverType::GMRES: return KSPGMRES;
    case LinearSolverType::FGMRES: return KSPFGMRES;
    }
    throw std::invalid_argument("unsupported PETSc Krylov solver type");
}

bool hasPreconditioner(PreconditionerType type) {
    return type != PreconditionerType::None;
}

}  // namespace

PetscLinearSystem::PetscLinearSystem(
    std::string identity, const PetscIndexMap& index_map,
    const PetscAssemblyPlan& plan, const ParallelContext& parallel,
    std::size_t right_hand_sides)
    : m_identity(std::move(identity)), m_index_map(&index_map), m_plan(&plan),
      m_parallel(parallel), m_right_hand_sides(right_hand_sides) {
    if (m_identity.empty()) throw std::invalid_argument("PETSc equation identity is empty");
    if (right_hand_sides == 0 || right_hand_sides > m_rhs.size())
        throw std::invalid_argument("PETSc linear system supports one or three right-hand sides");
    ensurePetscSession();
    const MPI_Comm communicator = parallel.communicator == MPI_COMM_NULL
        ? PETSC_COMM_SELF : parallel.communicator;
    checkPetsc(MatCreate(communicator, &m_matrix), "MatCreate");
    checkPetsc(MatSetSizes(m_matrix, index_map.localSize(), index_map.localSize(),
                           index_map.globalSize(), index_map.globalSize()), "MatSetSizes");
    checkPetsc(MatSetType(m_matrix, MATAIJ), "MatSetType(MATAIJ)");
    const auto& rows = plan.rows();
    const auto& columns = plan.columns();
    std::vector<PetscInt> coo_rows = rows;
    std::vector<PetscInt> coo_columns = columns;
    checkPetsc(MatSetPreallocationCOO(m_matrix, plan.nonzeroContributions(),
                                      coo_rows.data(), coo_columns.data()),
               "MatSetPreallocationCOO");
    checkPetsc(MatSetOption(m_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE),
               "MatSetOption(MAT_NEW_NONZERO_ALLOCATION_ERR)");
    checkPetsc(MatSetOption(m_matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE),
               "MatSetOption(MAT_KEEP_NONZERO_PATTERN)");
    checkPetsc(KSPCreate(communicator, &m_ksp), "KSPCreate");
    checkPetsc(KSPSetOperators(m_ksp, m_matrix, m_matrix), "KSPSetOperators");
    checkPetsc(KSPSetConvergenceTest(m_ksp, &PetscLinearSystem::convergenceTest,
                                     &m_convergence, nullptr), "KSPSetConvergenceTest");
    for (std::size_t i = 0; i < m_right_hand_sides; ++i) {
        checkPetsc(VecCreateMPI(communicator, index_map.localSize(), index_map.globalSize(),
                               &m_rhs[i]), "VecCreateMPI(rhs)");
        checkPetsc(VecDuplicate(m_rhs[i], &m_solution[i]), "VecDuplicate(solution)");
    }
    checkPetsc(VecDuplicate(m_rhs[0], &m_residual), "VecDuplicate(residual)");
    m_work_solution.resize(static_cast<std::size_t>(index_map.localSize()));
}

PetscLinearSystem::~PetscLinearSystem() { release(); }

void PetscLinearSystem::release() noexcept {
    if (m_residual) VecDestroy(&m_residual);
    for (std::size_t i = 0; i < m_right_hand_sides; ++i) {
        if (m_solution[i]) VecDestroy(&m_solution[i]);
        if (m_rhs[i]) VecDestroy(&m_rhs[i]);
    }
    if (m_ksp) KSPDestroy(&m_ksp);
    if (m_matrix) MatDestroy(&m_matrix);
}

PetscErrorCode PetscLinearSystem::convergenceTest(
    KSP, PetscInt, PetscReal norm, KSPConvergedReason* reason, void* context) {
    const auto* convergence = static_cast<const ConvergenceContext*>(context);
    if (!std::isfinite(static_cast<double>(norm))) {
        *reason = KSP_DIVERGED_NANORINF;
    } else if (norm <= convergence->target) {
        *reason = KSP_CONVERGED_RTOL;
    } else {
        *reason = KSP_CONVERGED_ITERATING;
    }
    return 0;
}

void PetscLinearSystem::configure(const LinearSolverConfig& config) {
    config.validate();
    if (m_configured && config == m_config) return;
    checkPetsc(KSPSetType(m_ksp, kspType(config.solver)), "KSPSetType");
    checkPetsc(KSPSetNormType(m_ksp, KSP_NORM_UNPRECONDITIONED), "KSPSetNormType");
    // BabelSim's convergence callback enforces its original true-residual
    // contract; PETSc's built-in tolerance test is disabled.
    checkPetsc(KSPSetTolerances(m_ksp, 0.0, 0.0, PETSC_DEFAULT,
                                static_cast<PetscInt>(config.max_iterations)),
               "KSPSetTolerances");
    checkPetsc(KSPSetInitialGuessNonzero(
                   m_ksp, config.warm_start ? PETSC_TRUE : PETSC_FALSE),
               "KSPSetInitialGuessNonzero");
    setupPreconditioner(config);
    m_config = config;
    m_configured = true;
    m_setup_pending = true;
}

void PetscLinearSystem::setupPreconditioner(const LinearSolverConfig& config) {
    PC pc = nullptr;
    checkPetsc(KSPGetPC(m_ksp, &pc), "KSPGetPC");
    switch (config.preconditioner) {
    case PreconditionerType::None:
        checkPetsc(PCSetType(pc, PCNONE), "PCSetType(PCNONE)");
        break;
    case PreconditionerType::IncompleteCholesky:
        if (m_parallel.distributed())
            throw std::invalid_argument(
                "PETSc ICC is sequential; choose blockJacobi, hypre, or gamg for MPI");
        checkPetsc(PCSetType(pc, PCICC), "PCSetType(PCICC)");
        break;
    case PreconditionerType::Hypre:
        checkPetsc(PCSetType(pc, PCHYPRE), "PCSetType(PCHYPRE)");
        checkPetsc(PCHYPRESetType(pc, "boomeramg"), "PCHYPRESetType(boomeramg)");
        break;
    case PreconditionerType::GAMG:
        checkPetsc(PCSetType(pc, PCGAMG), "PCSetType(PCGAMG)");
        checkPetsc(PCGAMGSetNlevels(pc, config.amg_max_levels), "PCGAMGSetNlevels");
        checkPetsc(PCGAMGSetNSmooths(pc, config.amg_smoothing_steps), "PCGAMGSetNSmooths");
        checkPetsc(PCGAMGSetCoarseEqLim(pc, config.amg_coarse_size), "PCGAMGSetCoarseEqLim");
        break;
    case PreconditionerType::BlockJacobi:
        checkPetsc(PCSetType(pc, PCBJACOBI), "PCSetType(PCBJACOBI)");
        break;
    case PreconditionerType::ASM:
        checkPetsc(PCSetType(pc, PCASM), "PCSetType(PCASM)");
        break;
    case PreconditionerType::Jacobi:
        checkPetsc(PCSetType(pc, PCJACOBI), "PCSetType(PCJACOBI)");
        break;
    }
}

void PetscLinearSystem::setVectorValues(
    Vec vector, const std::vector<PetscScalar>& values) {
    PetscInt local_size = 0;
    checkPetsc(VecGetLocalSize(vector, &local_size), "VecGetLocalSize");
    if (values.size() != static_cast<std::size_t>(local_size))
        throw std::invalid_argument("PETSc vector values do not match local ownership");
    PetscScalar* array = nullptr;
    checkPetsc(VecGetArray(vector, &array), "VecGetArray");
    std::copy(values.begin(), values.end(), array);
    checkPetsc(VecRestoreArray(vector, &array), "VecRestoreArray");
}

double PetscLinearSystem::vectorNorm(Vec vector) const {
    PetscReal norm = 0.0;
    checkPetsc(VecNorm(vector, NORM_2, &norm), "VecNorm");
    return static_cast<double>(norm);
}

void PetscLinearSystem::computeTrueResidual(Vec rhs, Vec solution, double& norm) {
    checkPetsc(MatMult(m_matrix, solution, m_residual), "MatMult(true residual)");
    checkPetsc(VecAYPX(m_residual, -1.0, rhs), "VecAYPX(true residual)"); // r = b - A*x
    norm = vectorNorm(m_residual);
}

PetscMatrixUpdateMetrics PetscLinearSystem::updateMatrix(
    const std::vector<PetscScalar>& matrix_values) {
    if (matrix_values.size() != m_plan->rows().size())
        throw std::invalid_argument("PETSc matrix values do not match fixed COO pattern");
    PetscMatrixUpdateMetrics metrics;
    const Clock::time_point update_start = Clock::now();
    const bool local_same = m_matrix_ready && sameValues(matrix_values, m_previous_matrix_values);
    const bool globally_same = m_parallel.size == 1
        ? local_same : m_parallel.sum(local_same ? 1 : 0) == m_parallel.size;
    if (!globally_same) {
        checkPetsc(MatSetValuesCOO(m_matrix, matrix_values.data(), INSERT_VALUES),
                   "MatSetValuesCOO(matrix update)");
        checkPetsc(KSPSetOperators(m_ksp, m_matrix, m_matrix), "KSPSetOperators(matrix update)");
        m_previous_matrix_values = matrix_values;
        m_matrix_ready = true;
        m_setup_pending = true;
        ++m_matrix_updates;
        metrics.changed = true;
    } else {
        ++m_rhs_only_solves;
    }
    metrics.seconds = secondsSince(update_start);
    return metrics;
}

PetscSetupMetrics PetscLinearSystem::prepare(const LinearSolverConfig& config) {
    configure(config);
    PetscSetupMetrics metrics;
    if (m_setup_pending) {
        const Clock::time_point setup_start = Clock::now();
        checkPetsc(KSPSetUp(m_ksp), "KSPSetUp");
        configureBlockJacobiSubsolvers(config);
        metrics.seconds = secondsSince(setup_start);
        metrics.preconditioner_was_setup = hasPreconditioner(config.preconditioner);
        if (metrics.preconditioner_was_setup) ++m_preconditioner_setups;
        m_setup_pending = false;
    }
    return metrics;
}

void PetscLinearSystem::configureBlockJacobiSubsolvers(
    const LinearSolverConfig& config) {
    // PETSc's default BJacobi sub-PC is ILU. CG requires an SPD
    // preconditioner, and the default ILU pivot handling is not invariant to
    // very small coefficient scales. Give each local block an ICC solve for
    // CG; retain PETSc's default ILU path for nonsymmetric Krylov methods.
    if (config.preconditioner != PreconditionerType::BlockJacobi ||
        config.solver != LinearSolverType::ConjugateGradient) return;

    PC pc = nullptr;
    checkPetsc(KSPGetPC(m_ksp, &pc), "KSPGetPC(block Jacobi)");
    PetscInt local_blocks = 0;
    PetscInt first_block = 0;
    KSP* subsolvers = nullptr;
    checkPetsc(PCBJacobiGetSubKSP(pc, &local_blocks, &first_block, &subsolvers),
               "PCBJacobiGetSubKSP");
    for (PetscInt i = 0; i < local_blocks; ++i) {
        checkPetsc(KSPSetType(subsolvers[i], KSPPREONLY),
                   "KSPSetType(block Jacobi subsolver)");
        PC sub_pc = nullptr;
        checkPetsc(KSPGetPC(subsolvers[i], &sub_pc),
                   "KSPGetPC(block Jacobi subsolver)");
        checkPetsc(PCSetType(sub_pc, PCICC),
                   "PCSetType(block Jacobi ICC)");
        checkPetsc(KSPSetUp(subsolvers[i]),
                   "KSPSetUp(block Jacobi ICC)");
    }
}

PetscSolveMetrics PetscLinearSystem::solveRhs(
    const std::vector<PetscScalar>& rhs_values,
    const std::vector<double>& initial_guess,
    std::vector<double>& solution,
    std::size_t right_hand_side) {
    if (!m_matrix_ready || !m_configured || m_setup_pending)
        throw std::logic_error("PETSc linear system must update and prepare before solving RHS");
    if (right_hand_side >= m_right_hand_sides)
        throw std::out_of_range("PETSc RHS index is outside this equation system");
    if (initial_guess.size() != static_cast<std::size_t>(m_index_map->localSize()))
        throw std::invalid_argument("PETSc initial guess does not match local ownership");
    PetscSolveMetrics metrics;
    const Clock::time_point rhs_start = Clock::now();
    Vec rhs = m_rhs[right_hand_side];
    Vec x = m_solution[right_hand_side];
    setVectorValues(rhs, rhs_values);
    if (m_config.warm_start) {
        PetscScalar* x_array = nullptr;
        checkPetsc(VecGetArray(x, &x_array), "VecGetArray(initial guess)");
        for (std::size_t i = 0; i < initial_guess.size(); ++i)
            x_array[i] = static_cast<PetscScalar>(initial_guess[i]);
        checkPetsc(VecRestoreArray(x, &x_array), "VecRestoreArray(initial guess)");
    } else {
        checkPetsc(VecSet(x, 0.0), "VecSet(zero initial guess)");
    }
    metrics.rhs_seconds = secondsSince(rhs_start);

    const double rhs_norm = vectorNorm(rhs);
    double initial_residual = 0.0;
    computeTrueResidual(rhs, x, initial_residual);
    const double scale = residualScale(initial_residual, rhs_norm);
    const double target = residualTarget(m_config, scale);
    m_convergence.target = static_cast<PetscReal>(target);

    const Clock::time_point solve_start = Clock::now();
    checkPetsc(KSPSolve(m_ksp, rhs, x), "KSPSolve");
    metrics.solve_seconds = secondsSince(solve_start);
    PetscInt iterations = 0;
    PetscInt max_iterations = 0;
    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    checkPetsc(KSPGetIterationNumber(m_ksp, &iterations), "KSPGetIterationNumber");
    checkPetsc(KSPGetConvergedReason(m_ksp, &reason), "KSPGetConvergedReason");
    max_iterations = static_cast<PetscInt>(m_config.max_iterations);
    const Clock::time_point residual_start = Clock::now();
    double final_residual = 0.0;
    computeTrueResidual(rhs, x, final_residual);
    ++m_true_residual_checks;
    metrics.residual_seconds = secondsSince(residual_start);

    SolveStatus status = SolveStatus::NumericalFailure;
    if (residualConverged(final_residual, target)) status = SolveStatus::Converged;
    else if (reason == KSP_DIVERGED_ITS || iterations >= max_iterations)
        status = SolveStatus::MaxIterations;
    metrics.result = SolveResult(status, static_cast<int>(iterations), initial_residual,
                                 final_residual, final_residual / scale);

    PetscScalar* x_array = nullptr;
    checkPetsc(VecGetArray(x, &x_array), "VecGetArray(solution)");
    solution.resize(static_cast<std::size_t>(m_index_map->localSize()));
    for (std::size_t i = 0; i < solution.size(); ++i)
        solution[i] = static_cast<double>(PetscRealPart(x_array[i]));
    checkPetsc(VecRestoreArray(x, &x_array), "VecRestoreArray(solution)");
    return metrics;
}

}  // namespace babelsim::detail
