#include "internal/petsc_linear_system.h"

#include "internal/petsc_session.h"

#include <petscpc.h>
#include <petscoptions.h>

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
        // skip its numerical update. NaN coefficients conservatively force it.
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

// Hypre's advanced controls have no public PETSc setters. Attach a private
// options database only while parsing them, then restore the original database
// before destroying ours (PetscObjectSetOptions does not retain a reference).
class ScopedPcOptions {
public:
    explicit ScopedPcOptions(PC pc) : m_pc(pc) {
        checkPetsc(PetscObjectGetOptions(reinterpret_cast<PetscObject>(pc), &m_previous),
                   "PetscObjectGetOptions(PC)");
        checkPetsc(PetscOptionsCreate(&m_options), "PetscOptionsCreate(PC)");
    }
    ~ScopedPcOptions() {
        if (m_attached)
            PetscObjectSetOptions(reinterpret_cast<PetscObject>(m_pc), m_previous);
        PetscOptionsDestroy(&m_options);
    }
    ScopedPcOptions(const ScopedPcOptions&) = delete;
    ScopedPcOptions& operator=(const ScopedPcOptions&) = delete;

    void set(const char* key, const char* value) {
        checkPetsc(PetscOptionsSetValue(m_options, key, value), "PetscOptionsSetValue(PC)");
    }
    void apply() {
        checkPetsc(PetscObjectSetOptions(reinterpret_cast<PetscObject>(m_pc), m_options),
                   "PetscObjectSetOptions(PC)");
        m_attached = true;
        checkPetsc(PCSetFromOptions(m_pc), "PCSetFromOptions(Hypre)");
    }

private:
    PC m_pc;
    PetscOptions m_previous = nullptr;
    PetscOptions m_options = nullptr;
    bool m_attached = false;
};

void configureHypreForCpu(PC pc, LinearSolverType solver, bool distributed) {
    ScopedPcOptions options(pc);
    // Parallel coarsening and bounded distance-two interpolation reduce the
    // work and communication of the hierarchy on general sparse operators.
    // Keep the symmetric smoother, one V cycle, and zero PC tolerance used
    // by the CG-compatible PETSc defaults. BabelSim still checks b-A*x.
    options.set("-pc_hypre_boomeramg_coarsen_type", "HMIS");
    options.set("-pc_hypre_boomeramg_interp_type", "ext+i");
    options.set("-pc_hypre_boomeramg_P_max", "4");
    // BiCGStab's fixed shadow residual can break down with the aggressively
    // coarsened hierarchy even on a nonsingular diffusion operator. Keep the
    // ordinary hierarchy for that Krylov method; CG and the GMRES variants
    // retain the cheaper aggressive level. This is a solver policy, independent
    // of geometry, equation identity, matrix size, and physical model.
    options.set("-pc_hypre_boomeramg_agg_nl",
                solver == LinearSolverType::BiCGSTAB ? "0" : "1");
    // MPI benefits from natural ordering for CG/GMRES. BiCGStab also needs
    // the ordinary CF ordering to avoid the observed shadow breakdown.
    // Set both values explicitly so a Krylov configuration change on the
    // same PC cannot retain the previous method's relaxation order.
    options.set("-pc_hypre_boomeramg_no_CF",
                distributed && solver != LinearSolverType::BiCGSTAB ? "true" : "false");
    options.apply();
}

}  // namespace

PetscLinearSystem::PetscLinearSystem(
    std::string identity, const PetscIndexMap& index_map,
    const PetscAssemblyPlan& plan, const ParallelContext& parallel,
    std::size_t right_hand_sides)
    : m_identity(std::move(identity)), m_index_map(&index_map), m_plan(&plan),
      m_parallel(parallel), m_right_hand_sides(right_hand_sides) {
    if (m_identity.empty()) throw std::invalid_argument("PETSc equation identity is empty");
    if (right_hand_sides == 0 || right_hand_sides > m_fresh_iterations.size())
        throw std::invalid_argument("PETSc linear system supports one or three right-hand sides");
    ensurePetscSession();
    const MPI_Comm communicator = parallel.communicator == MPI_COMM_NULL
        ? PETSC_COMM_SELF : parallel.communicator;
    checkPetsc(MatCreate(communicator, &m_matrix), "MatCreate");
    checkPetsc(MatSetSizes(m_matrix, index_map.localSize(), index_map.localSize(),
                           index_map.globalSize(), index_map.globalSize()), "MatSetSizes");
    checkPetsc(MatSetType(m_matrix, MATAIJ), "MatSetType(MATAIJ)");
    std::vector<PetscInt> coo_rows, coo_columns;
    plan.cooPattern(coo_rows, coo_columns);
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
    checkPetsc(VecCreateMPI(communicator, index_map.localSize(), index_map.globalSize(),
                           &m_rhs), "VecCreateMPI(rhs)");
    checkPetsc(VecDuplicate(m_rhs, &m_solution), "VecDuplicate(solution)");
    checkPetsc(VecDuplicate(m_rhs, &m_residual), "VecDuplicate(residual)");
}

PetscLinearSystem::~PetscLinearSystem() { release(); }

void PetscLinearSystem::release() noexcept {
    if (m_residual) VecDestroy(&m_residual);
    if (m_solution) VecDestroy(&m_solution);
    if (m_rhs) VecDestroy(&m_rhs);
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
    // Merge CG's two inner products on distributed systems. PETSc owns the
    // recurrence and work vectors; the final b-A*x acceptance and restart
    // budget below remain unchanged. Serial CG avoids the extra work vectors.
    if (config.solver == LinearSolverType::ConjugateGradient)
        checkPetsc(KSPCGUseSingleReduction(m_ksp,
            m_parallel.distributed() ? PETSC_TRUE : PETSC_FALSE),
            "KSPCGUseSingleReduction");
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
    m_force_rebuild = true;
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
        configureHypreForCpu(pc, config.solver, m_parallel.distributed());
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
    std::vector<PetscScalar>& matrix_values) {
    if (matrix_values.size() != static_cast<std::size_t>(m_plan->nonzeroContributions()))
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
        // PETSc has consumed the COO values. Exchange scratch storage instead
        // of copying every coefficient into the retained equality snapshot.
        m_previous_matrix_values.swap(matrix_values);
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
    if (!m_setup_pending && !m_force_rebuild) return {};
    const bool reuse = hasPreconditioner(config.preconditioner) &&
        m_preconditioner_setups != 0 && !m_force_rebuild &&
        m_matrix_updates - m_preconditioner_matrix_update <
            static_cast<std::uint64_t>(config.preconditioner_rebuild_interval) &&
        m_preconditioner_rhs_solves <
            static_cast<std::uint64_t>(config.preconditioner_rebuild_interval);
    // Many RHS solves amortize a fresh factorization/hierarchy. Limit lagged
    // work as well as matrix age, without aging unchanged-matrix RHS-only
    // calls into rebuilds. Counts are collective, so no extra reduction or
    // timing-dependent branch is needed.
    return setup(reuse);
}

PetscSetupMetrics PetscLinearSystem::setup(bool reuse) {
    PetscSetupMetrics metrics;
    const Clock::time_point setup_start = Clock::now();
    checkPetsc(KSPSetReusePreconditioner(m_ksp, reuse ? PETSC_TRUE : PETSC_FALSE),
               "KSPSetReusePreconditioner");
    // Also leave PETSc's RHS-only stage when a failed/growing lagged solve
    // requests a rebuild without another matrix value update.
    checkPetsc(KSPSetOperators(m_ksp, m_matrix, m_matrix), "KSPSetOperators(setup)");
    checkPetsc(KSPSetUp(m_ksp), "KSPSetUp");
    if (!reuse) {
        configureBlockJacobiSubsolvers(m_config);
        // Factor local blocks here rather than charging their lazy setup to
        // the first velocity component's KSPSolve time.
        checkPetsc(KSPSetUpOnBlocks(m_ksp), "KSPSetUpOnBlocks");
        metrics.preconditioner_was_setup = hasPreconditioner(m_config.preconditioner);
        if (metrics.preconditioner_was_setup) ++m_preconditioner_setups;
        m_preconditioner_matrix_update = m_matrix_updates;
        m_preconditioner_rhs_solves = 0;
        m_fresh_iterations.fill(0);
    }
    metrics.seconds = secondsSince(setup_start);
    m_setup_pending = false;
    m_force_rebuild = false;
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
    if (m_config.warm_start &&
        initial_guess.size() != static_cast<std::size_t>(m_index_map->localSize()))
        throw std::invalid_argument("PETSc initial guess does not match local ownership");
    PetscSolveMetrics metrics;
    const Clock::time_point rhs_start = Clock::now();
    Vec rhs = m_rhs;
    Vec x = m_solution;
    setVectorValues(rhs, rhs_values);
    const auto reset_guess = [&]() {
        if (m_config.warm_start) {
            PetscScalar* x_array = nullptr;
            checkPetsc(VecGetArray(x, &x_array), "VecGetArray(initial guess)");
            for (std::size_t i = 0; i < initial_guess.size(); ++i)
                x_array[i] = static_cast<PetscScalar>(initial_guess[i]);
            checkPetsc(VecRestoreArray(x, &x_array), "VecRestoreArray(initial guess)");
        } else {
            checkPetsc(VecSet(x, 0.0), "VecSet(zero initial guess)");
        }
    };
    reset_guess();
    metrics.rhs_seconds = secondsSince(rhs_start);

    double rhs_norm = 0.0;
    double initial_residual = 0.0;
    if (m_config.warm_start) {
        // Both norms are needed for the original residual scale. Queue them
        // together so distributed solves pay for one collective, not two.
        checkPetsc(MatMult(m_matrix, x, m_residual), "MatMult(initial residual)");
        checkPetsc(VecAYPX(m_residual, -1.0, rhs), "VecAYPX(initial residual)");
        PetscReal b_norm = 0.0, r_norm = 0.0;
        checkPetsc(VecNormBegin(rhs, NORM_2, &b_norm), "VecNormBegin(rhs)");
        checkPetsc(VecNormBegin(m_residual, NORM_2, &r_norm), "VecNormBegin(initial residual)");
        checkPetsc(VecNormEnd(rhs, NORM_2, &b_norm), "VecNormEnd(rhs)");
        checkPetsc(VecNormEnd(m_residual, NORM_2, &r_norm), "VecNormEnd(initial residual)");
        rhs_norm = static_cast<double>(b_norm);
        initial_residual = static_cast<double>(r_norm);
    } else {
        // For x0=0, r0=b exactly; no SpMV or second norm is necessary.
        initial_residual = rhs_norm = vectorNorm(rhs);
    }
    const double scale = residualScale(initial_residual, rhs_norm);
    const double target = residualTarget(m_config, scale);
    m_convergence.target = static_cast<PetscReal>(target);

    metrics.preconditioner_reused = hasPreconditioner(m_config.preconditioner) &&
        m_preconditioner_matrix_update < m_matrix_updates;
    if (residualConverged(initial_residual, target)) {
        // r0 is already the true residual (or exactly b for a zero guess).
        // Avoid KSPSolve and its duplicate norm/matvec for an accepted x0.
        metrics.result = SolveResult(SolveStatus::Converged, 0, initial_residual,
                                    initial_residual, initial_residual / scale);
        if (m_config.warm_start) solution = initial_guess;
        else solution.assign(rhs_values.size(), 0.0);
        ++m_true_residual_checks;
        ++metrics.true_residual_checks;
        return metrics;
    }
    int total_iterations = 0;
    for (int attempt = 0; attempt < 2; ++attempt) {
        int attempt_iterations = 0;
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        double final_residual = 0.0;
        for (int restart = 0; restart < 3; ++restart) {
            checkPetsc(KSPSetTolerances(m_ksp, 0.0, 0.0, PETSC_DEFAULT,
                static_cast<PetscInt>(m_config.max_iterations - attempt_iterations)),
                "KSPSetTolerances(remaining iterations)");
            const Clock::time_point solve_start = Clock::now();
            checkPetsc(KSPSolve(m_ksp, rhs, x), "KSPSolve");
            metrics.solve_seconds += secondsSince(solve_start);
            PetscInt iterations = 0;
            checkPetsc(KSPGetIterationNumber(m_ksp, &iterations), "KSPGetIterationNumber");
            checkPetsc(KSPGetConvergedReason(m_ksp, &reason), "KSPGetConvergedReason");
            attempt_iterations += static_cast<int>(iterations);
            total_iterations += static_cast<int>(iterations);
            const Clock::time_point residual_start = Clock::now();
            computeTrueResidual(rhs, x, final_residual);
            ++m_true_residual_checks;
            ++metrics.true_residual_checks;
            metrics.residual_seconds += secondsSince(residual_start);
            if (residualConverged(final_residual, target) || reason <= 0
                || !std::isfinite(final_residual) || attempt_iterations >= m_config.max_iterations)
                break;
            // A converged recursive Krylov norm may lie just below the true
            // norm. Continue from this x using a tighter internal target;
            // keep the original scale/acceptance target and iteration budget.
            // VecNorm and KSP reason are collective, including zero-owned ranks.
            m_convergence.target *= 0.5;
            checkPetsc(KSPSetInitialGuessNonzero(m_ksp, PETSC_TRUE),
                       "KSPSetInitialGuessNonzero(true-residual restart)");
        }
        checkPetsc(KSPSetInitialGuessNonzero(m_ksp,
            m_config.warm_start ? PETSC_TRUE : PETSC_FALSE), "KSPSetInitialGuessNonzero(restore)");
        checkPetsc(KSPSetTolerances(m_ksp, 0.0, 0.0, PETSC_DEFAULT,
            static_cast<PetscInt>(m_config.max_iterations)), "KSPSetTolerances(restore)");

        SolveStatus status = SolveStatus::NumericalFailure;
        if (residualConverged(final_residual, target)) status = SolveStatus::Converged;
        else if (reason == KSP_DIVERGED_ITS || attempt_iterations >= m_config.max_iterations)
            status = SolveStatus::MaxIterations;
        metrics.result = SolveResult(status, total_iterations, initial_residual,
                                     final_residual, final_residual / scale);
        if (status != SolveStatus::Converged && attempt == 0 &&
            metrics.preconditioner_reused) {
            // Vec norms and iteration counts are collective, so all ranks
            // take this branch together. Restore the original x0 and target;
            // a failed stale PC never weakens BabelSim's residual contract.
            metrics.recovery_setup = setup(false);
            ++metrics.retries;
            reset_guess();
            m_convergence.target = static_cast<PetscReal>(target);
            continue;
        }
        if (status == SolveStatus::Converged) {
            int& fresh = m_fresh_iterations[right_hand_side];
            if (m_preconditioner_matrix_update == m_matrix_updates)
                fresh = std::max(fresh, attempt_iterations);
            else if (attempt_iterations > m_config.preconditioner_reuse_iteration_growth *
                                      static_cast<double>(std::max(fresh, 1)))
                m_force_rebuild = true;
        }
        break;
    }
    if (total_iterations > 0) ++m_preconditioner_rhs_solves;

    PetscScalar* x_array = nullptr;
    checkPetsc(VecGetArray(x, &x_array), "VecGetArray(solution)");
    solution.resize(static_cast<std::size_t>(m_index_map->localSize()));
    for (std::size_t i = 0; i < solution.size(); ++i)
        solution[i] = static_cast<double>(PetscRealPart(x_array[i]));
    checkPetsc(VecRestoreArray(x, &x_array), "VecRestoreArray(solution)");
    return metrics;
}

}  // namespace babelsim::detail
