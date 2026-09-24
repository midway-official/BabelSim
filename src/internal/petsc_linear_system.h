#pragma once

#include "babelsim/parallel.h"
#include "babelsim/solver_control.h"
#include "internal/petsc_assembly_plan.h"

#include <petscksp.h>

#include <array>
#include <string>
#include <vector>

namespace babelsim::detail {

struct PetscSolveMetrics {
    SolveResult result;
    double rhs_seconds = 0.0;
    double solve_seconds = 0.0;
    double residual_seconds = 0.0;
};

struct PetscMatrixUpdateMetrics { bool changed = false; double seconds = 0.0; };
struct PetscSetupMetrics { bool preconditioner_was_setup = false; double seconds = 0.0; };

// A long-lived PETSc matrix/vector/KSP set for one logical BabelSim equation.
class PetscLinearSystem {
public:
    PetscLinearSystem(std::string identity, const PetscIndexMap& index_map,
                      const PetscAssemblyPlan& plan, const ParallelContext& parallel,
                      std::size_t right_hand_sides);
    ~PetscLinearSystem();
    PetscLinearSystem(const PetscLinearSystem&) = delete;
    PetscLinearSystem& operator=(const PetscLinearSystem&) = delete;

    PetscMatrixUpdateMetrics updateMatrix(const std::vector<PetscScalar>& matrix_values);
    PetscSetupMetrics prepare(const LinearSolverConfig& config);
    PetscSolveMetrics solveRhs(const std::vector<PetscScalar>& rhs_values,
                               const std::vector<double>& initial_guess,
                               std::vector<double>& solution,
                               std::size_t right_hand_side = 0);

    const std::string& identity() const { return m_identity; }
    std::uint64_t matrixUpdates() const { return m_matrix_updates; }
    std::uint64_t rhsOnlySolves() const { return m_rhs_only_solves; }
    std::uint64_t preconditionerSetups() const { return m_preconditioner_setups; }
    std::uint64_t trueResidualChecks() const { return m_true_residual_checks; }

private:
    struct ConvergenceContext { PetscReal target = 0.0; };
    static PetscErrorCode convergenceTest(KSP, PetscInt iteration, PetscReal norm,
                                          KSPConvergedReason* reason, void* context);

    void configure(const LinearSolverConfig& config);
    void setupPreconditioner(const LinearSolverConfig& config);
    void configureBlockJacobiSubsolvers(const LinearSolverConfig& config);
    void setVectorValues(Vec vector, const std::vector<PetscScalar>& values);
    double vectorNorm(Vec vector) const;
    void computeTrueResidual(Vec rhs, Vec solution, double& norm);
    void release() noexcept;

    std::string m_identity;
    const PetscIndexMap* m_index_map;
    const PetscAssemblyPlan* m_plan;
    ParallelContext m_parallel;
    Mat m_matrix = nullptr;
    KSP m_ksp = nullptr;
    std::array<Vec, 3> m_rhs{{nullptr, nullptr, nullptr}};
    std::array<Vec, 3> m_solution{{nullptr, nullptr, nullptr}};
    Vec m_residual = nullptr;
    std::vector<PetscScalar> m_previous_matrix_values;
    std::vector<double> m_work_solution;
    std::size_t m_right_hand_sides = 0;
    bool m_matrix_ready = false;
    bool m_setup_pending = true;
    bool m_configured = false;
    LinearSolverConfig m_config;
    ConvergenceContext m_convergence;
    std::uint64_t m_matrix_updates = 0;
    std::uint64_t m_rhs_only_solves = 0;
    std::uint64_t m_preconditioner_setups = 0;
    std::uint64_t m_true_residual_checks = 0;
};

}  // namespace babelsim::detail
