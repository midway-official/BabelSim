#include "internal/compute_backend.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include "internal/petsc_assembly_plan.h"
#include "internal/petsc_index_map.h"
#include "internal/petsc_linear_system.h"
#include "internal/petsc_session.h"

#include "babelsim/parallel.h"

#include <array>
#include <chrono>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace babelsim::detail {
namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

const char* kspName(LinearSolverType type) {
    switch (type) {
    case LinearSolverType::ConjugateGradient: return "cg";
    case LinearSolverType::BiCGSTAB: return "bcgs";
    case LinearSolverType::GMRES: return "gmres";
    case LinearSolverType::FGMRES: return "fgmres";
    }
    return "unknown";
}

const char* pcName(PreconditionerType type) {
    switch (type) {
    case PreconditionerType::None: return "none";
    case PreconditionerType::IncompleteCholesky: return "icc";
    case PreconditionerType::Hypre: return "hypre-boomeramg";
    case PreconditionerType::GAMG: return "gamg";
    case PreconditionerType::BlockJacobi: return "bjacobi";
    case PreconditionerType::ASM: return "asm";
    case PreconditionerType::Jacobi: return "jacobi";
    }
    return "unknown";
}

const char* solveStatusName(SolveStatus status) {
    switch (status) {
    case SolveStatus::Converged: return "converged";
    case SolveStatus::MaxIterations: return "maxIterations";
    case SolveStatus::NumericalFailure: return "numericalFailure";
    }
    return "unknown";
}

class PetscBackend final : public ComputeBackend {
public:
    PetscBackend(const Mesh& mesh, ParallelContext parallel)
        : m_mesh(&mesh), m_parallel(std::move(parallel)) {
        mesh.validate();
        m_parallel.validate();
        // A disconnected component can own every local cell and have no
        // ghosts while still being part of a distributed global mesh.
        const bool mesh_is_distributed = ownedCellCount(mesh) < mesh.globalCellCount();
        if (m_parallel.distributed() != mesh_is_distributed)
            throw std::invalid_argument(
                "PETSc backend and mesh ownership are inconsistent (rank=" +
                std::to_string(m_parallel.rank) + ", ranks=" +
                std::to_string(m_parallel.size) + ", owned=" +
                std::to_string(ownedCellCount(mesh)) + ", local=" +
                std::to_string(mesh.cellCount()) + ")");
        ensurePetscSession();
        m_index_map = std::make_unique<PetscIndexMap>(mesh, m_parallel);
        m_assembly_plan = std::make_unique<PetscAssemblyPlan>(mesh, *m_index_map);
        if (m_parallel.distributed()) m_halo = std::make_unique<HaloExchange>(mesh, m_parallel);
    }

    void synchronize(ScalarField& field) override { synchronizeField(field); }
    void synchronize(VectorField& field) override { synchronizeField(field); }
    void synchronize(TensorField& field) override { synchronizeField(field); }

    void sum(const double* local, double* global, int count) const override {
        const Clock::time_point start = Clock::now();
        m_parallel.sum(local, global, count);
        ++m_performance.global_reductions;
        m_performance.global_reduction_seconds += secondsSince(start);
    }

    void maximum(const double* local, double* global, int count) const override {
        const Clock::time_point start = Clock::now();
        m_parallel.maximum(local, global, count);
        ++m_performance.global_reductions;
        m_performance.global_reduction_seconds += secondsSince(start);
    }

    bool all(bool local_condition) const override {
        const Clock::time_point start = Clock::now();
        const bool result = m_parallel.sum(local_condition ? 1 : 0) == m_parallel.size;
        ++m_performance.global_reductions;
        m_performance.global_reduction_seconds += secondsSince(start);
        return result;
    }

    PerformanceCounters performance() const override {
        PerformanceCounters result = m_performance;
        result.equation_systems.clear();
        result.equation_systems.reserve(m_equation_performance.size());
        for (const auto& entry : m_equation_performance)
            result.equation_systems.push_back(entry.second);
        std::sort(result.equation_systems.begin(), result.equation_systems.end(),
            [](const EquationPerformanceCounters& left, const EquationPerformanceCounters& right) {
                return left.identity < right.identity;
            });
        return result;
    }

    SolveResult solve(const ScalarDiscreteEquation& equation, ScalarField& unknown,
                      const LinearSolverConfig& config) override {
        config.validate();
        equation.validateStorage();
        const std::string identity = equationIdentity(equation.numerical_identity, unknown.name(), false);
        PetscLinearSystem& system = getSystem(identity, 1);
        recordEquationConfig(identity, config);
        std::vector<PetscScalar> matrix_values;
        std::vector<PetscScalar> rhs_values;
        const Clock::time_point assembly_start = Clock::now();
        m_assembly_plan->matrixValues(equation, matrix_values);
        m_assembly_plan->rhsValues(equation, rhs_values);
        ++m_performance.equation_assemblies;
        m_performance.assembly_seconds += secondsSince(assembly_start);

        const PetscMatrixUpdateMetrics update = system.updateMatrix(matrix_values);
        recordMatrixUpdate(identity, update);
        const PetscSetupMetrics setup = system.prepare(config);
        recordPreconditionerSetup(identity, setup);

        std::vector<double> initial(static_cast<std::size_t>(m_index_map->localSize()));
        const double* field = fieldData(unknown);
        for (Index cell : meshData(*m_mesh).owned_cells)
            initial[static_cast<std::size_t>(ownedIndex(*m_mesh, cell))] = field[cell];
        std::vector<double> solution;
        const PetscSolveMetrics metrics = system.solveRhs(rhs_values, initial, solution);
        recordSolve(identity, metrics);
        double* target = fieldData(unknown);
        for (Index cell : meshData(*m_mesh).owned_cells)
            target[cell] = solution[static_cast<std::size_t>(ownedIndex(*m_mesh, cell))];
        synchronize(unknown);
        return metrics.result;
    }

    std::array<SolveResult, 3> solve(const VectorDiscreteEquation& equation,
                                     VectorField& unknown,
                                     const LinearSolverConfig& config) override {
        config.validate();
        equation.validateStorage();
        const std::string identity = equationIdentity(equation.numerical_identity, unknown.name(), true);
        PetscLinearSystem& system = getSystem(identity, 3);
        recordEquationConfig(identity, config);
        std::vector<PetscScalar> matrix_values;
        const Clock::time_point assembly_start = Clock::now();
        m_assembly_plan->matrixValues(equation, matrix_values);
        ++m_performance.equation_assemblies;
        m_performance.assembly_seconds += secondsSince(assembly_start);
        const PetscMatrixUpdateMetrics update = system.updateMatrix(matrix_values);
        recordMatrixUpdate(identity, update);
        const PetscSetupMetrics setup = system.prepare(config);
        recordPreconditionerSetup(identity, setup);

        std::array<SolveResult, 3> results;
        Vec3* values = fieldData(unknown);
        for (std::size_t component = 0; component < 3; ++component) {
            std::vector<PetscScalar> rhs_values;
            m_assembly_plan->rhsValues(equation, component, rhs_values);
            std::vector<double> initial(static_cast<std::size_t>(m_index_map->localSize()));
            for (Index cell : meshData(*m_mesh).owned_cells) {
                const Vec3& value = values[cell];
                const double component_value = component == 0 ? value.x : component == 1 ? value.y : value.z;
                initial[static_cast<std::size_t>(ownedIndex(*m_mesh, cell))] = component_value;
            }
            std::vector<double> solution;
            const PetscSolveMetrics metrics = system.solveRhs(rhs_values, initial, solution, component);
            results[component] = metrics.result;
            recordSolve(identity, metrics);
            for (Index cell : meshData(*m_mesh).owned_cells) {
                Vec3& value = values[cell];
                const double solved = solution[static_cast<std::size_t>(ownedIndex(*m_mesh, cell))];
                if (component == 0) value.x = solved;
                else if (component == 1) value.y = solved;
                else value.z = solved;
            }
        }
        synchronize(unknown);
        return results;
    }

private:
    static std::string equationIdentity(const std::string& equation_name,
                                       const std::string& field_name, bool vector) {
        return std::string(vector ? "vector:" : "scalar:") +
            (equation_name.empty() ? field_name : equation_name) + ":" + field_name;
    }

    PetscLinearSystem& getSystem(const std::string& identity, std::size_t right_hand_sides) {
        auto found = m_systems.find(identity);
        if (found == m_systems.end()) {
            const Clock::time_point build_start = Clock::now();
            auto system = std::make_unique<PetscLinearSystem>(
                identity, *m_index_map, *m_assembly_plan, m_parallel, right_hand_sides);
            found = m_systems.emplace(identity, std::move(system)).first;
            const double build_seconds = secondsSince(build_start);
            ++m_performance.matrix_pattern_builds;
            auto& equation = equationMetrics(identity);
            ++equation.matrix_pattern_builds;
            equation.pattern_build_seconds += build_seconds;
            m_performance.pattern_build_seconds += build_seconds;
        }
        return *found->second;
    }

    EquationPerformanceCounters& equationMetrics(const std::string& identity) {
        EquationPerformanceCounters& metrics = m_equation_performance[identity];
        metrics.identity = identity;
        return metrics;
    }

    void recordEquationConfig(const std::string& identity, const LinearSolverConfig& config) {
        EquationPerformanceCounters& metrics = equationMetrics(identity);
        metrics.ksp_type = kspName(config.solver);
        metrics.pc_type = pcName(config.preconditioner);
    }

    void recordMatrixUpdate(const std::string& identity,
                            const PetscMatrixUpdateMetrics& metrics) {
        EquationPerformanceCounters& equation = equationMetrics(identity);
        m_performance.matrix_update_seconds += metrics.seconds;
        equation.matrix_update_seconds += metrics.seconds;
        if (metrics.changed) {
            ++m_performance.matrix_value_updates;
            ++equation.matrix_value_updates;
        } else {
            ++m_performance.rhs_only_solves;
            ++equation.rhs_only_solves;
        }
    }

    void recordPreconditionerSetup(const std::string& identity,
                                   const PetscSetupMetrics& metrics) {
        EquationPerformanceCounters& equation = equationMetrics(identity);
        m_performance.preconditioner_seconds += metrics.seconds;
        equation.preconditioner_seconds += metrics.seconds;
        if (metrics.preconditioner_was_setup) {
            ++m_performance.preconditioner_setups;
            ++equation.preconditioner_setups;
        }
    }

    void recordSolve(const std::string& identity, const PetscSolveMetrics& metrics) {
        EquationPerformanceCounters& equation = equationMetrics(identity);
        ++m_performance.linear_solves;
        m_performance.krylov_iterations += static_cast<std::uint64_t>(metrics.result.iterations);
        m_performance.linear_solve_seconds += metrics.solve_seconds;
        m_performance.rhs_seconds += metrics.rhs_seconds;
        m_performance.residual_check_seconds += metrics.residual_seconds;
        ++m_performance.true_residual_checks;
        ++equation.linear_solves;
        equation.krylov_iterations += static_cast<std::uint64_t>(metrics.result.iterations);
        ++equation.true_residual_checks;
        equation.last_status = solveStatusName(metrics.result.status);
        equation.last_initial_residual = metrics.result.initial_residual;
        equation.last_final_residual = metrics.result.final_residual;
        equation.last_relative_residual = metrics.result.relative_residual;
        equation.rhs_seconds += metrics.rhs_seconds;
        equation.solve_seconds += metrics.solve_seconds;
        equation.residual_check_seconds += metrics.residual_seconds;
    }

    template <typename T>
    void synchronizeField(Field<T>& field) {
        field.validateStorage();
        if (&field.mesh() != m_mesh ||
            (field.location() != FieldLocation::Cell && field.location() != FieldLocation::Face))
            throw std::invalid_argument("field does not belong to the PETSc backend mesh");
        // Cache validity is local: a different owner may have modified the
        // field. Agree before skipping this collective, including on ranks
        // with no neighbours. Count the decision in reduction diagnostics.
        if (m_halo ? all(haloValid(field)) : haloValid(field)) return;
        if (m_halo) {
            const Clock::time_point start = Clock::now();
            m_halo->exchange(field);
            ++m_performance.halo_exchanges;
            std::size_t components = 1;
            if constexpr (std::is_same_v<T, Vec3>) components = 3;
            else if constexpr (std::is_same_v<T, Tensor3>) components = 9;
            m_performance.halo_bytes += m_halo->plannedBytes(
                components, false, field.location() == FieldLocation::Face);
            m_performance.halo_seconds += secondsSince(start);
        }
        markHaloValid(field);
    }

    const Mesh* m_mesh;
    ParallelContext m_parallel;
    std::unique_ptr<PetscIndexMap> m_index_map;
    std::unique_ptr<PetscAssemblyPlan> m_assembly_plan;
    std::unique_ptr<HaloExchange> m_halo;
    std::unordered_map<std::string, std::unique_ptr<PetscLinearSystem>> m_systems;
    std::unordered_map<std::string, EquationPerformanceCounters> m_equation_performance;
    mutable PerformanceCounters m_performance;
};

}  // namespace

std::unique_ptr<ComputeBackend> makeComputeBackend(const Mesh& mesh, ParallelContext parallel) {
    return std::make_unique<PetscBackend>(mesh, std::move(parallel));
}

}  // namespace babelsim::detail
