#include "internal/compute_backend.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"

#include "babelsim/assembly.h"
#include "babelsim/distributed_solver.h"
#include "babelsim/linear_solver.h"
#include "babelsim/parallel.h"

#include <Eigen/Core>

#include <array>
#include <memory>
#include <stdexcept>
#include <utility>

namespace babelsim::detail {
namespace {

template <typename Equation>
void prepare(
    const Eigen::SparseMatrix<double>& matrix,
    const Equation& equation,
    bool& pattern_ready,
    PreparedLinearSolver& serial_solver,
    DistributedLinearSolver* distributed_solver)
{
    if (distributed_solver != nullptr) {
        if (pattern_ready) distributed_solver->factorize(matrix, equation);
        else distributed_solver->compute(matrix, equation);
    } else if (pattern_ready) {
        serial_solver.factorize(matrix);
    } else {
        serial_solver.compute(matrix);
    }
    pattern_ready = true;
}

class EigenMpiBackend final : public ComputeBackend {
public:
    EigenMpiBackend(
        const Mesh& mesh,
        const LinearSolverConfig& scalar_config,
        const LinearSolverConfig& vector_config,
        ParallelContext parallel)
        : m_mesh(&mesh),
          m_parallel(std::move(parallel)),
          m_scalar_assembly(mesh),
          m_vector_assembly(mesh),
          m_scalar_solver(scalar_config),
          m_vector_solver(vector_config),
          m_scalar_source(Eigen::VectorXd::Zero(ownedCellCount(mesh))),
          m_scalar_solution(Eigen::VectorXd::Zero(ownedCellCount(mesh)))
    {
        mesh.validate();
        m_parallel.validate();
        if (m_parallel.distributed() != (ownedCellCount(mesh) < mesh.cellCount())) {
            throw std::invalid_argument("compute backend and mesh ownership are inconsistent");
        }
        if (m_parallel.distributed()) {
            m_halo = std::make_unique<HaloExchange>(mesh, m_parallel);
            m_scalar_distributed = std::make_unique<DistributedLinearSolver>(
                mesh, m_parallel, scalar_config);
            m_vector_distributed = std::make_unique<DistributedLinearSolver>(
                mesh, m_parallel, vector_config);
        }
        for (std::size_t component = 0; component < 3; ++component) {
            m_vector_source[component].resize(ownedCellCount(mesh));
            m_vector_solution[component].resize(ownedCellCount(mesh));
        }
    }

    void synchronize(ScalarField& field) override { synchronizeField(field); }
    void synchronize(VectorField& field) override { synchronizeField(field); }
    void synchronize(TensorField& field) override { synchronizeField(field); }

    void sum(const double* local, double* global, int count) const override {
        m_parallel.sum(local, global, count);
    }

    void maximum(const double* local, double* global, int count) const override {
        m_parallel.maximum(local, global, count);
    }

    bool all(bool local_condition) const override {
        return m_parallel.sum(local_condition ? 1 : 0) == m_parallel.size;
    }

    SolveResult solve(
        const ScalarDiscreteEquation& equation, ScalarField& unknown) override
    {
        m_scalar_assembly.update(equation);
        assembleSource(equation, m_scalar_source);
        prepare(
            m_scalar_assembly.matrix(), equation, m_scalar_pattern_ready,
            m_scalar_solver, m_scalar_distributed.get());

        for (Index cell : meshData(*m_mesh).owned_cells) {
            m_scalar_solution[ownedIndex(*m_mesh, cell)] = fieldData(unknown)[cell];
        }
        const SolveResult result = m_scalar_distributed
            ? m_scalar_distributed->solve(m_scalar_source, m_scalar_solution)
            : m_scalar_solver.solve(m_scalar_source, m_scalar_solution);
        for (Index cell : meshData(*m_mesh).owned_cells) {
            fieldData(unknown)[cell] = m_scalar_solution[ownedIndex(*m_mesh, cell)];
        }
        synchronize(unknown);
        return result;
    }

    std::array<SolveResult, 3> solve(
        const VectorDiscreteEquation& equation, VectorField& unknown) override
    {
        m_vector_assembly.update(equation);
        assembleSource(equation, m_vector_source);
        prepare(
            m_vector_assembly.matrix(), equation, m_vector_pattern_ready,
            m_vector_solver, m_vector_distributed.get());

        for (Index cell : meshData(*m_mesh).owned_cells) {
            for (std::size_t component = 0; component < 3; ++component) {
                m_vector_solution[component][ownedIndex(*m_mesh, cell)] =
                    fieldData(unknown)[cell][component];
            }
        }
        std::array<SolveResult, 3> results;
        for (std::size_t component = 0; component < 3; ++component) {
            results[component] = m_vector_distributed
                ? m_vector_distributed->solve(
                      m_vector_source[component], m_vector_solution[component])
                : m_vector_solver.solve(
                      m_vector_source[component], m_vector_solution[component]);
        }
        for (Index cell : meshData(*m_mesh).owned_cells) {
            for (std::size_t component = 0; component < 3; ++component) {
                fieldData(unknown)[cell][component] =
                    m_vector_solution[component][ownedIndex(*m_mesh, cell)];
            }
        }
        synchronize(unknown);
        return results;
    }

private:
    template <typename T>
    void synchronizeField(Field<T>& field) {
        field.validateStorage();
        if (&field.mesh() != m_mesh ||
            (field.location() != FieldLocation::Cell &&
             field.location() != FieldLocation::Face)) {
            throw std::invalid_argument("field does not belong to the compute backend mesh");
        }
        if (m_halo) m_halo->exchange(field);
    }

    const Mesh* m_mesh;
    ParallelContext m_parallel;
    std::unique_ptr<HaloExchange> m_halo;
    SparseAssembly m_scalar_assembly;
    SparseAssembly m_vector_assembly;
    PreparedLinearSolver m_scalar_solver;
    PreparedLinearSolver m_vector_solver;
    std::unique_ptr<DistributedLinearSolver> m_scalar_distributed;
    std::unique_ptr<DistributedLinearSolver> m_vector_distributed;
    Eigen::VectorXd m_scalar_source;
    Eigen::VectorXd m_scalar_solution;
    std::array<Eigen::VectorXd, 3> m_vector_source;
    std::array<Eigen::VectorXd, 3> m_vector_solution;
    bool m_scalar_pattern_ready = false;
    bool m_vector_pattern_ready = false;
};

}  // 匿名命名空间

std::unique_ptr<ComputeBackend> makeComputeBackend(
    const Mesh& mesh,
    const LinearSolverConfig& scalar_solver,
    const LinearSolverConfig& vector_solver,
    ParallelContext parallel)
{
    return std::make_unique<EigenMpiBackend>(
        mesh, scalar_solver, vector_solver, std::move(parallel));
}

}  // babelsim::detail 命名空间
