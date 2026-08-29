#pragma once

#include "babelsim/distributed_solver.h"
#include "babelsim/incompressible_operators.h"
#include "babelsim/linear_solver.h"
#include "babelsim/methods.h"
#include "babelsim/operators.h"

#include <array>
#include <vector>

namespace babelsim {

struct FluidProperties {
    double density = 1.0;
    double dynamic_viscosity = 1e-3;

    void validate() const;
};

struct SimpleControl {
    int max_iterations = 1000;
    int non_orthogonal_corrections = 1;
    double velocity_relaxation = 0.7;
    double pressure_relaxation = 0.3;
    double continuity_tolerance = 1e-8;
    double velocity_tolerance = 1e-7;
    LinearSolverConfig velocity_solver{};
    LinearSolverConfig pressure_solver{
        LinearSolverType::ConjugateGradient,
        PreconditionerType::IncompleteCholesky,
        1e-12,
        1e-8,
        1000,
        false,
    };

    void validate() const;
};

struct TimeState {
    double dt = 0.0;
    const VectorField* previous = nullptr;
    const VectorField* older = nullptr;
};

struct IncompressibleFields {
    explicit IncompressibleFields(const Mesh& mesh)
        : velocity(mesh, FieldLocation::Cell, "U"),
          pressure(mesh, FieldLocation::Cell, "p"),
          face_flux(mesh, FieldLocation::Face, "phi")
    {}

    VectorField velocity;
    ScalarField pressure;
    ScalarField face_flux;
};

struct ContinuityMetrics {
    double l1 = 0.0;
    double l2 = 0.0;
    double maximum = 0.0;
    double relative = 0.0;
};

struct SimpleIterationResult {
    std::array<SolveResult, 3> velocity;
    SolveResult pressure;
    ContinuityMetrics continuity;
    double relative_velocity_change = 0.0;
    double relative_pressure_correction = 0.0;
    bool healthy = false;
    // 本次迭代的所有线性系统均达到内层容差时为 true。MaxIterations
    // 仍是可诊断的内层状态，不应改变所有 rank 的外迭代停止时刻。
    bool linear_converged = false;
    // 外迭代的物理收敛状态，由全局归约后的连续性和速度变化量决定。
    bool converged = false;
};

class SimpleSolver {
public:
    SimpleSolver(
        IncompressibleFields& fields,
        FluidProperties fluid,
        Methods methods,
        SimpleControl control);
    SimpleSolver(
        IncompressibleFields& fields,
        FluidProperties fluid,
        Methods methods,
        SimpleControl control,
        ParallelContext parallel);

    SimpleIterationResult iterate(const TimeState& time = {});

private:
    void assembleMomentum(const TimeState& time);
    std::array<SolveResult, 3> solveMomentum();
    void momentumInterpolation();
    void assemblePressureCorrection();
    void correctPressureAndVelocity();
    ContinuityMetrics continuity() const;

    IncompressibleFields& fields_;
    const Mesh& mesh_;
    FluidProperties fluid_;
    Methods methods_;
    SimpleControl control_;
    ParallelContext parallel_;
    std::unique_ptr<HaloExchange> halo_;
    SparseAssembly momentum_assembly_;
    SparseAssembly pressure_assembly_;
    PreparedLinearSolver velocity_linear_solver_;
    PreparedLinearSolver pressure_linear_solver_;
    std::unique_ptr<DistributedLinearSolver> distributed_velocity_solver_;
    std::unique_ptr<DistributedLinearSolver> distributed_pressure_solver_;
    VectorEquation momentum_;
    ScalarEquation pressure_equation_;
    ScalarField pressure_correction_;
    VectorField pressure_gradient_;
    VectorField correction_gradient_;
    ScalarField mass_flux_;
    ScalarField mobility_;
    std::vector<Vec3> previous_velocity_;
    std::array<Eigen::VectorXd, 3> momentum_source_;
    std::array<Eigen::VectorXd, 3> velocity_solution_;
    Eigen::VectorXd pressure_source_;
    Eigen::VectorXd pressure_solution_;
    bool velocity_pattern_ready_ = false;
    bool pressure_pattern_ready_ = false;
    bool has_fixed_pressure_ = false;
};

}  // babelsim 命名空间
