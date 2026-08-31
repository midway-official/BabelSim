#pragma once

#include "babelsim/solver.h"

#include <array>

namespace babelsim {
// 不可压缩牛顿流体的常物性模型。若未来引入变黏度模型，仍由该层提供等价物性场，
// 而不是让 SIMPLE 接触底层存储或 Case 解析。
struct FluidProperties {
    double density = 1.0;
    double dynamic_viscosity = 1e-3;

    void validate() const;
};

// SIMPLE 只拥有算法控制；空间方法、时间方法和线性后端属于 RunTime 配置。
struct SimpleControl {
    int max_iterations = 1000;
    int non_orthogonal_corrections = 1;
    double velocity_relaxation = 0.7;
    double pressure_relaxation = 0.3;
    double continuity_tolerance = 1e-8;
    double velocity_tolerance = 1e-7;
    void validate() const;
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

struct SimpleIterationResult {
    std::array<SolveResult, 3> velocity;
    SolveResult pressure;
    FluxBalance continuity;
    double relative_velocity_change = 0.0;
    double relative_pressure_correction = 0.0;
    bool healthy = false;
    bool linear_converged = false;
    bool converged = false;
};

}  // babelsim 命名空间
