#pragma once

#include "babelsim/config.h"
#include "babelsim/solver.h"

#include <cmath>
#include <stdexcept>

namespace babelsim {

class Case;

// 稳态与瞬态 SIMPLE 共用的私有物理量和算法控制，不属于 BabelSim Solver SDK。
struct FluidProperties {
    double density = 1.0;
    double dynamic_viscosity = 1e-3;

    void validate() const {
        if (!(density > 0.0) || !std::isfinite(density) ||
            !(dynamic_viscosity > 0.0) || !std::isfinite(dynamic_viscosity)) {
            throw std::invalid_argument("density and dynamic viscosity must be positive");
        }
    }
};

struct SimpleControl {
    int max_iterations = 1000;
    int non_orthogonal_corrections = 1;
    double velocity_relaxation = 0.7;
    double pressure_relaxation = 0.3;
    double continuity_tolerance = 1e-8;
    double velocity_tolerance = 1e-7;
    // p' 是本轮用于更新 p 的未松弛压力修正。仅检查速度变化会让不同
    // 分区在压力仍变化时过早停止，因此它必须有独立的外迭代门槛。
    double pressure_correction_tolerance = 1e-6;

    void validate() const {
        if (max_iterations <= 0 || non_orthogonal_corrections < 0 ||
            non_orthogonal_corrections > 20 ||
            !(velocity_relaxation > 0.0 && velocity_relaxation <= 1.0) ||
            !(pressure_relaxation > 0.0 && pressure_relaxation <= 1.0) ||
            !(continuity_tolerance > 0.0) || !std::isfinite(continuity_tolerance) ||
            !(velocity_tolerance > 0.0) || !std::isfinite(velocity_tolerance) ||
            !(pressure_correction_tolerance > 0.0) ||
            !std::isfinite(pressure_correction_tolerance)) {
            throw std::invalid_argument("SIMPLE controls are invalid");
        }
    }
};

inline SimpleControl readSimpleControl(const Parameters& settings) {
    SimpleControl result;
    result.max_iterations = settings.integer("maxIterations", result.max_iterations);
    result.non_orthogonal_corrections = settings.integer(
        "nonOrthogonalCorrections", result.non_orthogonal_corrections);
    result.velocity_relaxation = settings.number("velocityRelaxation", result.velocity_relaxation);
    result.pressure_relaxation = settings.number("pressureRelaxation", result.pressure_relaxation);
    result.continuity_tolerance = settings.number("continuityTolerance", result.continuity_tolerance);
    result.velocity_tolerance = settings.number("velocityTolerance", result.velocity_tolerance);
    result.pressure_correction_tolerance = settings.number(
        "pressureCorrectionTolerance", result.pressure_correction_tolerance);
    result.validate();
    return result;
}

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
    SolveResult velocity;
    SolveResult pressure;
    SolveResult turbulence{SolveStatus::Converged, 0, 0.0, 0.0, 0.0};
    FluxBalance continuity;
    double relative_velocity_change = 0.0;
    double relative_pressure_correction = 0.0;
    double relative_turbulence_change = 0.0;
    bool turbulence_active = false;
    bool healthy = false;
    bool linear_converged = false;
    bool converged = false;
};

// 稳态和瞬态 SIMPLE 共享的私有湍流耦合边界。Model 的具体类型、输运变量与
// 工作场全部留在 physics/RANS；SIMPLE 只触发校正并读取统一的收敛语义。
namespace rans {
class Model;
Model* create(
    Case& problem,
    const VectorField& velocity,
    const ScalarField& face_flux,
    ScalarField& effective_viscosity,
    double density,
    double molecular_viscosity);
void destroy(Model* model) noexcept;
SolveResult correct(Model& model);
double relativeChange(const Model& model);
double tolerance(const Model& model);
const char* name(const Model& model);
}  // rans 命名空间

}  // babelsim 命名空间
