#pragma once

#include "babelsim/field.h"
#include "babelsim/methods.h"
#include "babelsim/runtime.h"

namespace babelsim {

// Rhie-Chow 面重构的持久工作场。它有明确的数值含义，并由 SIMPLE 创建一次、
// 在每个外迭代复用，避免原先每次动量插值分配三个完整 Field。
struct MomentumInterpolationWorkspace {
    explicit MomentumInterpolationWorkspace(const Mesh& mesh)
        : pressure_response(mesh, FieldLocation::Cell, "rAUGradP"),
          face_pressure_response(mesh, FieldLocation::Face, "rAUGradPFace"),
          face_mobility(mesh, FieldLocation::Face, "rAUFace")
    {}

    VectorField pressure_response;
    VectorField face_pressure_response;
    ScalarField face_mobility;
};

// 保留为具名算法组件的 CFD 专用算子。它从同位单元速度与压力场生成对压力
// 稳定的面通量。
struct MomentumInterpolation {
    // 兼容早期正交调用；新求解器应显式传入三类 Method。
    static void apply(
        RunTime& run_time,
        const Mesh& mesh,
        const VectorField& velocity,
        const ScalarField& pressure,
        const ScalarField& mobility,
        const VectorField& pressure_gradient,
        ScalarField& face_flux,
        MomentumInterpolationWorkspace& workspace);
    static void apply(
        RunTime& run_time,
        const Mesh& mesh,
        const VectorField& velocity,
        const ScalarField& pressure,
        const ScalarField& mobility,
        const VectorField& pressure_gradient,
        ScalarField& face_flux,
        MomentumInterpolationWorkspace& workspace,
        InterpolationMethod interpolation_method,
        GradientMethod gradient_method,
        DiffusionMethod diffusion_method);
};

// 压力修正同样是领域算子，而不是通用拉普拉斯算子。其方程通过 fvm::laplacian
// 和 fvc::div 进入 RunTime；专用部分只保留 SIMPLE 的参考压强、非正交循环和一致
// 的压力/速度/通量修正。
struct PressureCorrection {
    static SolveResult solve(
        RunTime& run_time,
        ScalarField& correction,
        const Mesh& mesh,
        const ScalarField& face_flux,
        const ScalarField& mobility,
        const ScalarField& pressure,
        ScalarField& face_mobility,
        ScalarField& divergence,
        bool has_fixed_pressure,
        bool reset_correction = true);

    static void apply(
        RunTime& run_time,
        const Mesh& mesh,
        double pressure_relaxation,
        ScalarField& pressure,
        VectorField& velocity,
        ScalarField& face_flux,
        const ScalarField& correction,
        const ScalarField& mobility,
        VectorField& correction_gradient,
        DiffusionMethod diffusion_method);
};

}  // babelsim 命名空间
