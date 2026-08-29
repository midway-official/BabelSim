#pragma once

#include "babelsim/equation.h"
#include "babelsim/field.h"
#include "babelsim/methods.h"
#include "babelsim/parallel.h"

namespace babelsim {

// 保留为具名算法组件的 CFD 专用算子。它从同位单元速度与压力场生成对压力
// 稳定的面通量。
struct MomentumInterpolation {
    // 兼容早期正交调用；新求解器应显式传入三类 Method。
    static void apply(
        const Mesh& mesh,
        const VectorField& velocity,
        const ScalarField& pressure,
        const ScalarField& mobility,
        const VectorField& pressure_gradient,
        ScalarField& face_flux);
    static void apply(
        const Mesh& mesh,
        const VectorField& velocity,
        const ScalarField& pressure,
        const ScalarField& mobility,
        const VectorField& pressure_gradient,
        ScalarField& face_flux,
        InterpolationMethod interpolation_method,
        GradientMethod gradient_method,
        DiffusionMethod diffusion_method);
};

// 压力修正同样是领域算子，而不是通用拉普拉斯算子。组装与施加共享相同 mobility 和面
// 几何，使 SIMPLE 不重复实现修正代数。
struct PressureCorrection {
    // 兼容早期正交调用；非正交求解使用下面带梯度和 DiffusionMethod 的接口。
    static void assemble(
        ScalarEquation& equation,
        ScalarField& correction,
        const Mesh& mesh,
        const ScalarField& face_flux,
        const ScalarField& mobility,
        const ScalarField& pressure,
        bool has_fixed_pressure,
        const ParallelContext& parallel);
    static void assemble(
        ScalarEquation& equation,
        ScalarField& correction,
        const Mesh& mesh,
        const ScalarField& face_flux,
        const ScalarField& mobility,
        const ScalarField& pressure,
        const VectorField* correction_gradient,
        DiffusionMethod diffusion_method,
        bool has_fixed_pressure,
        const ParallelContext& parallel);

    static void apply(
        const Mesh& mesh,
        double pressure_relaxation,
        ScalarField& pressure,
        VectorField& velocity,
        ScalarField& face_flux,
        const ScalarField& correction,
        const VectorField& correction_gradient,
        const ScalarField& mobility);
    static void apply(
        const Mesh& mesh,
        double pressure_relaxation,
        ScalarField& pressure,
        VectorField& velocity,
        ScalarField& face_flux,
        const ScalarField& correction,
        const VectorField& correction_gradient,
        const ScalarField& mobility,
        DiffusionMethod diffusion_method);
};

}  // babelsim 命名空间
