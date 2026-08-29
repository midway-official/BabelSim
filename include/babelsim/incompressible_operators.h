#pragma once

#include "babelsim/equation.h"
#include "babelsim/field.h"
#include "babelsim/parallel.h"

namespace babelsim {

// 保留为具名算法组件的 CFD 专用算子。它从 collocated cell 速度与压力场生成对压力
// 稳定的面通量。
struct MomentumInterpolation {
    static void apply(
        const Mesh& mesh,
        const VectorField& velocity,
        const ScalarField& pressure,
        const ScalarField& mobility,
        const VectorField& pressure_gradient,
        ScalarField& face_flux);
};

// 压力修正同样是领域算子，而不是通用 Laplacian。组装与施加共享相同 mobility 和面
// 几何，使 SIMPLE 不重复实现修正代数。
struct PressureCorrection {
    static void assemble(
        ScalarEquation& equation,
        ScalarField& correction,
        const Mesh& mesh,
        const ScalarField& face_flux,
        const ScalarField& mobility,
        const ScalarField& pressure,
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
};

}  // babelsim 命名空间
