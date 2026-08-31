#pragma once

#include "babelsim/field.h"
#include "babelsim/methods.h"

namespace babelsim::detail {

// SIMPLE/Rhie--Chow 的面通量修正核。该接口只供不可压算法的数值工作区调用，
// 以便 Solver 层不接触 face topology 或逐面存储。
void applyMomentumInterpolation(
    const ScalarField& pressure,
    const VectorField& pressure_gradient,
    const ScalarField& face_mobility,
    const VectorField& face_pressure_response,
    ScalarField& predicted_flux,
    DiffusionMethod method);

}  // babelsim::detail 命名空间
