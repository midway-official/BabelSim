#include "babelsim/thermal.h"

#include "babelsim/fvm.h"

#include <cmath>
#include <stdexcept>

namespace babelsim {

namespace {

void requireHeatFields(
    const RunTime& run_time,
    const ScalarField& temperature,
    const ThermalFieldProperties& material)
{
    const Mesh& mesh = run_time.mesh();
    if (&temperature.mesh() != &mesh || temperature.location() != FieldLocation::Cell ||
        &material.volumetric_heat_capacity.mesh() != &mesh ||
        material.volumetric_heat_capacity.location() != FieldLocation::Cell ||
        &material.conductivity.mesh() != &mesh ||
        material.conductivity.location() != FieldLocation::Cell ||
        &material.volumetric_source.mesh() != &mesh ||
        material.volumetric_source.location() != FieldLocation::Cell) {
        throw std::invalid_argument("thermal Fields must be cell fields on the run mesh");
    }
}

}  // 匿名命名空间

void ThermalProperties::validate() const {
    if (!(density > 0.0) || !(heat_capacity > 0.0) || !(conductivity >= 0.0) ||
        !std::isfinite(density) || !std::isfinite(heat_capacity) ||
        !std::isfinite(conductivity)) {
        throw std::invalid_argument("thermal properties must be finite and physically valid");
    }
}

HeatResult solveTransientHeat(
    RunTime& run_time,
    ScalarField& temperature,
    const ThermalProperties& material,
    double volumetric_source)
{
    material.validate();
    if (!std::isfinite(volumetric_source)) {
        throw std::invalid_argument("heat source must be finite");
    }
    if (&temperature.mesh() != &run_time.mesh() ||
        temperature.location() != FieldLocation::Cell) {
        throw std::invalid_argument("temperature must be a cell field on the run mesh");
    }

    HeatResult result;
    while (run_time.loop()) {
        result.linear = solve(
            fvm::ddt(material.volumetricHeatCapacity(), temperature) ==
                fvm::laplacian(material.conductivity, temperature) +
                fvm::source(volumetric_source));
        ++result.steps;
        if (!result.linear.converged()) return result;
    }
    result.converged = result.steps > 0 && result.linear.converged();
    return result;
}

SolveResult solveHeatStep(
    RunTime& run_time,
    ScalarField& temperature,
    const ThermalFieldProperties& material)
{
    requireHeatFields(run_time, temperature, material);
    return solve(
        fvm::ddt(material.volumetric_heat_capacity, temperature) ==
            fvm::laplacian(material.conductivity, temperature) +
            fvm::source(material.volumetric_source));
}

HeatResult solveTransientHeat(
    RunTime& run_time,
    ScalarField& temperature,
    const ThermalFieldProperties& material)
{
    requireHeatFields(run_time, temperature, material);
    HeatResult result;
    while (run_time.loop()) {
        result.linear = solveHeatStep(run_time, temperature, material);
        ++result.steps;
        if (!result.linear.converged()) return result;
    }
    result.converged = result.steps > 0 && result.linear.converged();
    return result;
}

}  // babelsim 命名空间
