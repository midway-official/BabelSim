#include "internal/case_runner.h"

#include "babelsim/field_io.h"
#include "babelsim/mesh_io.h"
#include "babelsim/parallel_writer.h"
#include "babelsim/thermal.h"
#include "babelsim/thermal_io.h"

#include <iostream>
#include <stdexcept>

namespace babelsim {

int runHeatCase(
    const CaseDefinition& definition,
    const ParallelContext& parallel,
    const std::string& requested_time)
{
    Mesh mesh = readDistributedMesh(definition.mesh_file, parallel);
    ScalarField temperature(mesh, FieldLocation::Cell, "T");
    readFieldFile(definition.fields_directory / "T.field", temperature);

    const ThermalCaseControl controls = readThermalCase(definition);
    RunTime run_time = RunTime::forMesh(mesh, controls.runtime);
    const HeatResult result = solveTransientHeat(
        run_time, temperature, controls.material, controls.volumetric_source);
    if (!result.converged) throw std::runtime_error("heat solver did not converge");

    const auto time_directory = outputTimeDirectory(definition, requested_time);
    writeOwnedFieldCsv(time_directory, temperature, parallel);
    writeOwnedResultMetadata(
        time_directory, mesh, parallel, time_directory.filename().string(),
        {{temperature.name(), "scalar", temperature.location()}});
    if (run_time.primary()) {
        std::cout << "Heat completed steps=" << result.steps << " residual="
                  << result.linear.relative_residual << '\n';
    }
    return 0;
}

}  // babelsim 命名空间
