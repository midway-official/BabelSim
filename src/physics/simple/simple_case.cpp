#include "internal/case_runner.h"

#include "babelsim/field_io.h"
#include "babelsim/mesh_io.h"
#include "babelsim/parallel_writer.h"
#include "babelsim/simple.h"
#include "babelsim/simple_io.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace babelsim {

int runSimpleCase(
    const CaseDefinition& definition,
    const ParallelContext& parallel,
    const std::string& requested_time)
{
    Mesh mesh = readDistributedMesh(definition.mesh_file, parallel);
    IncompressibleFields fields(mesh);
    readFieldFile(definition.fields_directory / "U.field", fields.velocity);
    readFieldFile(definition.fields_directory / "p.field", fields.pressure);
    fields.face_flux.fill(0.0);

    const SimpleCaseControl controls = readSimpleCase(definition);
    RunTime run_time = RunTime::forMesh(
        mesh, simpleRunTimeControl(controls.runtime, controls.simple));
    SimpleSolver solver(run_time, fields, controls.fluid, controls.simple);
    SimpleIterationResult result;
    int completed = 0;
    for (int iteration = 1; iteration <= controls.simple.max_iterations; ++iteration) {
        result = solver.iterate();
        completed = iteration;
        if (run_time.primary() &&
            (iteration == 1 || iteration % 100 == 0 || result.converged || !result.healthy)) {
            std::cout << "SIMPLE " << std::setw(5) << iteration << std::scientific
                      << std::setprecision(4) << " mass=" << result.continuity.relative
                      << " dU=" << result.relative_velocity_change
                      << " linP=" << result.pressure.relative_residual
                      << " linear=" << (result.linear_converged ? "ok" : "inexact") << '\n';
        }
        if (!result.healthy || result.converged) break;
    }
    if (!result.healthy) throw std::runtime_error("SIMPLE encountered a numerical failure");

    const auto time_directory = outputTimeDirectory(definition, requested_time);
    writeOwnedFieldCsv(time_directory, fields.velocity, parallel);
    writeOwnedFieldCsv(time_directory, fields.pressure, parallel);
    writeOwnedResultMetadata(
        time_directory, mesh, parallel, time_directory.filename().string(),
        {{fields.velocity.name(), "vector", fields.velocity.location()},
         {fields.pressure.name(), "scalar", fields.pressure.location()}});
    if (run_time.primary()) {
        std::cout << "SIMPLE completed=" << completed
                  << " converged=" << (result.converged ? "true" : "false")
                  << " linear=" << (result.linear_converged ? "ok" : "inexact")
                  << " mass=" << std::scientific << result.continuity.relative
                  << " dU=" << result.relative_velocity_change << '\n';
    }
    return result.converged ? 0 : 2;
}

}  // babelsim 命名空间
