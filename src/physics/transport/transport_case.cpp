#include "internal/case_runner.h"

#include "babelsim/field_io.h"
#include "babelsim/fvc.h"
#include "babelsim/mesh_io.h"
#include "babelsim/parallel_writer.h"
#include "babelsim/transport.h"
#include "babelsim/transport_io.h"

#include <iostream>
#include <stdexcept>

namespace babelsim {

int runTransportCase(
    const CaseDefinition& definition,
    const ParallelContext& parallel,
    const std::string& requested_time)
{
    Mesh mesh = readDistributedMesh(definition.mesh_file, parallel);
    ScalarField scalar(mesh, FieldLocation::Cell, "C");
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    ScalarField face_flux(mesh, FieldLocation::Face, "phi");
    readFieldFile(definition.fields_directory / "C.field", scalar);
    readFieldFile(definition.fields_directory / "U.field", velocity);

    const TransportCaseControl controls = readTransportCase(definition);
    RunTime run_time = RunTime::forMesh(mesh, controls.runtime);
    fvc::evaluate(fvc::flux(velocity), face_flux);
    const ScalarTransportResult result = solveTransientScalarTransport(
        run_time, scalar, face_flux, controls.storage, controls.diffusivity, controls.source);
    if (!result.converged) throw std::runtime_error("transport solver did not converge");

    const auto time_directory = outputTimeDirectory(definition, requested_time);
    writeOwnedFieldCsv(time_directory, scalar, parallel);
    writeOwnedResultMetadata(
        time_directory, mesh, parallel, time_directory.filename().string(),
        {{scalar.name(), "scalar", scalar.location()}});
    if (run_time.primary()) {
        std::cout << "Transport completed steps=" << result.steps << " residual="
                  << result.linear.relative_residual << '\n';
    }
    return 0;
}

}  // babelsim 命名空间
