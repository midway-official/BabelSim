#include "babelsim/transport.h"

#include "babelsim/fvm.h"

#include <cmath>
#include <stdexcept>

namespace babelsim {

ScalarTransportResult solveTransientScalarTransport(
    RunTime& run_time,
    ScalarField& scalar,
    const ScalarField& face_flux,
    double storage,
    double diffusivity,
    double source)
{
    if (&scalar.mesh() != &run_time.mesh() || &face_flux.mesh() != &run_time.mesh() ||
        scalar.location() != FieldLocation::Cell || face_flux.location() != FieldLocation::Face ||
        !(storage > 0.0) || !(diffusivity >= 0.0) || !std::isfinite(source)) {
        throw std::invalid_argument("scalar transport inputs are invalid");
    }

    ScalarTransportResult result;
    while (run_time.loop()) {
        result.linear = solve(
            fvm::ddt(storage, scalar) + fvm::div(face_flux, scalar) ==
                fvm::laplacian(diffusivity, scalar) + fvm::source(source));
        ++result.steps;
        if (!result.linear.converged()) return result;
    }
    result.converged = result.steps > 0 && result.linear.converged();
    return result;
}

}  // babelsim 命名空间
