#include "babelsim/mesh_io.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const auto mesh = babelsim::readMeshFile(argv[1]);
    mesh.validate();
    double total = 0., minimum = std::numeric_limits<double>::infinity();
    double maximum = 0.;
    for (babelsim::Index c = 0; c < mesh.cellCount(); ++c) {
        const double volume = mesh.cellVolume(c);
        if (!(volume > 0.) || !std::isfinite(volume))
            throw std::runtime_error("nonpositive/nonfinite cell volume");
        minimum = std::min(minimum, volume);
        maximum = std::max(maximum, volume);
        total += volume;
    }
    for (babelsim::Index f = 0; f < mesh.faceCount(); ++f)
        if (!(mesh.faceArea(f) > 0.) || !std::isfinite(mesh.faceArea(f)))
            throw std::runtime_error("nonpositive/nonfinite face area");
    if (std::abs(total - 1.) > 5e-10)
        throw std::runtime_error("mesh volumes do not cover the unit cube");
    std::cout << std::setprecision(17) << "{\"cells\":" << mesh.cellCount()
              << ",\"faces\":" << mesh.faceCount()
              << ",\"vertices\":" << mesh.vertexCount()
              << ",\"volume\":" << total << ",\"minimumVolume\":" << minimum
              << ",\"maximumVolume\":" << maximum << ",\"validated\":true}\n";
}
