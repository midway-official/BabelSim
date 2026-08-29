#pragma once

#include "babelsim/mesh.h"

#include <filesystem>

namespace babelsim {

// Reads BabelSim's native structured-mesh text format. The format stores
// dimensions, either Cartesian bounds or explicit vertices, and six logical
// side patch records; it does not contain physics-specific fields.
Mesh readMeshFile(const std::filesystem::path& path);

}  // namespace babelsim
