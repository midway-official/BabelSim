#pragma once

#include "babelsim/mesh.h"

#include <filesystem>

namespace babelsim {

// Read BabelSim's face-based polyhedral mesh format (BABELSIM_MESH 3).  Each
// face record stores a variable-length vertex ring, one owner, and either one
// neighbour or a boundary patch; the file does not contain physics fields.
Mesh readMeshFile(const std::filesystem::path& path);

}  // babelsim 命名空间
