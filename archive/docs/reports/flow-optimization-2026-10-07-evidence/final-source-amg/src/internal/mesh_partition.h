#pragma once

#include "babelsim/mesh.h"
#include <vector>

namespace babelsim::detail {
// Deterministic graph-based ownership with exact cell-count quotas.
std::vector<Index> partitionOwners(const Mesh& mesh, int partitions);
}
