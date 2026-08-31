#pragma once

#include "babelsim/result.h"

#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace babelsim {

ResultData readParallelResults(
    const std::filesystem::path& time_directory,
    Index global_cell_count);

}  // babelsim 命名空间
