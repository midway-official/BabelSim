#pragma once

#include "babelsim/parallel_writer.h"

#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace babelsim {

// 从 rank-local owned-cell 文件重构出的全局有序值。每个实体内分量连续：
// values[id * components + c]。
struct ResultField {
    FieldOutputInfo info;
    int components = 0;
    std::vector<double> values;
};

struct ResultData {
    std::string time_name;
    std::array<Index, 3> global_dimensions{};
    std::vector<ResultField> fields;
};

ResultData readParallelResults(
    const std::filesystem::path& time_directory,
    Index global_cell_count);

}  // babelsim 命名空间
