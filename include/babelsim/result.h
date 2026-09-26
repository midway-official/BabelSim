#pragma once

#include "babelsim/field.h"

#include <array>
#include <string>
#include <vector>

namespace babelsim {
// 结果文件的公开数据契约，与 MPI 通信和求解器存储无关。
struct FieldOutputInfo {
    std::string name;
    std::string type;
    FieldLocation location = FieldLocation::Cell;
};

struct ResultField {
    FieldOutputInfo info;
    int components = 0;
    // 结果格式采用实体内分量连续的顺序，不等于计算 Field 的内存布局。
    std::vector<double> values;
};

struct ResultData {
    std::string time_name;
    Index global_cell_count = 0;
    std::vector<ResultField> fields;
    // Version 3: variable-length cell vertex provenance indexed by global cell id.
    // The old fixed-size vector is retained solely for reading version 2 result
    // directories; mesh input never uses this representation.
    std::vector<std::vector<Vec3>> cell_geometry;
    std::vector<std::array<Vec3, 8>> cell_vertices;
};
}  // babelsim 命名空间
