#pragma once

#include "babelsim/field.h"

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
    // Version 2: ordered hexahedron vertices at output time, indexed by global cell id.
    // Empty for legacy version 1 files, which cannot prove geometric provenance.
    std::vector<std::array<Vec3, 8>> cell_vertices;
};
}  // babelsim 命名空间
