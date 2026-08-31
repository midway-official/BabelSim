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
    std::array<Index, 3> global_dimensions{};
    std::vector<ResultField> fields;
};
}  // babelsim 命名空间
