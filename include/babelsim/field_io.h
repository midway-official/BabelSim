#pragma once

#include "babelsim/field.h"

#include <filesystem>

namespace babelsim {

// Field 文件只描述数值和边界条件。网格拓扑与几何保留在 mesh 文件中，Field 仍是
// 普通连续数组。
void readFieldFile(const std::filesystem::path& path, ScalarField& field);
void readFieldFile(const std::filesystem::path& path, VectorField& field);

}  // babelsim 命名空间
