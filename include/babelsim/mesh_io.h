#pragma once

#include "babelsim/mesh.h"

#include <filesystem>

namespace babelsim {

// 读取 BabelSim 原生结构化网格文本格式。格式存储尺寸、笛卡尔边界或显式顶点，
// 以及六个逻辑边的 patch 记录；它不包含物理专属 Field。
Mesh readMeshFile(const std::filesystem::path& path);

}  // babelsim 命名空间
