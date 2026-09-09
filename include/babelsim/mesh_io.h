#pragma once

#include "babelsim/mesh.h"

#include <filesystem>

namespace babelsim {

// 读取 BabelSim 原生非结构六面体文本格式（BABELSIM_MESH 2）。格式显式存储
// 顶点、六面体连接和每个 patch 的边界四边形；它不包含物理专属 Field。
Mesh readMeshFile(const std::filesystem::path& path);

}  // babelsim 命名空间
