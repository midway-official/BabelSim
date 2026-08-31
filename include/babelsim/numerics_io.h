#pragma once

#include "babelsim/config.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"
#include "babelsim/time.h"

#include <filesystem>

namespace babelsim {

// Case 的通用数值输入。它把时间控制、空间/时间离散与线性求解配置分开，Physics
// Case reader 只解析本领域的物性和算法控制。
Methods readMethodsFile(const std::filesystem::path& path);
TimeControl readTimeControlFile(const std::filesystem::path& path);
void readLinearSolverLine(
    const std::filesystem::path& path,
    const ConfigLine& line,
    LinearSolverConfig& result);

}  // babelsim 命名空间
