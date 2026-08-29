#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace babelsim {

// 供 BabelSim 行式字典复用的小型解析器。值以空白分隔，'#' 之后均为注释。
struct ConfigLine {
    std::size_t number = 0;
    std::vector<std::string> tokens;
};

std::vector<ConfigLine> readConfigLines(const std::filesystem::path& path);

struct ConfigToken {
    std::size_t line = 0;
    std::string text;
};

// token 形式只用于嵌套值能显著提高可读性的场合（当前为 Field 文件），不是通用脚本层。
std::vector<ConfigToken> readConfigTokens(const std::filesystem::path& path);

}  // babelsim 命名空间
