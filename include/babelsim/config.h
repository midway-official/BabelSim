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

// Read-only metadata used by diagnostics and configuration reports.  Querying
// this state never marks an entry as consumed.
struct ParameterInfo {
    bool configured = false;
    bool consumed = false;
    std::size_t line = 0;
};

// token 形式只用于嵌套值能显著提高可读性的场合（当前为 Field 文件），不是通用脚本层。
std::vector<ConfigToken> readConfigTokens(const std::filesystem::path& path);

// 命名字典的只读数值接口。求解器只声明自己需要的物性/算法参数，不再编写解析器。
// 记录已读取项，运行前拒绝拼错或无人使用的设置；错误包含文件和行号。
class Parameters {
public:
    explicit Parameters(const std::filesystem::path& path);
    bool contains(const std::string& key) const;
    const std::filesystem::path& sourcePath() const { return m_path; }
    ParameterInfo inspect(const std::string& key) const;
    std::string word(const std::string& key) const;
    std::string word(const std::string& key, const std::string& fallback) const;
    bool boolean(const std::string& key) const;
    bool boolean(const std::string& key, bool fallback) const;
    double number(const std::string& key) const;
    double number(const std::string& key, double fallback) const;
    double positive(const std::string& key) const;
    double positive(const std::string& key, double fallback) const;
    double fraction(const std::string& key, double fallback) const;
    double nonnegative(const std::string& key) const;
    double nonnegative(const std::string& key, double fallback) const;
    int integer(const std::string& key, int fallback) const;
    int integer(const std::string& key) const;
    int integer(const std::string& key, int fallback, int minimum, int maximum) const;
    void requireAllUsed() const;

    // 字典基础设施入口，普通 Solver 不需要访问 token。
    const ConfigLine& entry(const std::string& key) const;

private:
    std::filesystem::path m_path;
    std::vector<ConfigLine> m_lines;
    mutable std::vector<bool> m_used;
    [[noreturn]] void invalid(const ConfigLine& line, const std::string& message) const;
};

}  // babelsim 命名空间
