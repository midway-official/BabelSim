#pragma once

#include <filesystem>
#include <string>

namespace babelsim {

// 与求解器无关的案例入口。物理专属字典分离存放，使启动器能先选择求解器再读取它们。
struct CaseDefinition {
    std::filesystem::path root;
    std::string solver;
    std::filesystem::path mesh_file;
    std::filesystem::path fields_directory;
    std::filesystem::path physics_file;
    std::filesystem::path methods_file;
    std::filesystem::path solution_file;
    std::filesystem::path control_file;
    std::filesystem::path output_file;
};

struct OutputControl {
    std::filesystem::path directory = "results";
    std::string time_name = "final";
};

CaseDefinition readCase(const std::filesystem::path& case_directory);
OutputControl readOutputControl(const CaseDefinition& definition);
std::filesystem::path outputTimeDirectory(
    const CaseDefinition& definition,
    const std::string& requested_time = {});

}  // babelsim 命名空间
