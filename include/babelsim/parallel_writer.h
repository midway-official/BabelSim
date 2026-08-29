#pragma once

#include "babelsim/parallel.h"

#include <filesystem>
#include <string>
#include <vector>

namespace babelsim {

struct FieldOutputInfo {
    std::string name;
    std::string type;
    FieldLocation location = FieldLocation::Cell;
};

// 每个进程仅在 <time>/rank-0000/ 下写出自己拥有的实体；独立后处理器按全局 ID
// 合并这些文件。
void writeOwnedFieldCsv(
    const std::filesystem::path& time_directory,
    const ScalarField& field,
    const ParallelContext& parallel);
void writeOwnedFieldCsv(
    const std::filesystem::path& time_directory,
    const VectorField& field,
    const ParallelContext& parallel);
void writeOwnedFieldCsv(
    const std::filesystem::path& time_directory,
    const TensorField& field,
    const ParallelContext& parallel);

void writeOwnedResultMetadata(
    const std::filesystem::path& time_directory,
    const Mesh& mesh,
    const ParallelContext& parallel,
    const std::string& time_name,
    const std::vector<FieldOutputInfo>& fields);

}  // babelsim 命名空间
