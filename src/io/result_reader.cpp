#include "babelsim/result_reader.h"

#include "babelsim/config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace babelsim {
namespace {

struct Metadata {
    std::string time_name;
    int ranks = 0;
    std::array<Index, 3> dimensions{};
    std::vector<FieldOutputInfo> fields;
};

[[noreturn]] void invalid(const std::filesystem::path& path, const std::string& message) {
    throw std::runtime_error("invalid BabelSim result " + path.string() + ": " + message);
}

int components(const std::string& type, const std::filesystem::path& path) {
    if (type == "scalar") return 1;
    if (type == "vector") return 3;
    if (type == "tensor") return 9;
    invalid(path, "unknown field type " + type);
}

FieldLocation location(const std::string& name, const std::filesystem::path& path) {
    if (name == "cell") return FieldLocation::Cell;
    if (name == "face") return FieldLocation::Face;
    if (name == "vertex") return FieldLocation::Vertex;
    invalid(path, "unknown field location " + name);
}

int integer(const std::string& value, const std::filesystem::path& path) {
    try {
        std::size_t consumed = 0;
        const int result = std::stoi(value, &consumed);
        if (consumed == value.size()) return result;
    } catch (const std::exception&) {
    }
    invalid(path, "expected integer " + value);
}

Metadata readMetadata(const std::filesystem::path& path) {
    Metadata result;
    bool format = false;
    bool time = false;
    bool ranks = false;
    bool dimensions = false;
    for (const ConfigLine& line : readConfigLines(path)) {
        const std::string& key = line.tokens.front();
        if (key == "format" && line.tokens.size() == 3 && !format &&
            line.tokens[1] == "babelsim_result" && line.tokens[2] == "1") {
            format = true;
        } else if (key == "time" && line.tokens.size() == 2 && !time) {
            result.time_name = line.tokens[1];
            time = true;
        } else if (key == "ranks" && line.tokens.size() == 2 && !ranks) {
            result.ranks = integer(line.tokens[1], path);
            ranks = true;
        } else if (key == "global_dimensions" && line.tokens.size() == 4 && !dimensions) {
            result.dimensions = {
                integer(line.tokens[1], path), integer(line.tokens[2], path),
                integer(line.tokens[3], path)};
            dimensions = true;
        } else if (key == "field" && line.tokens.size() == 4) {
            const FieldOutputInfo info{line.tokens[1], line.tokens[2], location(line.tokens[3], path)};
            if (info.location != FieldLocation::Cell || components(info.type, path) == 0) {
                invalid(path, "only cell field results are currently supported");
            }
            if (std::find_if(result.fields.begin(), result.fields.end(), [&](const FieldOutputInfo& other) {
                    return other.name == info.name;
                }) != result.fields.end()) {
                invalid(path, "duplicate field " + info.name);
            }
            result.fields.push_back(info);
        } else if (key != "rank" && key != "owned_cells") {
            invalid(path, "invalid metadata record at line " + std::to_string(line.number));
        }
    }
    if (!format || !time || !ranks || !dimensions || result.ranks <= 0 || result.fields.empty()) {
        invalid(path, "metadata is incomplete");
    }
    return result;
}

std::vector<std::string> csv(const std::string& line) {
    std::vector<std::string> values;
    std::stringstream input(line);
    for (std::string value; std::getline(input, value, ',');) values.push_back(std::move(value));
    return values;
}

double number(const std::string& value, const std::filesystem::path& path) {
    try {
        std::size_t consumed = 0;
        const double result = std::stod(value, &consumed);
        if (consumed == value.size() && std::isfinite(result)) return result;
    } catch (const std::exception&) {
    }
    invalid(path, "expected finite number " + value);
}

bool sameFields(const std::vector<FieldOutputInfo>& a, const std::vector<FieldOutputInfo>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t field = 0; field < a.size(); ++field) {
        if (a[field].name != b[field].name || a[field].type != b[field].type ||
            a[field].location != b[field].location) return false;
    }
    return true;
}

}  // 匿名命名空间

ResultData readParallelResults(
    const std::filesystem::path& time_directory,
    Index global_cell_count)
{
    if (global_cell_count <= 0 || !std::filesystem::is_directory(time_directory)) {
        throw std::invalid_argument("result directory or global cell count is invalid");
    }
    std::vector<std::filesystem::path> rank_directories;
    for (const auto& entry : std::filesystem::directory_iterator(time_directory)) {
        if (entry.is_directory() && entry.path().filename().string().rfind("rank-", 0) == 0) {
            rank_directories.push_back(entry.path());
        }
    }
    std::sort(rank_directories.begin(), rank_directories.end());
    if (rank_directories.empty()) invalid(time_directory, "contains no rank directories");

    const Metadata first = readMetadata(rank_directories.front() / "metadata.bs");
    if (static_cast<int>(rank_directories.size()) != first.ranks) {
        invalid(time_directory, "rank directory count does not match metadata");
    }
    ResultData result;
    result.time_name = first.time_name;
    result.global_dimensions = first.dimensions;
    result.fields.reserve(first.fields.size());
    std::vector<std::vector<bool>> seen;
    for (const FieldOutputInfo& info : first.fields) {
        const int count = components(info.type, time_directory);
        result.fields.push_back({info, count,
            std::vector<double>(static_cast<std::size_t>(global_cell_count) * count)});
        seen.emplace_back(static_cast<std::size_t>(global_cell_count), false);
    }

    for (const auto& rank_directory : rank_directories) {
        const Metadata metadata = readMetadata(rank_directory / "metadata.bs");
        if (metadata.time_name != first.time_name || metadata.ranks != first.ranks ||
            metadata.dimensions != first.dimensions || !sameFields(metadata.fields, first.fields)) {
            invalid(rank_directory, "metadata does not match the other ranks");
        }
        for (std::size_t field = 0; field < result.fields.size(); ++field) {
            const ResultField& descriptor = result.fields[field];
            const auto path = rank_directory / (descriptor.info.name + ".csv");
            std::ifstream input(path);
            if (!input) invalid(path, "cannot open field data");
            std::string line;
            if (!std::getline(input, line)) invalid(path, "missing CSV header");
            const std::vector<std::string> header = csv(line);
            if (header.size() != static_cast<std::size_t>(4 + descriptor.components) ||
                header[0] != "global_id" || header[1] != "x" || header[2] != "y" || header[3] != "z") {
                invalid(path, "unexpected CSV header");
            }
            while (std::getline(input, line)) {
                const std::vector<std::string> values = csv(line);
                if (values.size() != header.size()) invalid(path, "malformed CSV record");
                const int id = integer(values[0], path);
                if (id < 0 || id >= global_cell_count || seen[field][static_cast<std::size_t>(id)]) {
                    invalid(path, "duplicate or out-of-range global cell id");
                }
                number(values[1], path);
                number(values[2], path);
                number(values[3], path);
                const std::size_t begin = static_cast<std::size_t>(id) * descriptor.components;
                for (int component = 0; component < descriptor.components; ++component) {
                    result.fields[field].values[begin + component] = number(values[4 + component], path);
                }
                seen[field][static_cast<std::size_t>(id)] = true;
            }
        }
    }
    for (const auto& field_seen : seen) {
        if (std::find(field_seen.begin(), field_seen.end(), false) != field_seen.end()) {
            invalid(time_directory, "rank files do not cover every global cell exactly once");
        }
    }
    return result;
}

}  // babelsim 命名空间
