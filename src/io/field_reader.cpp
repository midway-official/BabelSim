#include "internal/mesh_access.h"
#include "internal/field_access.h"
#include "babelsim/field_io.h"

#include "babelsim/config.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace babelsim {
namespace {

class Reader {
public:
    Reader(const std::filesystem::path& path, std::vector<ConfigToken> tokens)
        : m_path(path), m_tokens(std::move(tokens)) {}

    bool done() const { return m_position == m_tokens.size(); }
    bool nextIs(const char* text) const {
        return !done() && m_tokens[m_position].text == text;
    }
    const ConfigToken& take(const char* expected = nullptr) {
        if (done()) fail("unexpected end of file");
        const ConfigToken& token = m_tokens[m_position++];
        if (expected != nullptr && token.text != expected) {
            fail(token, std::string("expected ") + expected + ", got " + token.text);
        }
        return token;
    }
    [[noreturn]] void fail(const std::string& message) const {
        const std::size_t line = done() ? 0 : m_tokens[m_position].line;
        throw std::runtime_error("invalid " + m_path.string() + ":" +
                                 std::to_string(line) + ": " + message);
    }
    [[noreturn]] void fail(const ConfigToken& token, const std::string& message) const {
        throw std::runtime_error("invalid " + m_path.string() + ":" +
                                 std::to_string(token.line) + ": " + message);
    }
    double number() {
        const ConfigToken& token = take();
        try {
            std::size_t consumed = 0;
            const double result = std::stod(token.text, &consumed);
            if (consumed == token.text.size() && std::isfinite(result)) return result;
        } catch (const std::exception&) {
        }
        fail(token, "expected a number");
    }

private:
    const std::filesystem::path& m_path;
    std::vector<ConfigToken> m_tokens;
    std::size_t m_position = 0;
};

Index patchIndex(const Mesh& mesh, const std::string& name) {
    for (Index patch = 0; patch < static_cast<Index>(detail::meshData(mesh).patches.size()); ++patch) {
        if (detail::meshData(mesh).patches[static_cast<std::size_t>(patch)].name == name) return patch;
    }
    return invalid_index;
}

BoundaryType boundaryType(Reader& input) {
    const ConfigToken& token = input.take();
    if (token.text == "fixedValue" || token.text == "fixed_value" ||
        token.text == "dirichlet") return BoundaryType::FixedValue;
    if (token.text == "fixedGradient" || token.text == "fixed_gradient" ||
        token.text == "neumann") return BoundaryType::FixedGradient;
    if (token.text == "zeroGradient" || token.text == "zero_gradient") {
        return BoundaryType::ZeroGradient;
    }
    if (token.text == "symmetry" || token.text == "mirror") return BoundaryType::Symmetry;
    if (token.text == "inletOutlet" || token.text == "inlet_outlet") {
        return BoundaryType::InletOutlet;
    }
    input.fail(token, "unknown boundary condition " + token.text);
}

template <typename T>
T readValue(Reader& input);

template <>
double readValue<double>(Reader& input) {
    input.take("(");
    const double result = input.number();
    input.take(")");
    return result;
}

template <>
Vec3 readValue<Vec3>(Reader& input) {
    input.take("(");
    const Vec3 result{input.number(), input.number(), input.number()};
    input.take(")");
    return result;
}

template <>
Tensor3 readValue<Tensor3>(Reader& input) {
    input.take("(");
    Tensor3 result;
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column) result[row][column] = input.number();
    input.take(")");
    return result;
}

template <typename T>
bool readDataValue(std::istringstream& row, T& value);

template <>
bool readDataValue<double>(std::istringstream& row, double& value) {
    return static_cast<bool>(row >> value);
}

template <>
bool readDataValue<Vec3>(std::istringstream& row, Vec3& value) {
    return static_cast<bool>(row >> value.x >> value.y >> value.z);
}

template <>
bool readDataValue<Tensor3>(std::istringstream& row, Tensor3& value) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (!(row >> value[i][j])) return false;
    return true;
}

template <typename T>
constexpr int componentCount() {
    if constexpr (std::is_same_v<T, Vec3>) return 3;
    if constexpr (std::is_same_v<T, Tensor3>) return 9;
    return 1;
}

template <typename T>
void readNonuniformCellValues(const std::filesystem::path& data_path, Field<T>& field) {
    const Index global_cells = field.mesh().globalCellCount();
    if (global_cells <= 0) {
        throw std::runtime_error("nonuniform field requires a partitioned mesh with global cell IDs");
    }
    if (!std::filesystem::exists(data_path))
        throw std::runtime_error("cannot open nonuniform field data: " + data_path.string());
    std::vector<T> values(static_cast<std::size_t>(global_cells));
    std::vector<bool> seen(static_cast<std::size_t>(global_cells), false);
    const auto store = [&](Index global_id, const T& value,
                           const std::filesystem::path& source, std::size_t line_number) {
        if (global_id < 0 || global_id >= global_cells) {
            throw std::runtime_error("invalid global cell ID in " + source.string() + ":" +
                                     std::to_string(line_number));
        }
        const std::size_t index = static_cast<std::size_t>(global_id);
        if (seen[index]) {
            throw std::runtime_error("duplicate global cell ID in " + source.string() + ":" +
                                     std::to_string(line_number));
        }
        values[index] = value;
        seen[index] = true;
    };
    const auto read_csv = [&](const std::filesystem::path& csv_path) {
        std::ifstream csv(csv_path);
        if (!csv) throw std::runtime_error("cannot open restart field CSV: " + csv_path.string());
        std::string line;
        if (!std::getline(csv, line))
            throw std::runtime_error("restart field CSV is empty: " + csv_path.string());
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::ostringstream expected_header;
        expected_header << "global_id,x,y,z";
        for (int component = 0; component < componentCount<T>(); ++component)
            expected_header << ",value" << component;
        if (line != expected_header.str())
            throw std::runtime_error("unexpected restart field CSV header in " + csv_path.string());
        std::size_t line_number = 1;
        while (std::getline(csv, line)) {
            ++line_number;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream row(line);
            Index global_id = invalid_index;
            double x = 0.0, y = 0.0, z = 0.0;
            T value{};
            if (!(row >> global_id >> x >> y >> z) || !readDataValue(row, value))
                throw std::runtime_error("invalid restart field row in " + csv_path.string() + ":" +
                                         std::to_string(line_number));
            std::string trailing;
            if (row >> trailing)
                throw std::runtime_error("unexpected data after restart field value in " + csv_path.string() + ":" +
                                         std::to_string(line_number));
            store(global_id, value, csv_path, line_number);
        }
    };

    if (std::filesystem::is_directory(data_path)) {
        std::vector<std::filesystem::path> rank_files;
        const std::string filename = field.name() + ".csv";
        for (const auto& entry : std::filesystem::directory_iterator(data_path)) {
            if (!entry.is_directory() || entry.path().filename().string().rfind("rank-", 0) != 0) continue;
            const auto candidate = entry.path() / filename;
            if (std::filesystem::is_regular_file(candidate)) rank_files.push_back(candidate);
        }
        std::sort(rank_files.begin(), rank_files.end());
        if (rank_files.empty())
            throw std::runtime_error("no rank-sharded restart CSVs for field " + field.name() +
                                     " under " + data_path.string());
        for (const auto& csv_path : rank_files) read_csv(csv_path);
    } else if (data_path.extension() == ".csv") {
        read_csv(data_path);
    } else {
        std::ifstream data(data_path);
        if (!data) throw std::runtime_error("cannot open nonuniform field data: " + data_path.string());
        std::string text;
        std::size_t line_number = 0;
        while (std::getline(data, text)) {
            ++line_number;
            const std::size_t comment = text.find('#');
            if (comment != std::string::npos) text.resize(comment);
            std::istringstream row(text);
            Index global_id = invalid_index;
            if (!(row >> global_id)) continue;
            T value{};
            if (!readDataValue(row, value)) {
                throw std::runtime_error("invalid field value in " + data_path.string() + ":" +
                                         std::to_string(line_number));
            }
            std::string trailing;
            if (row >> trailing) {
                throw std::runtime_error("unexpected data after field value in " + data_path.string() + ":" +
                                         std::to_string(line_number));
            }
            store(global_id, value, data_path, line_number);
        }
    }
    if (std::find(seen.begin(), seen.end(), false) != seen.end()) {
        throw std::runtime_error("nonuniform field data is missing global cells: " + data_path.string());
    }
    auto* local = detail::fieldData(field);
    for (Index cell = 0; cell < field.mesh().cellCount(); ++cell) {
        const Index global_id = detail::globalCellId(field.mesh(), cell);
        local[static_cast<std::size_t>(cell)] = values[static_cast<std::size_t>(global_id)];
    }
}

template <typename T>
void read(const std::filesystem::path& path, Field<T>& field, const char* type_name) {
    if (field.location() != FieldLocation::Cell) {
        throw std::invalid_argument("field files currently initialize cell fields only");
    }
    Reader input(path, readConfigTokens(path));
    input.take("field");
    if (input.take().text != field.name()) input.fail("field name must be " + field.name());
    input.take("{");
    bool type = false;
    bool location = false;
    bool internal = false;
    std::vector<bool> configured(detail::meshData(field.mesh()).patches.size(), false);
    while (!input.nextIs("}")) {
        const ConfigToken& entry = input.take();
        if (entry.text == "type") {
            if (type || input.take().text != type_name) input.fail(entry, "field type is invalid");
            type = true;
        } else if (entry.text == "location") {
            if (location || input.take().text != "cell") input.fail(entry, "field location must be cell");
            location = true;
        } else if (entry.text == "internal") {
            if (internal) input.fail(entry, "internal value may only be configured once");
            const ConfigToken& mode = input.take();
            if (mode.text == "uniform") {
                field.fill(readValue<T>(input));
            } else if (mode.text == "file") {
                const ConfigToken& file = input.take();
                readNonuniformCellValues(path.parent_path() / file.text, field);
            } else {
                input.fail(mode, "internal value must be uniform or file-backed");
            }
            internal = true;
        } else if (entry.text == "boundary") {
            input.take("{");
            while (!input.nextIs("}")) {
                const ConfigToken& patch_name = input.take();
                const Index patch = patchIndex(field.mesh(), patch_name.text);
                if (patch == invalid_index || configured[static_cast<std::size_t>(patch)]) {
                    input.fail(patch_name, "unknown or duplicate boundary patch " + patch_name.text);
                }
                input.take("{");
                input.take("type");
                const BoundaryType kind = boundaryType(input);
                BoundaryCondition<T> condition;
                condition.type = kind;
                if (kind == BoundaryType::FixedValue || kind == BoundaryType::FixedGradient ||
                    kind == BoundaryType::InletOutlet) {
                    input.take("value");
                    condition.value = readValue<T>(input);
                }
                input.take("}");
                field.setBoundary(patch, condition);
                configured[static_cast<std::size_t>(patch)] = true;
            }
            input.take("}");
        } else {
            input.fail(entry, "unknown entry " + entry.text);
        }
    }
    input.take("}");
    if (!input.done()) input.fail("unexpected trailing tokens");
    // 分区网格会附加 processor patch；其值由 halo exchange 提供，输入文件无需
    // 为每个 rank 重复保存一份。物理边界仍必须显式配置。
    for (Index patch = 0;
         patch < static_cast<Index>(configured.size()); ++patch) {
        if (!configured[static_cast<std::size_t>(patch)] &&
            detail::meshData(field.mesh()).patches[static_cast<std::size_t>(patch)].kind ==
                PatchKind::Processor) {
            field.setBoundary(patch, BoundaryCondition<T>::zeroGradient());
            configured[static_cast<std::size_t>(patch)] = true;
        }
    }
    if (!type || !location || !internal ||
        std::find(configured.begin(), configured.end(), false) != configured.end()) {
        throw std::runtime_error("incomplete field file: " + path.string());
    }
}

}  // 匿名命名空间

void readFieldFile(const std::filesystem::path& path, ScalarField& field) {
    read(path, field, "scalar");
}

void readFieldFile(const std::filesystem::path& path, VectorField& field) {
    read(path, field, "vector");
}

void readFieldFile(const std::filesystem::path& path, TensorField& field) {
    read(path, field, "tensor");
}

}  // babelsim 命名空间
