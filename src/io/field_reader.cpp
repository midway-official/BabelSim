#include "internal/mesh_access.h"
#include "babelsim/field_io.h"

#include "babelsim/config.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
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
            if (internal || input.take().text != "uniform") {
                input.fail(entry, "internal value must be uniform");
            }
            field.fill(readValue<T>(input));
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
