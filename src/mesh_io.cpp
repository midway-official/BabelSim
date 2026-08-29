#include "babelsim/mesh_io.h"

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace babelsim {
namespace {

[[noreturn]] void invalidFile(
    const std::filesystem::path& path,
    const std::string& message)
{
    throw std::runtime_error("invalid BabelSim mesh file " + path.string() +
                             ": " + message);
}

template <typename T>
T read(std::istream& input, const std::filesystem::path& path, const char* what) {
    T value{};
    if (!(input >> value)) {
        invalidFile(path, std::string("missing ") + what);
    }
    return value;
}

Index sideIndex(const std::string& name, const std::filesystem::path& path) {
    static constexpr std::array<const char*, 6> names = {
        "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
    for (Index side = 0; side < static_cast<Index>(names.size()); ++side) {
        if (name == names[static_cast<std::size_t>(side)]) {
            return side;
        }
    }
    invalidFile(path, "unknown logical side " + name);
}

PatchKind patchKind(const std::string& name, const std::filesystem::path& path) {
    if (name == "generic") return PatchKind::Generic;
    if (name == "wall") return PatchKind::Wall;
    if (name == "inlet") return PatchKind::Inlet;
    if (name == "outlet") return PatchKind::Outlet;
    if (name == "symmetry" || name == "mirror") return PatchKind::Symmetry;
    if (name == "processor") return PatchKind::Processor;
    invalidFile(path, "unknown patch kind " + name);
}

}  // namespace

Mesh readMeshFile(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open BabelSim mesh file: " + path.string());
    }

    const std::string magic = read<std::string>(input, path, "file header");
    const int version = read<int>(input, path, "format version");
    if (magic != "BABELSIM_MESH" || version != 1) {
        invalidFile(path, "expected BABELSIM_MESH version 1");
    }

    if (read<std::string>(input, path, "dimensions keyword") != "dimensions") {
        invalidFile(path, "expected dimensions");
    }
    std::array<Index, 3> dimensions = {
        read<Index>(input, path, "nx"),
        read<Index>(input, path, "ny"),
        read<Index>(input, path, "nz")};
    if (dimensions[0] <= 0 || dimensions[1] <= 0 || dimensions[2] <= 0) {
        invalidFile(path, "dimensions must be positive");
    }

    if (read<std::string>(input, path, "geometry keyword") != "geometry") {
        invalidFile(path, "expected geometry");
    }
    const std::string geometry = read<std::string>(input, path, "geometry type");
    std::vector<Vec3> vertices;
    Vec3 minimum{};
    Vec3 maximum{};
    if (geometry == "cartesian") {
        if (read<std::string>(input, path, "bounds keyword") != "bounds") {
            invalidFile(path, "expected bounds after cartesian geometry");
        }
        minimum = {
            read<double>(input, path, "minimum x"),
            read<double>(input, path, "minimum y"),
            read<double>(input, path, "minimum z")};
        maximum = {
            read<double>(input, path, "maximum x"),
            read<double>(input, path, "maximum y"),
            read<double>(input, path, "maximum z")};
    } else if (geometry == "vertices") {
        const std::size_t count =
            static_cast<std::size_t>(dimensions[0] + 1) *
            static_cast<std::size_t>(dimensions[1] + 1) *
            static_cast<std::size_t>(dimensions[2] + 1);
        vertices.resize(count);
        for (Vec3& vertex : vertices) {
            vertex = {
                read<double>(input, path, "vertex x"),
                read<double>(input, path, "vertex y"),
                read<double>(input, path, "vertex z")};
        }
    } else {
        invalidFile(path, "geometry must be cartesian or vertices");
    }

    std::array<PatchSpec, 6> patches{};
    std::array<bool, 6> seen{};
    for (Index record = 0; record < 6; ++record) {
        if (read<std::string>(input, path, "patch keyword") != "patch") {
            invalidFile(path, "expected six patch records");
        }
        const std::string side_name = read<std::string>(input, path, "patch side");
        const std::string name = read<std::string>(input, path, "patch name");
        const std::string kind = read<std::string>(input, path, "patch kind");
        const Index side = sideIndex(side_name, path);
        if (seen[static_cast<std::size_t>(side)] || name.empty()) {
            invalidFile(path, "duplicate side or empty patch name");
        }
        seen[static_cast<std::size_t>(side)] = true;
        patches[static_cast<std::size_t>(side)] = {name, patchKind(kind, path)};
    }

    Mesh mesh = geometry == "cartesian"
        ? Mesh::cartesian(dimensions, minimum, maximum, patches)
        : Mesh::structured(dimensions, std::move(vertices), patches);
    std::string trailing;
    if (input >> trailing && trailing != "end") {
        invalidFile(path, "unexpected trailing token " + trailing);
    }
    return mesh;
}

}  // namespace babelsim
