#include "babelsim/mesh_io.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace babelsim {
namespace {

[[noreturn]] void invalidFile(const std::filesystem::path& path, const std::string& message) {
    throw std::runtime_error("invalid BabelSim mesh file " + path.string() + ": " + message);
}

template <typename T>
T read(std::istream& input, const std::filesystem::path& path, const char* what) {
    T value{};
    if (!(input >> value)) invalidFile(path, std::string("missing ") + what);
    return value;
}

Index count(std::istream& input, const std::filesystem::path& path, const char* what) {
    const Index value = read<Index>(input, path, what);
    if (value <= 0 || value == std::numeric_limits<Index>::max()) {
        invalidFile(path, std::string(what) + " must be positive");
    }
    return value;
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
    if (!input) throw std::runtime_error("cannot open BabelSim mesh file: " + path.string());

    const std::string magic = read<std::string>(input, path, "file header");
    const int version = read<int>(input, path, "format version");
    if (magic != "BABELSIM_MESH" || version != 2) {
        invalidFile(path, "only BABELSIM_MESH version 2 is supported");
    }
    if (read<std::string>(input, path, "vertices keyword") != "vertices") {
        invalidFile(path, "expected vertices");
    }
    const Index vertex_count = count(input, path, "vertex count");
    std::vector<Vec3> vertices(static_cast<std::size_t>(vertex_count));
    for (Vec3& vertex : vertices) {
        vertex = {read<double>(input, path, "vertex x"), read<double>(input, path, "vertex y"),
                  read<double>(input, path, "vertex z")};
        if (!isFinite(vertex)) invalidFile(path, "vertex must be finite");
    }

    if (read<std::string>(input, path, "cells keyword") != "cells") {
        invalidFile(path, "expected cells");
    }
    const Index cell_count = count(input, path, "cell count");
    std::vector<std::array<Index, 8>> cells(static_cast<std::size_t>(cell_count));
    for (auto& cell : cells) {
        for (Index& vertex : cell) vertex = read<Index>(input, path, "cell vertex index");
    }

    if (read<std::string>(input, path, "patches keyword") != "patches") {
        invalidFile(path, "expected patches");
    }
    const Index patch_count = count(input, path, "patch count");
    std::vector<PatchSpec> patches;
    std::vector<BoundaryFaceSpec> boundary_faces;
    patches.reserve(static_cast<std::size_t>(patch_count));
    for (Index patch = 0; patch < patch_count; ++patch) {
        if (read<std::string>(input, path, "patch keyword") != "patch") {
            invalidFile(path, "expected patch");
        }
        const std::string name = read<std::string>(input, path, "patch name");
        const PatchKind kind = patchKind(read<std::string>(input, path, "patch kind"), path);
        const Index face_count = count(input, path, "patch face count");
        if (name.empty()) invalidFile(path, "patch name must not be empty");
        patches.push_back({name, kind});
        for (Index face = 0; face < face_count; ++face) {
            BoundaryFaceSpec boundary;
            boundary.patch = patch;
            for (Index& vertex : boundary.vertices) {
                vertex = read<Index>(input, path, "boundary face vertex index");
            }
            boundary_faces.push_back(boundary);
        }
    }
    const std::string end = read<std::string>(input, path, "end marker");
    if (end != "end") invalidFile(path, "expected end");
    std::string trailing;
    if (input >> trailing) invalidFile(path, "unexpected trailing token " + trailing);

    try {
        return Mesh::unstructured(std::move(vertices), std::move(cells), std::move(patches),
                                  std::move(boundary_faces));
    } catch (const std::exception& error) {
        invalidFile(path, error.what());
    }
}

}  // namespace babelsim
