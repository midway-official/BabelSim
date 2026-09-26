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
    if (magic != "BABELSIM_MESH" || version != 3) {
        invalidFile(path, "only BABELSIM_MESH version 3 is supported");
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

    if (read<std::string>(input, path, "faces keyword") != "faces") {
        invalidFile(path, "expected faces");
    }
    const Index face_count = count(input, path, "face count");
    std::vector<PolyhedralFaceSpec> faces;
    faces.reserve(static_cast<std::size_t>(face_count));
    for (Index face = 0; face < face_count; ++face) {
        if (read<std::string>(input, path, "face keyword") != "face") {
            invalidFile(path, "expected face");
        }
        const Index vertex_count_face = count(input, path, "face vertex count");
        PolyhedralFaceSpec specification;
        specification.owner = read<Index>(input, path, "face owner");
        specification.neighbour = read<Index>(input, path, "face neighbour");
        specification.patch = read<Index>(input, path, "face patch");
        specification.vertices.resize(static_cast<std::size_t>(vertex_count_face));
        for (Index& vertex : specification.vertices) {
            vertex = read<Index>(input, path, "face vertex index");
        }
        faces.push_back(std::move(specification));
    }

    if (read<std::string>(input, path, "patches keyword") != "patches") {
        invalidFile(path, "expected patches");
    }
    const Index patch_count = count(input, path, "patch count");
    std::vector<PatchSpec> patches;
    patches.reserve(static_cast<std::size_t>(patch_count));
    for (Index patch = 0; patch < patch_count; ++patch) {
        if (read<std::string>(input, path, "patch keyword") != "patch") {
            invalidFile(path, "expected patch");
        }
        const std::string name = read<std::string>(input, path, "patch name");
        const PatchKind kind = patchKind(read<std::string>(input, path, "patch kind"), path);
        if (name.empty()) invalidFile(path, "patch name must not be empty");
        patches.push_back({name, kind});
    }
    const std::string end = read<std::string>(input, path, "end marker");
    if (end != "end") invalidFile(path, "expected end");
    std::string trailing;
    if (input >> trailing) invalidFile(path, "unexpected trailing token " + trailing);

    try {
        return Mesh::polyhedral(std::move(vertices), std::move(faces), std::move(patches));
    } catch (const std::exception& error) {
        invalidFile(path, error.what());
    }
}

}  // namespace babelsim
