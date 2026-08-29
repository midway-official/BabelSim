#include "babelsim/legacy_taiho.h"

#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace babelsim {
namespace {

template <typename T>
std::vector<T> readValues(
    const std::filesystem::path& path,
    std::size_t expected)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open Taiho mesh file: " + path.string());
    }
    std::vector<T> values;
    values.reserve(expected);
    T value{};
    while (input >> value) {
        values.push_back(value);
    }
    if (values.size() != expected) {
        std::ostringstream message;
        message << path << " contains " << values.size()
                << " values, expected " << expected;
        throw std::runtime_error(message.str());
    }
    return values;
}

struct PatchAndCondition {
    BoundaryPatch patch;
    ImportedPatchConditions condition;
};

PatchAndCondition convertPatch(
    int type,
    int zone,
    const std::vector<Vec3>& zone_velocity)
{
    if (zone < 0 || static_cast<std::size_t>(zone) >= zone_velocity.size()) {
        throw std::runtime_error("Taiho boundary zone is outside zoneuv.txt");
    }
    PatchAndCondition result;
    if (type > 0) {
        result.patch = {
            "wall_" + std::to_string(type) + '_' + std::to_string(zone),
            PatchKind::Wall,
            {},
        };
        result.condition.velocity =
            BoundaryCondition<Vec3>::fixedValue(zone_velocity[static_cast<std::size_t>(zone)]);
        result.condition.pressure = BoundaryCondition<double>::zeroGradient();
    } else if (type == -2) {
        result.patch = {
            "velocityInlet_" + std::to_string(zone),
            PatchKind::Inlet,
            {},
        };
        result.condition.velocity =
            BoundaryCondition<Vec3>::fixedValue(zone_velocity[static_cast<std::size_t>(zone)]);
        result.condition.pressure = BoundaryCondition<double>::zeroGradient();
    } else if (type == -1) {
        result.patch = {
            "pressureOutlet_" + std::to_string(zone),
            PatchKind::Outlet,
            {},
        };
        result.condition.velocity = BoundaryCondition<Vec3>::zeroGradient();
        result.condition.pressure = BoundaryCondition<double>::fixedValue(0.0);
    } else {
        throw std::runtime_error(
            "unsupported Taiho physical boundary type: " + std::to_string(type));
    }
    return result;
}

}  // namespace

ImportedTaihoMesh readTaihoMesh(const std::filesystem::path& directory) {
    if (!std::filesystem::is_directory(directory)) {
        throw std::runtime_error("Taiho mesh directory does not exist: " + directory.string());
    }
    std::ifstream parameters(directory / "params.txt");
    Index old_nx = 0;
    Index old_ny = 0;
    if (!(parameters >> old_nx >> old_ny) || old_nx < 3 || old_ny < 3) {
        throw std::runtime_error("Taiho params.txt must contain nx ny >= 3");
    }
    const Index nx = old_nx - 2;
    const Index ny = old_ny - 2;
    const std::size_t old_cells =
        static_cast<std::size_t>(old_nx) * static_cast<std::size_t>(old_ny);
    const std::size_t old_vertices =
        static_cast<std::size_t>(old_nx + 1) * static_cast<std::size_t>(old_ny + 1);
    const auto x = readValues<double>(directory / "x.dat", old_vertices);
    const auto y = readValues<double>(directory / "y.dat", old_vertices);
    const auto type = readValues<int>(directory / "bctype.dat", old_cells);
    const auto zone = readValues<int>(directory / "zoneid.dat", old_cells);

    std::ifstream zone_file(directory / "zoneuv.txt");
    if (!zone_file) {
        throw std::runtime_error("cannot open Taiho zoneuv.txt");
    }
    std::vector<Vec3> zone_velocity;
    double u = 0.0;
    double v = 0.0;
    while (zone_file >> u >> v) {
        zone_velocity.push_back({u, v, 0.0});
    }
    if (zone_velocity.empty()) {
        throw std::runtime_error("Taiho zoneuv.txt is empty");
    }

    const auto oldVertex = [old_nx](Index i, Index j) {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(old_nx + 1) +
            static_cast<std::size_t>(i);
    };
    std::vector<Vec3> points;
    points.reserve(
        static_cast<std::size_t>(nx + 1) * static_cast<std::size_t>(ny + 1) * 2U);
    for (Index k = 0; k <= 1; ++k) {
        for (Index j = 0; j <= ny; ++j) {
            for (Index i = 0; i <= nx; ++i) {
                const std::size_t old = oldVertex(i + 1, j + 1);
                points.push_back({x[old], y[old], static_cast<double>(k)});
            }
        }
    }

    auto temporary_patches = defaultPatches();
    temporary_patches[static_cast<std::size_t>(Side::ZMin)] =
        {"front", PatchKind::Symmetry};
    temporary_patches[static_cast<std::size_t>(Side::ZMax)] =
        {"back", PatchKind::Symmetry};
    Mesh mesh = Mesh::structured({nx, ny, 1}, std::move(points), temporary_patches);

    std::map<std::pair<int, int>, Index> imported_patch_ids;
    std::vector<BoundaryPatch> patches;
    std::vector<ImportedPatchConditions> conditions;
    std::vector<Index> new_face_patch(
        static_cast<std::size_t>(mesh.faceCount()), invalid_index);
    const auto oldCell = [old_nx](Index i, Index j) {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(old_nx) +
            static_cast<std::size_t>(i);
    };

    for (Index face = 0; face < mesh.faceCount(); ++face) {
        if (!mesh.boundaryFace(face)) {
            continue;
        }
        const auto f = static_cast<std::size_t>(face);
        const Index old_side = mesh.face_patch[f];
        if (old_side == static_cast<Index>(Side::ZMin) ||
            old_side == static_cast<Index>(Side::ZMax)) {
            continue;
        }
        const Index owner = mesh.face_owner[f];
        const Index i = owner % nx;
        const Index j = (owner / nx) % ny;
        std::size_t boundary_cell = 0;
        switch (static_cast<Side>(old_side)) {
            case Side::XMin:
                boundary_cell = oldCell(0, j + 1);
                break;
            case Side::XMax:
                boundary_cell = oldCell(old_nx - 1, j + 1);
                break;
            case Side::YMin:
                boundary_cell = oldCell(i + 1, 0);
                break;
            case Side::YMax:
                boundary_cell = oldCell(i + 1, old_ny - 1);
                break;
            case Side::ZMin:
            case Side::ZMax:
                throw std::logic_error("unexpected extrusion side");
        }
        const std::pair<int, int> key{type[boundary_cell], zone[boundary_cell]};
        auto [position, inserted] = imported_patch_ids.emplace(
            key, static_cast<Index>(patches.size()));
        if (inserted) {
            PatchAndCondition converted = convertPatch(
                key.first, key.second, zone_velocity);
            patches.push_back(std::move(converted.patch));
            conditions.push_back(std::move(converted.condition));
        }
        new_face_patch[f] = position->second;
        patches[static_cast<std::size_t>(position->second)].faces.push_back(face);
    }

    for (Side side : {Side::ZMin, Side::ZMax}) {
        const Index patch_id = static_cast<Index>(patches.size());
        BoundaryPatch patch{
            side == Side::ZMin ? "front" : "back",
            PatchKind::Symmetry,
            {},
        };
        for (Index face : mesh.patches[static_cast<std::size_t>(side)].faces) {
            patch.faces.push_back(face);
            new_face_patch[static_cast<std::size_t>(face)] = patch_id;
        }
        patches.push_back(std::move(patch));
        conditions.push_back({
            BoundaryCondition<Vec3>::symmetry(),
            BoundaryCondition<double>::symmetry(),
        });
    }

    mesh.patches = std::move(patches);
    mesh.face_patch = std::move(new_face_patch);
    mesh.validate();
    return {std::move(mesh), std::move(conditions)};
}

void applyImportedBoundaryConditions(
    const ImportedTaihoMesh& imported,
    VectorField& velocity,
    ScalarField& pressure)
{
    if (&velocity.mesh() != &imported.mesh ||
        &pressure.mesh() != &imported.mesh ||
        velocity.location() != FieldLocation::Cell ||
        pressure.location() != FieldLocation::Cell ||
        imported.conditions.size() != imported.mesh.patches.size()) {
        throw std::invalid_argument("imported boundary fields do not match the mesh");
    }
    for (Index patch = 0;
         patch < static_cast<Index>(imported.conditions.size()); ++patch) {
        velocity.setBoundary(
            patch, imported.conditions[static_cast<std::size_t>(patch)].velocity);
        pressure.setBoundary(
            patch, imported.conditions[static_cast<std::size_t>(patch)].pressure);
    }
}

}  // namespace babelsim

