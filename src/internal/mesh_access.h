#pragma once

#include "babelsim/mesh.h"

namespace babelsim::detail {

// 维护接口：仅网格构造、分区和离散实现可访问缓存/映射；普通 Solver 只读几何。
struct MeshAccess {
    static auto& data(Mesh& mesh) { return mesh.m_storage; }
    static const auto& data(const Mesh& mesh) { return mesh.m_storage; }
    static Index ownedCellCount(const Mesh& mesh) { return mesh.ownedCellCount(); }
    static bool isOwned(const Mesh& mesh, Index cell) { return mesh.isOwned(cell); }
    static Index ownedIndex(const Mesh& mesh, Index cell) { return mesh.ownedIndex(cell); }
    static Index globalCellId(const Mesh& mesh, Index cell) { return mesh.globalCellId(cell); }
    static void replace(Mesh& target, Mesh source) { target = std::move(source); }
    static void setOwnership(Mesh& mesh, std::array<Index, 3> global_cells,
                             Index offset, Index begin, Index end, Index layers) {
        mesh.setOwnership(global_cells, offset, begin, end, layers);
    }
    static void setPatches(Mesh& mesh, std::vector<BoundaryPatch> patches) {
        mesh.setPatches(std::move(patches));
    }
    static void addPatchFace(Mesh& mesh, Index patch, Index face) {
        mesh.addPatchFace(patch, face);
    }
};

inline auto& meshData(Mesh& mesh) { return MeshAccess::data(mesh); }
inline const auto& meshData(const Mesh& mesh) { return MeshAccess::data(mesh); }
inline Index ownedCellCount(const Mesh& mesh) { return MeshAccess::ownedCellCount(mesh); }
inline bool isOwned(const Mesh& mesh, Index cell) { return MeshAccess::isOwned(mesh, cell); }
inline Index ownedIndex(const Mesh& mesh, Index cell) { return MeshAccess::ownedIndex(mesh, cell); }
inline Index globalCellId(const Mesh& mesh, Index cell) { return MeshAccess::globalCellId(mesh, cell); }

}  // babelsim::detail 命名空间
