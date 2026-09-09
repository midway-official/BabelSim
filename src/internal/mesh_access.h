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
    static Index globalFaceId(const Mesh& mesh, Index face) { return mesh.globalFaceId(face); }
    static Index cellOwnerRank(const Mesh& mesh, Index cell) { return mesh.cellOwnerRank(cell); }
    static Index faceOwnerRank(const Mesh& mesh, Index face) { return mesh.faceOwnerRank(face); }
    static void setPartition(
        Mesh& mesh,
        Index global_cells,
        Index layers,
        std::vector<Index> cell_ids,
        std::vector<Index> cell_owners,
        std::vector<Index> cell_depths,
        std::vector<Index> face_ids,
        std::vector<Index> face_owners,
        Index local_rank)
    {
        mesh.setPartition(global_cells, layers, std::move(cell_ids), std::move(cell_owners),
                          std::move(cell_depths), std::move(face_ids),
                          std::move(face_owners), local_rank);
    }
};

inline auto& meshData(Mesh& mesh) { return MeshAccess::data(mesh); }
inline const auto& meshData(const Mesh& mesh) { return MeshAccess::data(mesh); }
inline Index ownedCellCount(const Mesh& mesh) { return MeshAccess::ownedCellCount(mesh); }
inline bool isOwned(const Mesh& mesh, Index cell) { return MeshAccess::isOwned(mesh, cell); }
inline Index ownedIndex(const Mesh& mesh, Index cell) { return MeshAccess::ownedIndex(mesh, cell); }
inline Index globalCellId(const Mesh& mesh, Index cell) { return MeshAccess::globalCellId(mesh, cell); }
inline Index globalFaceId(const Mesh& mesh, Index face) { return MeshAccess::globalFaceId(mesh, face); }
inline Index cellOwnerRank(const Mesh& mesh, Index cell) { return MeshAccess::cellOwnerRank(mesh, cell); }
inline Index faceOwnerRank(const Mesh& mesh, Index face) { return MeshAccess::faceOwnerRank(mesh, face); }

}  // babelsim::detail 命名空间
