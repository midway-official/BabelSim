#include "babelsim/parallel.h"
#include "babelsim/mesh_io.h"
#include "internal/mesh_access.h"
#include "internal/petsc_session.h"
#include "test_util.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <unistd.h>

using namespace babelsim;

namespace {
void verifyDistributedRead(const ParallelContext& parallel) {
    const auto global = readMeshFile("tests/data/babelsim_channel.mesh");
    const auto local = readDistributedMesh(parallel.rank == 0
        ? "tests/data/babelsim_channel.mesh" : "/not-a-worker-mesh-file", parallel);
    require(parallel.sum(detail::ownedCellCount(local)) == global.cellCount(),
            "distributed read lost owned cells");
    for (Index cell = 0; cell < local.cellCount(); ++cell) {
        const Index id = detail::globalCellId(local, cell);
        require(near(local.cellCentre(cell), global.cellCentre(id), 1e-13) &&
                near(local.cellVolume(cell), global.cellVolume(id), 1e-13),
                "distributed local geometry differs from the full mesh");
    }
    for (Index face = 0; face < local.faceCount(); ++face) {
        const Index id = detail::globalFaceId(local, face);
        const double orientation = detail::globalCellId(local, local.owner(face)) == global.owner(id) ? 1 : -1;
        require(near(local.faceCentre(face), global.faceCentre(id), 1e-13) &&
                near(local.faceAreaVector(face), orientation * global.faceAreaVector(id), 1e-13) &&
                local.facePoints(face).size() == global.facePoints(id).size(),
                "distributed local face lost its geometry or global orientation");
        if (global.neighbour(id) == invalid_index)
            require(local.patchName(local.boundaryPatch(face)) == global.patchName(global.boundaryPatch(id)),
                    "distributed read changed a physical patch");
    }
    HaloExchange halo(local, parallel);
    std::vector<double> values(static_cast<std::size_t>(local.cellCount()), -1.0);
    for (Index cell : detail::meshData(local).owned_cells) values[cell] = detail::globalCellId(local, cell) + 0.5;
    halo.exchange(values);
    for (Index cell = 0; cell < local.cellCount(); ++cell)
        require(values[cell] == detail::globalCellId(local, cell) + 0.5,
                "distributed read lost the authoritative halo mapping");

    bool caught = false;
    try { const auto missing = readDistributedMesh("/not-a-root-mesh-file", parallel); (void)missing; }
    catch (const std::runtime_error&) { caught = true; }
    require(parallel.sum(caught ? 1 : 0) == parallel.size,
            "distributed read did not report a root failure to every rank");

    // Independent concave geometry, disconnected components, and six-vertex
    // rings exercise payload lengths and quota balancing beyond hex meshes.
    const auto path = std::filesystem::path("/tmp") /
        ("babelsim-distribution-" + std::to_string(getpid()) + ".mesh");
    if (parallel.rank == 0) {
        const auto cell = readMeshFile("tests/data/babelsim_concave.mesh");
        std::ofstream out(path);
        out << "BABELSIM_MESH 3\nvertices\n" << 9 * cell.vertexCount() << '\n';
        for (Index copy = 0; copy < 9; ++copy)
            for (Index vertex = 0; vertex < cell.vertexCount(); ++vertex) {
                const auto point = cell.vertex(vertex);
                out << point.x + 4 * copy << ' ' << point.y << ' ' << point.z << '\n';
            }
        out << "faces\n" << 9 * cell.faceCount() << '\n';
        for (Index copy = 0; copy < 9; ++copy)
            for (Index face = 0; face < cell.faceCount(); ++face) {
                out << "face " << cell.facePoints(face).size() << ' ' << copy << " -1 0";
                for (Index vertex : cell.facePoints(face)) out << ' ' << vertex + copy * cell.vertexCount();
                out << '\n';
            }
        out << "patches\n1\npatch\nboundary\ngeneric\nend\n";
        require(out.good(), "cannot write distributed polygon fixture");
    }
    const auto polygons = readDistributedMesh(path, parallel);
    require(parallel.sum(detail::ownedCellCount(polygons)) == 9, "disconnected partition lost cells");
    for (Index cell = 0; cell < polygons.cellCount(); ++cell) {
        const Index id = detail::globalCellId(polygons, cell);
        require(near(polygons.cellVolume(cell), 3.0) &&
                near(polygons.cellCentre(cell), Vec3{4.0 * id + 5.0 / 6.0, 5.0 / 6.0, 0.5}),
                "distributed polygon reconstruction changed the independent geometry");
    }
    if (parallel.rank == 0) std::filesystem::remove(path);
}
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    const auto parallel = ParallelContext::world();
    try {
        require(parallel.size >= 2 && parallel.size <= 8,
                "partition test requires two to eight ranks");
        verifyDistributedRead(parallel);
        for (const auto shape : {std::array<Index, 3>{64, 64, 1},
                                 std::array<Index, 3>{17, 11, 3},
                                 std::array<Index, 3>{128, 2, 1}}) {
            const auto global = makeHexBox({shape[0], shape[1], shape[2]}, {0, 0, 0}, {1, 1, 0.01});
            const auto local = decompose(global, parallel);
            bool rejected_local = false;
            try { (void)decompose(local, ParallelContext{}); }
            catch (const std::invalid_argument&) { rejected_local = true; }
            require(rejected_local, "serial decomposition accepted an incomplete mesh");
            const auto owned = detail::ownedCellCount(local);
            require(parallel.sum(owned) == global.cellCount(), "partition lost cells");
            require(owned == global.cellCount() / parallel.size +
                        (parallel.rank < global.cellCount() % parallel.size ? 1 : 0),
                    "partition is not balanced for non-divisible cell counts");
            std::vector<int> ownership(static_cast<std::size_t>(global.cellCount()), 0);
            for (Index cell : detail::meshData(local).owned_cells)
                ++ownership[detail::globalCellId(local, cell)];
            std::vector<int> total(ownership.size());
            parallel.sum(ownership.data(), total.data(), global.cellCount());
            for (int count : total) require(count == 1, "global cell ownership is not unique");
            HaloExchange halo(local, parallel);
            require(halo.plannedBytes(1, true) <= halo.plannedBytes(1),
                    "first-layer communication exceeds full halo");
            std::vector<double> values(static_cast<std::size_t>(local.cellCount()), -1.0);
            for (Index cell : detail::meshData(local).owned_cells)
                values[cell] = detail::globalCellId(local, cell) + 0.5;
            halo.exchange(values);
            for (Index cell = 0; cell < local.cellCount(); ++cell)
                require(values[cell] == detail::globalCellId(local, cell) + 0.5,
                        "halo owner mapping is inconsistent");
            int interfaces = 0;
            for (Index f = 0; f < local.faceCount(); ++f) {
                if (local.neighbour(f) != invalid_index &&
                    detail::isOwned(local, local.owner(f)) !=
                    detail::isOwned(local, local.neighbour(f))) ++interfaces;
            }
            const int cut = parallel.sum(interfaces) / 2;
            const int ghost_cells = parallel.sum(local.cellCount() - owned);
            const int max_local_faces = parallel.maximum(local.faceCount());
            const int total_local_faces = parallel.sum(local.faceCount());
            const int total_bytes = parallel.sum(static_cast<int>(halo.plannedBytes(1)));
            const int first_layer_bytes = parallel.sum(static_cast<int>(halo.plannedBytes(1, true)));
            if (parallel.rank == 0) std::cout << "partition ranks=" << parallel.size
                << " cells=" << global.cellCount() << " cutFaces=" << cut
                << " ghostCellsTotal=" << ghost_cells
                << " localFaceImbalance=" << double(max_local_faces) * parallel.size / total_local_faces
                << " haloBytesTotal=" << total_bytes
                << " firstLayerBytesTotal=" << first_layer_bytes << std::endl;
            // A regular plane should have boundary-sized communication, not a
            // space-filling collection of one-cell-wide growth paths.
            if (shape[0] == 64)
                require(cut <= 64 * parallel.size, "regular-grid partition has excessive interface area");
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        MPI_Abort(parallel.communicator, 1);
        return 1;
    }
    detail::finalizePetscSession();
    MPI_Finalize();
    return 0;
}
