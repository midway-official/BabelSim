#include "babelsim/parallel.h"
#include "internal/mesh_access.h"
#include "internal/petsc_session.h"
#include "test_util.h"
#include <mpi.h>
#include <iostream>

using namespace babelsim;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    const auto parallel = ParallelContext::world();
    try {
        require(parallel.size >= 2 && parallel.size <= 8,
                "partition test requires two to eight ranks");
        for (const auto shape : {std::array<Index, 3>{64, 64, 1},
                                 std::array<Index, 3>{17, 11, 3},
                                 std::array<Index, 3>{128, 2, 1}}) {
            const auto global = makeHexBox({shape[0], shape[1], shape[2]}, {0, 0, 0}, {1, 1, 0.01});
            const auto local = decompose(global, parallel);
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
