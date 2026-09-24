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
        require(parallel.size == 2 || parallel.size == 4,
                "partition test requires two or four ranks");
        const auto global = makeHexBox({64, 64, 1}, {0, 0, 0}, {1, 1, 0.01});
        const auto local = decompose(global, parallel);
        const auto owned = detail::ownedCellCount(local);
        require(parallel.sum(owned) == global.cellCount(), "partition lost cells");
        require(owned == global.cellCount() / parallel.size, "partition is not balanced");
        int interfaces = 0;
        for (Index f = 0; f < local.faceCount(); ++f) {
            if (local.neighbour(f) != invalid_index &&
                detail::isOwned(local, local.owner(f)) !=
                detail::isOwned(local, local.neighbour(f))) ++interfaces;
        }
        const int cut = parallel.sum(interfaces) / 2;
        if (parallel.rank == 0) std::cout << "partition ranks=" << parallel.size
            << " cutFaces=" << cut << std::endl;
        // A regular plane should have boundary-sized communication, not a
        // space-filling collection of one-cell-wide growth paths.
        require(cut <= 64 * parallel.size, "regular-grid partition has excessive interface area");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        MPI_Abort(parallel.communicator, 1);
        return 1;
    }
    detail::finalizePetscSession();
    MPI_Finalize();
    return 0;
}
