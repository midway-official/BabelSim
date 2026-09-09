#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include "babelsim/parallel.h"

#include "test_util.h"

#include <mpi.h>

#include <iostream>
#include <deque>
#include <utility>
#include <vector>

using namespace babelsim;

namespace {

Mesh checkerboardHexes() {
    const Mesh base = makeHexBox({2, 2, 2}, {0, 0, 0}, {1, 1, 1});
    const std::array<Index, 8> permutation{{0, 3, 5, 6, 1, 2, 4, 7}};
    std::vector<Vec3> vertices;
    std::vector<std::array<Index, 8>> cells;
    std::vector<PatchSpec> patches;
    std::vector<BoundaryFaceSpec> boundaries;
    for (Index vertex = 0; vertex < base.vertexCount(); ++vertex) vertices.push_back(base.vertex(vertex));
    for (Index source : permutation) cells.push_back(base.cellVertices(source));
    for (Index patch = 0; patch < base.patchCount(); ++patch) {
        patches.push_back({base.patchName(patch), base.patchKind(patch)});
    }
    for (Index face = 0; face < base.faceCount(); ++face) {
        if (base.boundaryFace(face)) {
            boundaries.push_back({base.faceVertices(face), detail::meshData(base).face_patch[face]});
        }
    }
    return Mesh::unstructured(std::move(vertices), std::move(cells), std::move(patches),
                              std::move(boundaries));
}

}  // namespace

int main(int argc, char* argv[]) {
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) return 1;
    const ParallelContext parallel = ParallelContext::world();
    try {
        require(parallel.size == 2, "parallel_unstructured_test requires two MPI ranks");
        const Mesh local = decompose(checkerboardHexes(), parallel, 1);
        require(detail::ownedCellCount(local) == 4 && local.cellCount() > 4,
                "topology partition did not create a balanced owned-plus-ghost mesh");
        std::vector<bool> reached(static_cast<std::size_t>(local.cellCount()), false);
        std::deque<Index> connected{detail::meshData(local).owned_cells.front()};
        reached[static_cast<std::size_t>(connected.front())] = true;
        while (!connected.empty()) {
            const Index cell = connected.front();
            connected.pop_front();
            for (Index neighbour : detail::meshData(local).cell_neighbours[static_cast<std::size_t>(cell)]) {
                if (neighbour != invalid_index && detail::isOwned(local, neighbour) &&
                    !reached[static_cast<std::size_t>(neighbour)]) {
                    reached[static_cast<std::size_t>(neighbour)] = true;
                    connected.push_back(neighbour);
                }
            }
        }
        Index connected_owned = 0;
        for (Index cell : detail::meshData(local).owned_cells) {
            if (reached[static_cast<std::size_t>(cell)]) ++connected_owned;
        }
        require(connected_owned == detail::ownedCellCount(local),
                "graph partition followed the scrambled input order instead of cell adjacency");
        Index remote_links = 0;
        for (Index cell : detail::meshData(local).owned_cells) {
            for (Index neighbour : detail::meshData(local).cell_neighbours[static_cast<std::size_t>(cell)]) {
                if (neighbour != invalid_index && !detail::isOwned(local, neighbour)) ++remote_links;
            }
        }
        require(remote_links > 0,
                "topology partition did not retain a remote graph interface");
        ScalarField values(local, FieldLocation::Cell, "globalId", -1.0);
        for (Index cell : detail::meshData(local).owned_cells) {
            detail::fieldData(values)[cell] = detail::globalCellId(local, cell) + 0.5;
        }
        HaloExchange halo(local, parallel);
        halo.exchange(values);
        for (Index cell = 0; cell < local.cellCount(); ++cell) {
            require(near(detail::fieldData(values)[cell], detail::globalCellId(local, cell) + 0.5),
                    "graph halo exchange did not resolve an arbitrary-topology neighbour");
        }
        if (parallel.rank == 0) {
            std::cout << "parallel_unstructured_test: scrambled hex graph partition and halo passed\n";
        }
    } catch (const std::exception& error) {
        std::cerr << "parallel_unstructured_test: " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
