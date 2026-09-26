#include "test_util.h"
#include "babelsim/mesh.h"
#include "internal/mesh_access.h"
#include <algorithm>
#include <chrono>
#include <deque>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace babelsim;
std::vector<Index> before(const Mesh& mesh, int partitions) {
    const Index cells = mesh.cellCount();
    if (partitions <= 0 || cells < partitions) {
        throw std::invalid_argument("graph partition has an invalid number of parts");
    }
    std::vector<Index> capacities(static_cast<std::size_t>(partitions));
    std::vector<Index> filled(static_cast<std::size_t>(partitions), 0);
    for (int part = 0; part < partitions; ++part) {
        capacities[static_cast<std::size_t>(part)] = cells / partitions +
            (part < cells % partitions ? 1 : 0);
    }

    // Select dispersed deterministic seeds using only graph distance.  The cell ID
    // decides ties, but never defines the partition boundary.
    const Index unreachable = std::numeric_limits<Index>::max();
    std::vector<Index> nearest(static_cast<std::size_t>(cells), unreachable);
    std::vector<Index> seeds;
    seeds.reserve(static_cast<std::size_t>(partitions));
    for (int part = 0; part < partitions; ++part) {
        Index seed = invalid_index;
        for (Index cell = 0; cell < cells; ++cell) {
            if (std::find(seeds.begin(), seeds.end(), cell) != seeds.end()) continue;
            if (seed == invalid_index ||
                nearest[static_cast<std::size_t>(cell)] > nearest[static_cast<std::size_t>(seed)] ||
                (nearest[static_cast<std::size_t>(cell)] == nearest[static_cast<std::size_t>(seed)] &&
                 detail::globalCellId(mesh, cell) < detail::globalCellId(mesh, seed))) {
                seed = cell;
            }
        }
        if (seed == invalid_index) throw std::logic_error("graph partition seed selection failed");
        seeds.push_back(seed);

        std::vector<Index> distance(static_cast<std::size_t>(cells), invalid_index);
        std::deque<Index> frontier{seed};
        distance[static_cast<std::size_t>(seed)] = 0;
        while (!frontier.empty()) {
            const Index cell = frontier.front();
            frontier.pop_front();
            for (Index neighbour : mesh.cellNeighbours(cell)) {
                if (neighbour == invalid_index ||
                    distance[static_cast<std::size_t>(neighbour)] != invalid_index) continue;
                distance[static_cast<std::size_t>(neighbour)] =
                    distance[static_cast<std::size_t>(cell)] + 1;
                frontier.push_back(neighbour);
            }
        }
        for (Index cell = 0; cell < cells; ++cell) {
            const Index path = distance[static_cast<std::size_t>(cell)];
            if (path != invalid_index) {
                nearest[static_cast<std::size_t>(cell)] = std::min(
                    nearest[static_cast<std::size_t>(cell)], path);
            }
        }
    }

    std::vector<Index> owners(static_cast<std::size_t>(cells), invalid_index);
    std::vector<std::deque<Index>> frontiers(static_cast<std::size_t>(partitions));
    for (int part = 0; part < partitions; ++part) {
        const Index seed = seeds[static_cast<std::size_t>(part)];
        owners[static_cast<std::size_t>(seed)] = part;
        ++filled[static_cast<std::size_t>(part)];
        frontiers[static_cast<std::size_t>(part)].push_back(seed);
    }

    Index remaining = cells - partitions;
    while (remaining > 0) {
        bool grew = false;
        for (int part = 0; part < partitions; ++part) {
            if (filled[static_cast<std::size_t>(part)] >= capacities[static_cast<std::size_t>(part)]) continue;
            bool assigned = false;
            std::deque<Index>& frontier = frontiers[static_cast<std::size_t>(part)];
            while (!frontier.empty() && !assigned) {
                const Index cell = frontier.front();
                frontier.pop_front();
                for (Index neighbour : mesh.cellNeighbours(cell)) {
                    if (neighbour == invalid_index ||
                        owners[static_cast<std::size_t>(neighbour)] != invalid_index) continue;
                    owners[static_cast<std::size_t>(neighbour)] = part;
                    ++filled[static_cast<std::size_t>(part)];
                    --remaining;
                    // The round-robin quota assigns one neighbour at a time.
                    // Keep this parent at the head until all its unassigned
                    // neighbours have been visited; otherwise growth discards
                    // branches and degenerates into thin paths and reseeding.
                    frontier.push_front(cell);
                    frontier.push_back(neighbour);
                    assigned = true;
                    grew = true;
                    break;
                }
            }
        }
        if (grew) continue;

        // Disconnected components have no frontier edge; seed the next component
        // in a non-full part, still without referring to any geometric axes.
        int part = 0;
        while (part < partitions &&
               filled[static_cast<std::size_t>(part)] >= capacities[static_cast<std::size_t>(part)]) {
            ++part;
        }
        Index seed = invalid_index;
        for (Index cell = 0; cell < cells; ++cell) {
            if (owners[static_cast<std::size_t>(cell)] == invalid_index &&
                (seed == invalid_index || detail::globalCellId(mesh, cell) < detail::globalCellId(mesh, seed))) {
                seed = cell;
            }
        }
        if (part == partitions || seed == invalid_index) {
            throw std::logic_error("graph partition growth failed");
        }
        owners[static_cast<std::size_t>(seed)] = part;
        ++filled[static_cast<std::size_t>(part)];
        --remaining;
        frontiers[static_cast<std::size_t>(part)].push_back(seed);
    }
    return owners;
}
std::vector<Index> after(const Mesh& mesh, int partitions) {
    const Index cells = mesh.cellCount();
    if (partitions <= 0 || cells < partitions) {
        throw std::invalid_argument("graph partition has an invalid number of parts");
    }
    std::vector<Index> capacities(static_cast<std::size_t>(partitions));
    std::vector<Index> filled(static_cast<std::size_t>(partitions), 0);
    for (int part = 0; part < partitions; ++part) {
        capacities[static_cast<std::size_t>(part)] = cells / partitions +
            (part < cells % partitions ? 1 : 0);
    }

    // Select dispersed deterministic seeds using only graph distance.  The cell ID
    // decides ties, but never defines the partition boundary.
    const Index unreachable = std::numeric_limits<Index>::max();
    std::vector<Index> nearest(static_cast<std::size_t>(cells), unreachable);
    std::vector<unsigned char> seeded(static_cast<std::size_t>(cells), 0);
    std::vector<Index> seeds;
    seeds.reserve(static_cast<std::size_t>(partitions));
    for (int part = 0; part < partitions; ++part) {
        Index seed = invalid_index;
        for (Index cell = 0; cell < cells; ++cell) {
            if (seeded[static_cast<std::size_t>(cell)]) continue;
            if (seed == invalid_index ||
                nearest[static_cast<std::size_t>(cell)] > nearest[static_cast<std::size_t>(seed)] ||
                (nearest[static_cast<std::size_t>(cell)] == nearest[static_cast<std::size_t>(seed)] &&
                 detail::globalCellId(mesh, cell) < detail::globalCellId(mesh, seed))) {
                seed = cell;
            }
        }
        if (seed == invalid_index) throw std::logic_error("graph partition seed selection failed");
        seeds.push_back(seed);
        seeded[static_cast<std::size_t>(seed)] = 1;

        // Relax only cells whose distance to the nearest seed improves.
        // An existing shorter path also bounds all downstream paths, so a
        // full graph traversal and a fresh N-cell distance array are needless.
        std::deque<Index> frontier{seed};
        nearest[static_cast<std::size_t>(seed)] = 0;
        while (!frontier.empty()) {
            const Index cell = frontier.front();
            frontier.pop_front();
            const Index next = nearest[static_cast<std::size_t>(cell)] + 1;
            for (Index neighbour : mesh.cellNeighbours(cell)) {
                if (neighbour == invalid_index ||
                    nearest[static_cast<std::size_t>(neighbour)] <= next) continue;
                nearest[static_cast<std::size_t>(neighbour)] = next;
                frontier.push_back(neighbour);
            }
        }
    }

    std::vector<Index> owners(static_cast<std::size_t>(cells), invalid_index);
    std::vector<std::deque<Index>> frontiers(static_cast<std::size_t>(partitions));
    for (int part = 0; part < partitions; ++part) {
        const Index seed = seeds[static_cast<std::size_t>(part)];
        owners[static_cast<std::size_t>(seed)] = part;
        ++filled[static_cast<std::size_t>(part)];
        frontiers[static_cast<std::size_t>(part)].push_back(seed);
    }

    Index remaining = cells - partitions;
    while (remaining > 0) {
        bool grew = false;
        for (int part = 0; part < partitions; ++part) {
            if (filled[static_cast<std::size_t>(part)] >= capacities[static_cast<std::size_t>(part)]) continue;
            bool assigned = false;
            std::deque<Index>& frontier = frontiers[static_cast<std::size_t>(part)];
            while (!frontier.empty() && !assigned) {
                const Index cell = frontier.front();
                frontier.pop_front();
                for (Index neighbour : mesh.cellNeighbours(cell)) {
                    if (neighbour == invalid_index ||
                        owners[static_cast<std::size_t>(neighbour)] != invalid_index) continue;
                    owners[static_cast<std::size_t>(neighbour)] = part;
                    ++filled[static_cast<std::size_t>(part)];
                    --remaining;
                    // The round-robin quota assigns one neighbour at a time.
                    // Keep this parent at the head until all its unassigned
                    // neighbours have been visited; otherwise growth discards
                    // branches and degenerates into thin paths and reseeding.
                    frontier.push_front(cell);
                    frontier.push_back(neighbour);
                    assigned = true;
                    grew = true;
                    break;
                }
            }
        }
        if (grew) continue;

        // Disconnected components have no frontier edge; seed the next component
        // in a non-full part, still without referring to any geometric axes.
        int part = 0;
        while (part < partitions &&
               filled[static_cast<std::size_t>(part)] >= capacities[static_cast<std::size_t>(part)]) {
            ++part;
        }
        Index seed = invalid_index;
        for (Index cell = 0; cell < cells; ++cell) {
            if (owners[static_cast<std::size_t>(cell)] == invalid_index &&
                (seed == invalid_index || detail::globalCellId(mesh, cell) < detail::globalCellId(mesh, seed))) {
                seed = cell;
            }
        }
        if (part == partitions || seed == invalid_index) {
            throw std::logic_error("graph partition growth failed");
        }
        owners[static_cast<std::size_t>(seed)] = part;
        ++filled[static_cast<std::size_t>(part)];
        --remaining;
        frontiers[static_cast<std::size_t>(part)].push_back(seed);
    }
    return owners;
}

int main() {
for (int n : {64,128,256}) {
auto mesh=makeHexBox({n,n,1},{0,0,0},{1,1,.01});
for(int p : {2,3,4,8,16}) {
auto a=std::chrono::steady_clock::now();auto x=before(mesh,p);
auto b=std::chrono::steady_clock::now();auto y=after(mesh,p);
auto c=std::chrono::steady_clock::now();
if(x!=y) return 1;
std::cout<<n<<","<<p<<","<<std::chrono::duration<double>(b-a).count()<<","<<std::chrono::duration<double>(c-b).count()<<",identical\n";
}}
}
