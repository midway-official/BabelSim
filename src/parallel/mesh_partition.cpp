#include "internal/mesh_partition.h"
#include "internal/mesh_access.h"

#ifdef BABELSIM_USE_METIS
#include <metis.h>
#endif

#include <algorithm>
#include <cstdint>
#include <deque>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace babelsim::detail {
namespace {

#ifndef BABELSIM_USE_METIS
std::vector<Index> graphGrowthOwners(const Mesh& mesh, int partitions) {
    const Index cells = mesh.cellCount();
    if (partitions <= 0 || cells < partitions) {
        throw std::invalid_argument("graph partition has an invalid number of parts");
    }
    std::vector<Index> capacities(static_cast<std::size_t>(partitions));
    std::vector<Index> filled(static_cast<std::size_t>(partitions), 0);
    for (int part = 0; part < partitions; ++part) {
        capacities[static_cast<std::size_t>(part)] =
            cells / partitions + (part < cells % partitions ? 1 : 0);
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
            if (seeded[static_cast<std::size_t>(cell)])
                continue;
            if (seed == invalid_index ||
                nearest[static_cast<std::size_t>(cell)] > nearest[static_cast<std::size_t>(seed)] ||
                (nearest[static_cast<std::size_t>(cell)] ==
                     nearest[static_cast<std::size_t>(seed)] &&
                 detail::globalCellId(mesh, cell) < detail::globalCellId(mesh, seed))) {
                seed = cell;
            }
        }
        if (seed == invalid_index)
            throw std::logic_error("graph partition seed selection failed");
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
                    nearest[static_cast<std::size_t>(neighbour)] <= next)
                    continue;
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
            if (filled[static_cast<std::size_t>(part)] >=
                capacities[static_cast<std::size_t>(part)])
                continue;
            bool assigned = false;
            std::deque<Index>& frontier = frontiers[static_cast<std::size_t>(part)];
            while (!frontier.empty() && !assigned) {
                const Index cell = frontier.front();
                frontier.pop_front();
                for (Index neighbour : mesh.cellNeighbours(cell)) {
                    if (neighbour == invalid_index ||
                        owners[static_cast<std::size_t>(neighbour)] != invalid_index)
                        continue;
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
        if (grew)
            continue;

        // Disconnected components have no frontier edge; seed the next component
        // in a non-full part, still without referring to any geometric axes.
        int part = 0;
        while (part < partitions && filled[static_cast<std::size_t>(part)] >=
                                        capacities[static_cast<std::size_t>(part)]) {
            ++part;
        }
        Index seed = invalid_index;
        for (Index cell = 0; cell < cells; ++cell) {
            if (owners[static_cast<std::size_t>(cell)] == invalid_index &&
                (seed == invalid_index ||
                 detail::globalCellId(mesh, cell) < detail::globalCellId(mesh, seed))) {
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

#else

// The policy is selected by the Python mesh-graph study: exact cell quotas,
// least cut-cost boundary transfers, and a bounded proof of safe removal.
// The proof only selects transfers; it never limits stencil/halo coverage.
class QuotaBalance {
public:
    QuotaBalance(const Mesh& mesh, int parts, std::vector<Index>& owners)
        : m_mesh(mesh), m_parts(parts), m_owners(owners), m_counts(parts, 0), m_capacity(parts) {
        for (Index owner : owners) {
            if (owner < 0 || owner >= parts)
                throw std::runtime_error("METIS returned invalid ownership");
            ++m_counts[owner];
        }
        for (int part = 0; part < parts; ++part)
            m_capacity[part] = mesh.cellCount() / parts + (part < mesh.cellCount() % parts ? 1 : 0);
    }

    void apply() {
        if (m_counts == m_capacity)
            return;
        m_recorded.assign(m_owners.size(), 0);
        m_visited.assign(m_owners.size(), 0);
        m_goals.assign(m_owners.size(), 0);
        m_mobility.assign(m_owners.size(), 0);
        m_epoch.assign(m_owners.size(), ~std::uint64_t{0});
        m_versions.assign(m_owners.size(), 0);
        m_targets.assign(m_parts, 0);
        m_queue.reserve(256);
        for (Index cell = 0; cell < m_mesh.cellCount(); ++cell) {
            for (Index neighbour : m_mesh.cellNeighbours(cell)) {
                if (m_owners[cell] != m_owners[neighbour]) {
                    record(cell);
                    break;
                }
            }
        }
        rebuildHeap();
        while (m_counts != m_capacity) {
            Index chosen = invalid_index, target = invalid_index;
            while (!m_heap.empty()) {
                const auto candidate = m_heap.top();
                m_heap.pop();
                const Index cell = std::get<2>(candidate), peer = std::get<3>(candidate);
                if (std::get<4>(candidate) != m_versions[cell] ||
                    m_counts[m_owners[cell]] <= m_capacity[m_owners[cell]] ||
                    m_counts[peer] >= m_capacity[peer] || !movable(cell))
                    continue;
                chosen = cell;
                target = peer;
                break;
            }
            if (chosen != invalid_index) {
                move(chosen, target, true);
                continue;
            }

            // A balanced part can separate a surplus from a deficit. One
            // witness per directed part edge suffices to find a safe route.
            std::vector<std::vector<Index>> adjacency(m_parts), safe(m_parts);
            for (Index cell : m_boundary) {
                const Index source = m_owners[cell];
                for (Index neighbour : m_mesh.cellNeighbours(cell)) {
                    const Index peer = m_owners[neighbour];
                    if (source == peer)
                        continue;
                    adjacency[source].push_back(peer);
                    if (std::find(safe[source].begin(), safe[source].end(), peer) ==
                            safe[source].end() &&
                        movable(cell))
                        safe[source].push_back(peer);
                }
            }
            for (auto* graph : {&adjacency, &safe}) {
                for (auto& peers : *graph) {
                    std::sort(peers.begin(), peers.end());
                    peers.erase(std::unique(peers.begin(), peers.end()), peers.end());
                }
            }
            auto path = balancePath(safe);
            if (path.empty())
                path = balancePath(adjacency);
            if (path.empty()) {
                // Disconnected input can require transfer without an
                // interface. Exact quotas and unique ownership take priority.
                target = 0;
                while (m_counts[target] >= m_capacity[target])
                    ++target;
                move(bestTransfer(invalid_index, target), target, false);
            } else {
                for (std::size_t step = 1; step < path.size(); ++step)
                    move(bestTransfer(path[step - 1], path[step]), path[step], false);
            }
            // Routes can temporarily turn a balanced part into a surplus.
            // Direct transfers update only the moved cell and its neighbours.
            rebuildHeap();
        }
    }

private:
    using Candidate = std::tuple<int, Index, Index, Index, std::uint64_t>;
    const Mesh& m_mesh;
    int m_parts;
    std::vector<Index>& m_owners;
    std::vector<Index> m_counts, m_capacity, m_boundary, m_queue;
    std::vector<unsigned char> m_recorded;
    std::vector<std::uint64_t> m_visited, m_goals;
    std::vector<unsigned char> m_mobility;
    std::vector<std::uint64_t> m_epoch, m_versions, m_targets;
    std::uint64_t m_state = 0, m_sequence = 0, m_target_sequence = 0;
    std::priority_queue<Candidate, std::vector<Candidate>, std::greater<Candidate>> m_heap;

    void record(Index cell) {
        if (!m_recorded[cell]) {
            m_recorded[cell] = 1;
            m_boundary.push_back(cell);
        }
    }

    int cost(Index cell, Index target) const {
        int value = 0;
        for (Index neighbour : m_mesh.cellNeighbours(cell))
            value += (m_owners[neighbour] == m_owners[cell]) - (m_owners[neighbour] == target);
        return value;
    }

    bool movable(Index cell) {
        if (m_epoch[cell] == m_state)
            return m_mobility[cell] != 0;
        const Index source = m_owners[cell];
        const std::uint64_t stamp = ++m_sequence;
        Index seed = invalid_index, remaining = 0;
        for (Index neighbour : m_mesh.cellNeighbours(cell)) {
            if (m_owners[neighbour] != source)
                continue;
            m_goals[neighbour] = stamp;
            seed = neighbour;
            ++remaining;
        }
        bool result = remaining <= 1;
        if (!result) {
            m_visited[cell] = m_visited[seed] = stamp;
            --remaining;
            m_queue.clear();
            m_queue.push_back(seed);
            for (std::size_t head = 0; head < m_queue.size() && head < 64 && !result; ++head) {
                for (Index neighbour : m_mesh.cellNeighbours(m_queue[head])) {
                    if (m_owners[neighbour] != source || m_visited[neighbour] == stamp)
                        continue;
                    m_visited[neighbour] = stamp;
                    if (m_goals[neighbour] == stamp && --remaining == 0) {
                        result = true;
                        break;
                    }
                    m_queue.push_back(neighbour);
                }
            }
        }
        m_epoch[cell] = m_state;
        m_mobility[cell] = static_cast<unsigned char>(result);
        return result;
    }

    void push(Index cell) {
        const Index source = m_owners[cell];
        if (m_counts[source] <= m_capacity[source])
            return;
        const std::uint64_t stamp = ++m_target_sequence;
        for (Index neighbour : m_mesh.cellNeighbours(cell)) {
            const Index target = m_owners[neighbour];
            if (source == target || m_counts[target] >= m_capacity[target] ||
                m_targets[target] == stamp)
                continue;
            m_targets[target] = stamp;
            m_heap.emplace(cost(cell, target), globalCellId(m_mesh, cell), cell, target,
                           m_versions[cell]);
        }
    }

    void rebuildHeap() {
        m_heap = decltype(m_heap){};
        for (Index cell : m_boundary)
            push(cell);
    }

    void move(Index cell, Index target, bool update_heap) {
        if (cell == invalid_index)
            throw std::logic_error("partition quota transfer failed");
        --m_counts[m_owners[cell]];
        ++m_counts[target];
        m_owners[cell] = target;
        ++m_state;
        ++m_versions[cell];
        record(cell);
        if (update_heap)
            push(cell);
        for (Index neighbour : m_mesh.cellNeighbours(cell)) {
            ++m_versions[neighbour];
            record(neighbour);
            if (update_heap)
                push(neighbour);
        }
    }

    Index bestTransfer(Index source, Index target) {
        Index chosen = invalid_index;
        std::pair<int, Index> best{std::numeric_limits<int>::max(),
                                   std::numeric_limits<Index>::max()};
        for (int attempt = 0; attempt < 2; ++attempt) {
            const auto consider = [&](Index cell) {
                if (source == invalid_index) {
                    if (m_counts[m_owners[cell]] <= m_capacity[m_owners[cell]])
                        return;
                } else {
                    if (m_owners[cell] != source)
                        return;
                    bool adjacent = false;
                    for (Index neighbour : m_mesh.cellNeighbours(cell))
                        adjacent |= m_owners[neighbour] == target;
                    if (!adjacent)
                        return;
                }
                const auto score = std::make_pair(cost(cell, target), globalCellId(m_mesh, cell));
                if (score < best && (attempt == 1 || movable(cell))) {
                    best = score;
                    chosen = cell;
                }
            };
            if (source == invalid_index) {
                for (Index cell = 0; cell < m_mesh.cellCount(); ++cell)
                    consider(cell);
            } else {
                for (Index cell : m_boundary)
                    consider(cell);
            }
            if (chosen != invalid_index)
                break;
        }
        return chosen;
    }

    std::vector<Index> balancePath(const std::vector<std::vector<Index>>& adjacency) const {
        std::vector<Index> parent(m_parts, invalid_index), queue;
        for (Index part = 0; part < m_parts; ++part) {
            if (m_counts[part] > m_capacity[part]) {
                parent[part] = part;
                queue.push_back(part);
            }
        }
        Index deficit = invalid_index;
        for (std::size_t head = 0; head < queue.size() && deficit == invalid_index; ++head) {
            for (Index target : adjacency[queue[head]]) {
                if (parent[target] != invalid_index)
                    continue;
                parent[target] = queue[head];
                if (m_counts[target] < m_capacity[target]) {
                    deficit = target;
                    break;
                }
                queue.push_back(target);
            }
        }
        std::vector<Index> path;
        if (deficit != invalid_index) {
            path.push_back(deficit);
            while (parent[path.back()] != path.back())
                path.push_back(parent[path.back()]);
            std::reverse(path.begin(), path.end());
        }
        return path;
    }
};

bool connectedGraph(const std::vector<idx_t>& offsets, const std::vector<idx_t>& neighbours) {
    std::vector<unsigned char> visited(offsets.size() - 1U, 0);
    std::vector<idx_t> queue{0};
    visited[0] = 1;
    for (std::size_t head = 0; head < queue.size(); ++head) {
        for (idx_t edge = offsets[queue[head]]; edge < offsets[queue[head] + 1]; ++edge) {
            const idx_t cell = neighbours[edge];
            if (!visited[cell]) {
                visited[cell] = 1;
                queue.push_back(cell);
            }
        }
    }
    return queue.size() == visited.size();
}

// Partition choices are sensitive to the last bit of these proportions.
// Keep the once-per-part normalization outside the solver's reciprocal-math
// optimization: the Python policy study uses a rounded quotient, not a rounded
// reciprocal followed by multiplication. Double intermediates also preserve
// integer cell counts beyond the single-precision exact-integer range.
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("no-fast-math"), noinline))
#endif
real_t
normalizedTarget(Index quota, Index cells) {
    return static_cast<real_t>(static_cast<double>(quota) / static_cast<double>(cells));
}
#endif

} // namespace

std::vector<Index> partitionOwners(const Mesh& mesh, int partitions) {
    if (partitions <= 0 || mesh.cellCount() < partitions)
        throw std::invalid_argument("graph partition has an invalid number of parts");
    if (partitions == 1)
        return std::vector<Index>(static_cast<std::size_t>(mesh.cellCount()), 0);
#ifndef BABELSIM_USE_METIS
    return graphGrowthOwners(mesh, partitions);
#else
    idx_t cells = mesh.cellCount(), constraints = 2, parts = partitions, objective = 0;
    std::vector<idx_t> offsets(static_cast<std::size_t>(cells) + 1U, 0), neighbours;
    std::vector<idx_t> weights(static_cast<std::size_t>(cells) * 2U),
        owners(static_cast<std::size_t>(cells));
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        for (Index neighbour : mesh.cellNeighbours(cell))
            neighbours.push_back(neighbour);
        if (neighbours.size() > static_cast<std::size_t>(std::numeric_limits<idx_t>::max()))
            throw std::overflow_error("mesh adjacency exceeds METIS index width");
        offsets[static_cast<std::size_t>(cell) + 1U] = static_cast<idx_t>(neighbours.size());
        weights[2U * static_cast<std::size_t>(cell)] = 1;
        weights[2U * static_cast<std::size_t>(cell) + 1U] =
            static_cast<idx_t>(mesh.cellFaces(cell).size());
    }
    if (neighbours.empty()) {
        // Independent cells have neither cut faces nor halo communication.
        std::vector<Index> result(static_cast<std::size_t>(cells));
        Index cell = 0;
        for (int part = 0; part < partitions; ++part) {
            const Index quota =
                mesh.cellCount() / partitions + (part < mesh.cellCount() % partitions ? 1 : 0);
            for (Index count = 0; count < quota; ++count)
                result[cell++] = part;
        }
        return result;
    }
    std::vector<real_t> targets(static_cast<std::size_t>(parts) * 2U);
    for (int part = 0; part < parts; ++part) {
        const Index quota = mesh.cellCount() / parts + (part < mesh.cellCount() % parts ? 1 : 0);
        targets[2U * static_cast<std::size_t>(part)] =
            targets[2U * static_cast<std::size_t>(part) + 1U] =
                normalizedTarget(quota, mesh.cellCount());
    }
    // Cell counts govern storage and algebra; face incidence approximates the
    // variable stencil/assembly work of general polyhedral cells. The exact
    // cell quotas are restored after METIS's tolerance-based optimization.
    real_t imbalance[2]{1.001, 1.05};
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_SEED] = 0;
    options[METIS_OPTION_NCUTS] = 2;
    options[METIS_OPTION_MINCONN] = 1;
    options[METIS_OPTION_CONTIG] = connectedGraph(offsets, neighbours) ? 1 : 0;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
    options[METIS_OPTION_NUMBERING] = 0;
    const int status = METIS_PartGraphKway(&cells, &constraints, offsets.data(), neighbours.data(),
                                           weights.data(), nullptr, nullptr, &parts, targets.data(),
                                           imbalance, options, &objective, owners.data());
    if (status != METIS_OK)
        throw std::runtime_error("METIS graph partition failed: " + std::to_string(status));
    std::vector<Index> result(owners.begin(), owners.end());
    // METIS scratch is no longer needed while restoring exact quotas.
    std::vector<idx_t>{}.swap(offsets);
    std::vector<idx_t>{}.swap(neighbours);
    std::vector<idx_t>{}.swap(weights);
    std::vector<idx_t>{}.swap(owners);
    std::vector<real_t>{}.swap(targets);
    QuotaBalance(mesh, partitions, result).apply();
    return result;
#endif
}

} // namespace babelsim::detail
