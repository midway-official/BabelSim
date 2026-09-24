#include "internal/petsc_index_map.h"

#include "internal/mesh_access.h"
#include "babelsim/mpi_support.h"

#include <algorithm>
#include <climits>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

namespace babelsim::detail {
namespace {

PetscInt checkedPetscIndex(std::uint64_t value, const char* what) {
    if (value > static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max()))
        throw std::overflow_error(std::string(what) + " exceeds PETSc PetscInt range");
    return static_cast<PetscInt>(value);
}

}  // namespace

PetscIndexMap::PetscIndexMap(const Mesh& mesh, const ParallelContext& parallel)
    : m_mesh(&mesh), m_local_size(checkedPetscIndex(
          static_cast<std::uint64_t>(ownedCellCount(mesh)), "local owned-cell count")),
      m_local_to_global(static_cast<std::size_t>(mesh.cellCount()), PetscInt{-1}) {
    mesh.validate();
    parallel.validate();
    std::vector<PetscInt> counts(static_cast<std::size_t>(parallel.size));
    if (parallel.size == 1) {
        counts[0] = m_local_size;
    } else {
        int local = mpiCount(static_cast<std::size_t>(m_local_size), "owned row count");
        std::vector<int> gathered(static_cast<std::size_t>(parallel.size));
        checkMpi(MPI_Allgather(&local, 1, MPI_INT, gathered.data(), 1, MPI_INT,
                               parallel.communicator), "MPI_Allgather(PETSc row sizes)");
        for (int rank = 0; rank < parallel.size; ++rank) {
            if (gathered[static_cast<std::size_t>(rank)] < 0)
                throw std::runtime_error("negative PETSc local row count");
            counts[static_cast<std::size_t>(rank)] = checkedPetscIndex(
                static_cast<std::uint64_t>(gathered[static_cast<std::size_t>(rank)]),
                "rank row count");
        }
    }

    std::vector<PetscInt> offsets(static_cast<std::size_t>(parallel.size) + 1U, 0);
    for (int rank = 0; rank < parallel.size; ++rank) {
        const std::uint64_t next = static_cast<std::uint64_t>(offsets[static_cast<std::size_t>(rank)]) +
            static_cast<std::uint64_t>(counts[static_cast<std::size_t>(rank)]);
        offsets[static_cast<std::size_t>(rank) + 1U] = checkedPetscIndex(next, "global row count");
    }
    m_row_begin = offsets[static_cast<std::size_t>(parallel.rank)];
    m_row_end = offsets[static_cast<std::size_t>(parallel.rank) + 1U];
    m_global_size = offsets.back();
    if (m_row_end - m_row_begin != m_local_size)
        throw std::logic_error("PETSc row ownership prefix is inconsistent");

    for (Index cell : meshData(mesh).owned_cells) {
        const PetscInt local = checkedPetscIndex(
            static_cast<std::uint64_t>(ownedIndex(mesh, cell)), "owned local row");
        m_local_to_global.at(static_cast<std::size_t>(cell)) = m_row_begin + local;
    }
    if (parallel.size == 1) return;

    std::vector<std::vector<Index>> ids_by_owner(static_cast<std::size_t>(parallel.size));
    std::vector<std::vector<Index>> cells_by_owner(static_cast<std::size_t>(parallel.size));
    for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
        if (isOwned(mesh, cell)) continue;
        const Index owner = cellOwnerRank(mesh, cell);
        if (owner < 0 || owner >= parallel.size)
            throw std::invalid_argument("ghost cell has invalid owner rank for PETSc mapping");
        ids_by_owner[static_cast<std::size_t>(owner)].push_back(globalCellId(mesh, cell));
        cells_by_owner[static_cast<std::size_t>(owner)].push_back(cell);
    }

    std::vector<int> request_counts(static_cast<std::size_t>(parallel.size), 0);
    std::vector<int> request_offsets(static_cast<std::size_t>(parallel.size), 0);
    std::vector<Index> requests;
    std::vector<Index> request_cells;
    for (int rank = 0; rank < parallel.size; ++rank) {
        auto& peer_ids = ids_by_owner[static_cast<std::size_t>(rank)];
        auto& peer_cells = cells_by_owner[static_cast<std::size_t>(rank)];
        if (peer_ids.size() != peer_cells.size())
            throw std::logic_error("PETSc ghost request buffers disagree");
        request_offsets[static_cast<std::size_t>(rank)] = mpiCount(requests.size(), "request offset");
        request_counts[static_cast<std::size_t>(rank)] = mpiCount(peer_ids.size(), "request count");
        requests.insert(requests.end(), peer_ids.begin(), peer_ids.end());
        request_cells.insert(request_cells.end(), peer_cells.begin(), peer_cells.end());
    }

    std::vector<int> incoming_counts(static_cast<std::size_t>(parallel.size), 0);
    checkMpi(MPI_Alltoall(request_counts.data(), 1, MPI_INT, incoming_counts.data(), 1,
                          MPI_INT, parallel.communicator), "MPI_Alltoall(PETSc map counts)");
    std::vector<int> incoming_offsets(static_cast<std::size_t>(parallel.size), 0);
    std::size_t incoming_total = 0;
    for (int rank = 0; rank < parallel.size; ++rank) {
        incoming_offsets[static_cast<std::size_t>(rank)] = mpiCount(incoming_total, "incoming request offset");
        if (incoming_counts[static_cast<std::size_t>(rank)] < 0)
            throw std::runtime_error("negative PETSc mapping request count");
        incoming_total += static_cast<std::size_t>(incoming_counts[static_cast<std::size_t>(rank)]);
    }
    mpiCount(incoming_total, "incoming request total");
    std::vector<Index> incoming_ids(incoming_total);
    Index dummy = 0;
    checkMpi(MPI_Alltoallv(
        requests.empty() ? &dummy : requests.data(), request_counts.data(), request_offsets.data(), MPI_INT,
        incoming_ids.empty() ? &dummy : incoming_ids.data(), incoming_counts.data(), incoming_offsets.data(), MPI_INT,
        parallel.communicator), "MPI_Alltoallv(PETSc ghost IDs)");

    std::map<Index, PetscInt> owned_ids;
    for (Index cell : meshData(mesh).owned_cells) {
        if (!owned_ids.emplace(globalCellId(mesh, cell), m_local_to_global[static_cast<std::size_t>(cell)]).second)
            throw std::logic_error("duplicate owned global cell ID in PETSc mapping");
    }
    std::vector<PetscInt> incoming_rows(incoming_total);
    for (std::size_t i = 0; i < incoming_ids.size(); ++i) {
        const auto found = owned_ids.find(incoming_ids[i]);
        if (found == owned_ids.end())
            throw std::runtime_error("PETSc row request reached a rank that does not own the cell");
        incoming_rows[i] = found->second;
    }
    std::vector<PetscInt> requested_rows(requests.size());
    checkMpi(MPI_Alltoallv(
        incoming_rows.empty() ? nullptr : incoming_rows.data(), incoming_counts.data(), incoming_offsets.data(), MPIU_INT,
        requested_rows.empty() ? nullptr : requested_rows.data(), request_counts.data(), request_offsets.data(), MPIU_INT,
        parallel.communicator), "MPI_Alltoallv(PETSc row IDs)");
    for (std::size_t i = 0; i < request_cells.size(); ++i)
        m_local_to_global.at(static_cast<std::size_t>(request_cells[i])) = requested_rows[i];
    if (std::find(m_local_to_global.begin(), m_local_to_global.end(), PetscInt{-1}) != m_local_to_global.end())
        throw std::logic_error("PETSc cell map left a local cell unresolved");
}

PetscInt PetscIndexMap::globalCell(Index local_cell) const {
    if (local_cell < 0 || static_cast<std::size_t>(local_cell) >= m_local_to_global.size())
        throw std::out_of_range("local cell index is outside PETSc map");
    return m_local_to_global[static_cast<std::size_t>(local_cell)];
}

PetscInt PetscIndexMap::row(Index owned_cell) const {
    if (!isOwned(*m_mesh, owned_cell))
        throw std::invalid_argument("PETSc row requested for a ghost cell");
    return globalCell(owned_cell);
}

}  // namespace babelsim::detail
