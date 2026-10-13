#include "babelsim/parallel.h"
#include "babelsim/mpi_support.h"
#include "internal/mesh_access.h"
#include "internal/field_access.h"

#include <algorithm>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace babelsim {

// Plans translate global requests to local indices once during setup. begin/finish
// reuse those plans and buffers; no mesh lookup or allocation is needed at steady size.
HaloExchange::HaloExchange(const Mesh& mesh, ParallelContext parallel)
    : m_mesh(&mesh), m_parallel(parallel) {
    m_parallel.validate();
    mesh.validate();
    if (!m_parallel.distributed())
        return;
    const auto build = [&](ExchangePlan& plan, bool faces, bool first_layer) {
        std::vector<std::vector<Index>> requested_ids(static_cast<std::size_t>(m_parallel.size));
        std::vector<std::vector<Index>> requested_indices(
            static_cast<std::size_t>(m_parallel.size));
        const Index count = faces ? mesh.faceCount() : mesh.cellCount();
        for (Index entity = 0; entity < count; ++entity) {
            const Index owner =
                faces ? detail::faceOwnerRank(mesh, entity) : detail::cellOwnerRank(mesh, entity);
            if (owner == m_parallel.rank)
                continue;
            if (!faces && first_layer &&
                detail::meshData(mesh).cell_ghost_depths[static_cast<std::size_t>(entity)] != 1)
                continue;
            if (!faces && !first_layer && detail::isOwned(mesh, entity))
                continue;
            if (owner < 0 || owner >= m_parallel.size) {
                throw std::invalid_argument("halo entity has an invalid owner rank");
            }
            requested_ids[static_cast<std::size_t>(owner)].push_back(
                faces ? detail::globalFaceId(mesh, entity) : detail::globalCellId(mesh, entity));
            requested_indices[static_cast<std::size_t>(owner)].push_back(entity);
        }
        plan.receive_counts.assign(static_cast<std::size_t>(m_parallel.size), 0);
        plan.receive_offsets.assign(static_cast<std::size_t>(m_parallel.size), 0);
        std::vector<Index> outgoing_ids;
        for (int peer = 0; peer < m_parallel.size; ++peer) {
            plan.receive_offsets[static_cast<std::size_t>(peer)] =
                detail::mpiCount(outgoing_ids.size(), "halo request offset");
            plan.receive_counts[static_cast<std::size_t>(peer)] = detail::mpiCount(
                requested_ids[static_cast<std::size_t>(peer)].size(), "halo request count");
            outgoing_ids.insert(outgoing_ids.end(),
                                requested_ids[static_cast<std::size_t>(peer)].begin(),
                                requested_ids[static_cast<std::size_t>(peer)].end());
            plan.receive_indices.insert(plan.receive_indices.end(),
                                        requested_indices[static_cast<std::size_t>(peer)].begin(),
                                        requested_indices[static_cast<std::size_t>(peer)].end());
        }
        plan.send_counts.assign(static_cast<std::size_t>(m_parallel.size), 0);
        detail::checkMpi(MPI_Alltoall(plan.receive_counts.data(), 1, MPI_INT,
                                      plan.send_counts.data(), 1, MPI_INT, m_parallel.communicator),
                         "MPI_Alltoall(halo request counts)");
        plan.send_offsets.assign(static_cast<std::size_t>(m_parallel.size), 0);
        std::size_t incoming_size = 0;
        for (int peer = 0; peer < m_parallel.size; ++peer) {
            plan.send_offsets[static_cast<std::size_t>(peer)] =
                detail::mpiCount(incoming_size, "halo request offset");
            incoming_size +=
                static_cast<std::size_t>(plan.send_counts[static_cast<std::size_t>(peer)]);
        }
        std::vector<Index> incoming_ids(incoming_size);
        Index dummy = 0;
        detail::checkMpi(MPI_Alltoallv(outgoing_ids.empty() ? &dummy : outgoing_ids.data(),
                                       plan.receive_counts.data(), plan.receive_offsets.data(),
                                       MPI_INT, incoming_ids.empty() ? &dummy : incoming_ids.data(),
                                       plan.send_counts.data(), plan.send_offsets.data(), MPI_INT,
                                       m_parallel.communicator),
                         "MPI_Alltoallv(halo requests)");
        std::map<Index, Index> owned;
        for (Index entity = 0; entity < count; ++entity) {
            const bool owner = faces ? detail::faceOwnerRank(mesh, entity) == m_parallel.rank
                                     : detail::isOwned(mesh, entity);
            if (owner)
                owned.emplace(faces ? detail::globalFaceId(mesh, entity)
                                    : detail::globalCellId(mesh, entity),
                              entity);
        }
        for (Index id : incoming_ids) {
            const auto local = owned.find(id);
            if (local == owned.end()) {
                throw std::runtime_error(std::string("halo request is not owned by rank ") +
                                         std::to_string(m_parallel.rank) + " (entity=" +
                                         std::to_string(id) + (faces ? ",faces)" : ",cells)"));
            }
            plan.send_indices.push_back(local->second);
        }
    };
    build(m_cells, false, false);
    build(m_first_layer_cells, false, true);
    build(m_faces, true, false);
    buildNeighbourCommunicator(m_cells);
    buildNeighbourCommunicator(m_first_layer_cells);
    buildNeighbourCommunicator(m_faces);
}

void HaloExchange::buildNeighbourCommunicator(ExchangePlan& plan) {
    if (!m_parallel.distributed())
        return;
    plan.outgoing_peers.clear();
    plan.incoming_peers.clear();
    for (int peer = 0; peer < m_parallel.size; ++peer) {
        const std::size_t p = static_cast<std::size_t>(peer);
        // send_counts describe requests received from a peer, therefore the
        // values flow out to that peer; receive_counts are the IDs this rank
        // requested and form the incoming-neighbour list.
        if (plan.send_counts[p] > 0)
            plan.outgoing_peers.push_back(peer);
        if (plan.receive_counts[p] > 0)
            plan.incoming_peers.push_back(peer);
    }
    detail::checkMpi(MPI_Dist_graph_create_adjacent(
                         m_parallel.communicator, static_cast<int>(plan.incoming_peers.size()),
                         plan.incoming_peers.empty() ? nullptr : plan.incoming_peers.data(),
                         MPI_UNWEIGHTED, static_cast<int>(plan.outgoing_peers.size()),
                         plan.outgoing_peers.empty() ? nullptr : plan.outgoing_peers.data(),
                         MPI_UNWEIGHTED, MPI_INFO_NULL, 0, &plan.graph_communicator),
                     "MPI_Dist_graph_create_adjacent(halo)");
}

HaloExchange::~HaloExchange() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
        return;
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized)
        return;
    if (m_cells.graph_communicator != MPI_COMM_NULL)
        MPI_Comm_free(&m_cells.graph_communicator);
    if (m_first_layer_cells.graph_communicator != MPI_COMM_NULL)
        MPI_Comm_free(&m_first_layer_cells.graph_communicator);
    if (m_faces.graph_communicator != MPI_COMM_NULL)
        MPI_Comm_free(&m_faces.graph_communicator);
}

HaloExchange::ExchangePlan::ValueLayout& HaloExchange::valueLayout(ExchangePlan& plan,
                                                                   std::size_t components) {
    for (auto& layout : plan.common_layouts) {
        if (layout.components == components)
            return layout;
    }
    for (auto& layout : plan.common_layouts) {
        if (layout.components != 0)
            continue;
        layout.components = components;
        layout.send_counts.resize(static_cast<std::size_t>(m_parallel.size));
        layout.send_offsets.resize(static_cast<std::size_t>(m_parallel.size));
        layout.receive_counts.resize(static_cast<std::size_t>(m_parallel.size));
        layout.receive_offsets.resize(static_cast<std::size_t>(m_parallel.size));
        for (int peer = 0; peer < m_parallel.size; ++peer) {
            const std::size_t p = static_cast<std::size_t>(peer);
            layout.send_counts[p] = detail::mpiCount(
                static_cast<std::size_t>(plan.send_counts[p]) * components, "halo send values");
            layout.send_offsets[p] = detail::mpiCount(
                static_cast<std::size_t>(plan.send_offsets[p]) * components, "halo send offset");
            layout.receive_counts[p] =
                detail::mpiCount(static_cast<std::size_t>(plan.receive_counts[p]) * components,
                                 "halo receive values");
            layout.receive_offsets[p] =
                detail::mpiCount(static_cast<std::size_t>(plan.receive_offsets[p]) * components,
                                 "halo receive offset");
        }
        layout.neighbour_send_counts.clear();
        layout.neighbour_send_offsets.clear();
        layout.neighbour_receive_counts.clear();
        layout.neighbour_receive_offsets.clear();
        for (const int peer : plan.outgoing_peers) {
            const std::size_t p = static_cast<std::size_t>(peer);
            layout.neighbour_send_counts.push_back(layout.send_counts[p]);
            layout.neighbour_send_offsets.push_back(layout.send_offsets[p]);
        }
        for (const int peer : plan.incoming_peers) {
            const std::size_t p = static_cast<std::size_t>(peer);
            layout.neighbour_receive_counts.push_back(layout.receive_counts[p]);
            layout.neighbour_receive_offsets.push_back(layout.receive_offsets[p]);
        }
        return layout;
    }
    throw std::invalid_argument("halo exchange supports at most three component layouts");
}

void HaloExchange::begin(double* values, std::size_t components, ExchangePlan& plan) {
    if (!m_parallel.distributed())
        return;
    if (values == nullptr || components == 0) {
        throw std::invalid_argument("halo exchange values are invalid");
    }
    if (plan.active)
        throw std::logic_error("halo exchange is already active");
    auto& layout = valueLayout(plan, components);
    const std::size_t send_size = plan.send_indices.size() * components;
    const std::size_t receive_size = plan.receive_indices.size() * components;
    if (plan.send_buffer.capacity() < send_size)
        plan.send_buffer.reserve(send_size);
    if (plan.receive_buffer.capacity() < receive_size)
        plan.receive_buffer.reserve(receive_size);
    plan.send_buffer.resize(send_size);
    plan.receive_buffer.resize(receive_size);
    for (std::size_t index = 0; index < plan.send_indices.size(); ++index) {
        const auto source =
            values + static_cast<std::size_t>(plan.send_indices[index]) * components;
        std::copy_n(source, components, plan.send_buffer.data() + index * components);
    }
    if (plan.graph_communicator != MPI_COMM_NULL) {
        detail::checkMpi(
            MPI_Ineighbor_alltoallv(
                plan.send_buffer.empty() ? &plan.dummy : plan.send_buffer.data(),
                layout.neighbour_send_counts.data(), layout.neighbour_send_offsets.data(),
                MPI_DOUBLE, plan.receive_buffer.empty() ? &plan.dummy : plan.receive_buffer.data(),
                layout.neighbour_receive_counts.data(), layout.neighbour_receive_offsets.data(),
                MPI_DOUBLE, plan.graph_communicator, &plan.request),
            "MPI_Ineighbor_alltoallv(halo values)");
    } else {
        detail::checkMpi(
            MPI_Ialltoallv(plan.send_buffer.empty() ? &plan.dummy : plan.send_buffer.data(),
                           layout.send_counts.data(), layout.send_offsets.data(), MPI_DOUBLE,
                           plan.receive_buffer.empty() ? &plan.dummy : plan.receive_buffer.data(),
                           layout.receive_counts.data(), layout.receive_offsets.data(), MPI_DOUBLE,
                           m_parallel.communicator, &plan.request),
            "MPI_Ialltoallv(halo values)");
    }
    plan.active = true;
    plan.active_components = components;
    plan.active_values = values;
}

void HaloExchange::finish(double* values, std::size_t components, ExchangePlan& plan) {
    if (!m_parallel.distributed())
        return;
    if (!plan.active || plan.active_values != values || plan.active_components != components) {
        throw std::logic_error("halo exchange finish does not match begin");
    }
    detail::checkMpi(MPI_Wait(&plan.request, MPI_STATUS_IGNORE), "MPI_Wait(halo values)");
    for (std::size_t index = 0; index < plan.receive_indices.size(); ++index) {
        const auto destination =
            values + static_cast<std::size_t>(plan.receive_indices[index]) * components;
        std::copy_n(plan.receive_buffer.data() + index * components, components, destination);
    }
    plan.active = false;
    plan.active_components = 0;
    plan.active_values = nullptr;
}

void HaloExchange::exchange(double* values, std::size_t components, ExchangePlan& plan) {
    begin(values, components, plan);
    finish(values, components, plan);
}

void HaloExchange::exchange(double* values, std::size_t components) {
    exchange(values, components, m_cells);
}

void HaloExchange::exchangeFaces(double* values, std::size_t components) {
    exchange(values, components, m_faces);
}

void HaloExchange::exchange(std::vector<double>& values) {
    if (m_mesh == nullptr || values.size() != static_cast<std::size_t>(m_mesh->cellCount())) {
        throw std::invalid_argument("raw halo field has the wrong size");
    }
    exchange(values.data(), 1);
}

void HaloExchange::exchangeFirstLayer(std::vector<double>& values) {
    if (m_mesh == nullptr || values.size() != static_cast<std::size_t>(m_mesh->cellCount())) {
        throw std::invalid_argument("raw first-layer halo field has the wrong size");
    }
    exchange(values.data(), 1, m_first_layer_cells);
}

void HaloExchange::beginFirstLayer(std::vector<double>& values) {
    if (m_mesh == nullptr || values.size() != static_cast<std::size_t>(m_mesh->cellCount())) {
        throw std::invalid_argument("raw first-layer halo field has the wrong size");
    }
    begin(values.data(), 1, m_first_layer_cells);
}

void HaloExchange::finishFirstLayer(std::vector<double>& values) {
    if (m_mesh == nullptr || values.size() != static_cast<std::size_t>(m_mesh->cellCount())) {
        throw std::invalid_argument("raw first-layer halo field has the wrong size");
    }
    finish(values.data(), 1, m_first_layer_cells);
}

std::size_t HaloExchange::plannedBytes(std::size_t components, bool first_layer, bool faces) const {
    if (components == 0)
        throw std::invalid_argument("halo component count is invalid");
    const ExchangePlan& plan = faces ? m_faces : (first_layer ? m_first_layer_cells : m_cells);
    const std::size_t values = plan.send_indices.size() + plan.receive_indices.size();
    if (values > (std::numeric_limits<std::size_t>::max() / components) / sizeof(double)) {
        throw std::overflow_error("halo payload size overflow");
    }
    return values * components * sizeof(double);
}

void HaloExchange::exchange(ScalarField& field) {
    if (&field.mesh() != m_mesh)
        throw std::invalid_argument("scalar halo field is incompatible");
    field.validateStorage();
    if (field.location() == FieldLocation::Cell)
        exchange(detail::fieldData(field), 1);
    else if (field.location() == FieldLocation::Face)
        exchangeFaces(detail::fieldData(field), 1);
    else
        throw std::invalid_argument("vertex scalar halo exchange is not supported");
}

void HaloExchange::exchange(VectorField& field) {
    if (&field.mesh() != m_mesh)
        throw std::invalid_argument("vector halo field is incompatible");
    field.validateStorage();
    if (field.location() == FieldLocation::Cell)
        exchange(&detail::fieldData(field)->x, 3);
    else if (field.location() == FieldLocation::Face)
        exchangeFaces(&detail::fieldData(field)->x, 3);
    else
        throw std::invalid_argument("vertex vector halo exchange is not supported");
}

void HaloExchange::exchange(TensorField& field) {
    if (&field.mesh() != m_mesh)
        throw std::invalid_argument("tensor halo field is incompatible");
    field.validateStorage();
    if (field.location() == FieldLocation::Cell)
        exchange(&detail::fieldData(field)->rows[0].x, 9);
    else if (field.location() == FieldLocation::Face)
        exchangeFaces(&detail::fieldData(field)->rows[0].x, 9);
    else
        throw std::invalid_argument("vertex tensor halo exchange is not supported");
}

} // namespace babelsim
