#include "babelsim/parallel.h"
#include "babelsim/mpi_support.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace babelsim {
namespace {

void requireCommunicator(MPI_Comm communicator, const char* operation) {
    if (communicator == MPI_COMM_NULL) {
        throw std::invalid_argument(std::string(operation) + " received MPI_COMM_NULL");
    }
    detail::requireMpiReady(operation);
}

void requireCollectiveContext(const ParallelContext& parallel, const char* operation) {
    if (parallel.size <= 0 || parallel.rank < 0 || parallel.rank >= parallel.size) {
        throw std::invalid_argument("parallel rank or size is invalid");
    }
    if (parallel.communicator == MPI_COMM_NULL) {
        if (parallel.size != 1 || parallel.rank != 0) {
            throw std::invalid_argument("distributed collective has MPI_COMM_NULL");
        }
        return;
    }
    detail::requireMpiReady(operation);
}

} // namespace

ParallelContext ParallelContext::world(MPI_Comm communicator_value) {
    requireCommunicator(communicator_value, "ParallelContext::world");
    ParallelContext result;
    result.communicator = communicator_value;
    detail::checkMpi(MPI_Comm_rank(communicator_value, &result.rank), "MPI_Comm_rank");
    detail::checkMpi(MPI_Comm_size(communicator_value, &result.size), "MPI_Comm_size");
    result.validate();
    return result;
}

void ParallelContext::validate() const {
    if (size <= 0 || rank < 0 || rank >= size) {
        throw std::invalid_argument("parallel communicator, rank, or size is invalid");
    }
    if (communicator == MPI_COMM_NULL) {
        if (size != 1 || rank != 0) {
            throw std::invalid_argument("distributed context cannot use MPI_COMM_NULL");
        }
        return;
    }
    requireCommunicator(communicator, "ParallelContext::validate");
    int actual_rank = -1;
    int actual_size = 0;
    detail::checkMpi(MPI_Comm_rank(communicator, &actual_rank), "MPI_Comm_rank");
    detail::checkMpi(MPI_Comm_size(communicator, &actual_size), "MPI_Comm_size");
    if (actual_rank != rank || actual_size != size) {
        throw std::logic_error("parallel context rank/size no longer matches communicator");
    }
}

void ParallelContext::sum(const double* local, double* global, int count) const {
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel sum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::sum");
    if (count == 0)
        return;
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(MPI_Allreduce(local, global, count, MPI_DOUBLE, MPI_SUM, communicator),
                     "MPI_Allreduce(sum)");
}

void ParallelContext::sum(const int* local, int* global, int count) const {
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel integer sum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::sum int array");
    if (count == 0)
        return;
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(MPI_Allreduce(local, global, count, MPI_INT, MPI_SUM, communicator),
                     "MPI_Allreduce(sum int array)");
}

void ParallelContext::maximum(const double* local, double* global, int count) const {
    if (count < 0 || (count > 0 && (local == nullptr || global == nullptr))) {
        throw std::invalid_argument("parallel maximum buffer is invalid");
    }
    requireCollectiveContext(*this, "ParallelContext::maximum");
    if (count == 0)
        return;
    if (!distributed()) {
        std::copy(local, local + count, global);
        return;
    }
    detail::checkMpi(MPI_Allreduce(local, global, count, MPI_DOUBLE, MPI_MAX, communicator),
                     "MPI_Allreduce(maximum)");
}

int ParallelContext::sum(int local) const {
    int global = 0;
    sum(&local, &global, 1);
    return global;
}

int ParallelContext::maximum(int local) const {
    requireCollectiveContext(*this, "ParallelContext::maximum int");
    if (!distributed())
        return local;
    int global = 0;
    detail::checkMpi(MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, communicator),
                     "MPI_Allreduce(maximum int)");
    return global;
}

void ParallelContext::barrier() const {
    validate();
    if (distributed())
        detail::checkMpi(MPI_Barrier(communicator), "MPI_Barrier");
}

} // namespace babelsim
