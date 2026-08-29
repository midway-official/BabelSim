#pragma once

#include <mpi.h>

#include <limits>
#include <stdexcept>
#include <string>

namespace babelsim::detail {

// MPI 函数即使在极少数实现中返回错误，也不能继续使用未定义的通信状态。
inline void checkMpi(int code, const char* operation) {
    if (code == MPI_SUCCESS) return;
    char message[MPI_MAX_ERROR_STRING]{};
    int length = 0;
    const int string_code = MPI_Error_string(code, message, &length);
    if (string_code != MPI_SUCCESS || length < 0 ||
        length >= MPI_MAX_ERROR_STRING) {
        throw std::runtime_error(
            std::string(operation) + " failed with MPI error code " +
            std::to_string(code));
    }
    throw std::runtime_error(
        std::string(operation) + " failed: " + std::string(message, length));
}

inline void requireMpiReady(const char* operation) {
    int initialized = 0;
    checkMpi(MPI_Initialized(&initialized), "MPI_Initialized");
    if (initialized == 0) {
        throw std::logic_error(
            std::string(operation) + " requires MPI_Init/MPI_Init_thread");
    }
    int finalized = 0;
    checkMpi(MPI_Finalized(&finalized), "MPI_Finalized");
    if (finalized != 0) {
        throw std::logic_error(std::string(operation) + " called after MPI_Finalize");
    }
}

inline int mpiCount(std::size_t count, const char* what) {
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error(std::string(what) + " exceeds the MPI int count limit");
    }
    return static_cast<int>(count);
}

}  // namespace babelsim::detail
