#include "babelsim/case.h"
#include "babelsim/parallel.h"
#include "babelsim/application.h"

#include <mpi.h>

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace babelsim {
namespace {

struct Arguments {
    std::filesystem::path case_directory;
    std::string time_name;
};

Arguments parseArguments(int argc, char* argv[]) {
    Arguments result;
    for (int index = 1; index < argc;) {
        const std::string option = argv[index++];
        if ((option != "-case" && option != "-time") || index == argc) {
            throw std::invalid_argument(
                "usage: babelsim-solve -case <case-directory> [-time <name>]");
        }
        if (option == "-case") result.case_directory = argv[index];
        else result.time_name = argv[index];
        ++index;
    }
    if (result.case_directory.empty()) {
        throw std::invalid_argument("babelsim-solve needs -case <case-directory>");
    }
    return result;
}

int run(const Arguments& arguments, const SolverEntry* solvers, std::size_t count) {
    if (solvers == nullptr || count == 0) throw std::invalid_argument("empty solver table");
    Case problem(arguments.case_directory, arguments.time_name);
    int status = 1;
    bool found = false;
    for (std::size_t index = 0; index < count; ++index) {
        if (solvers[index].name == nullptr || solvers[index].run == nullptr)
            throw std::invalid_argument("invalid solver entry");
        if (solvers[index].name[0] == '\0') throw std::invalid_argument("empty solver name");
        for (std::size_t previous = 0; previous < index; ++previous)
            if (std::string(solvers[previous].name) == solvers[index].name)
                throw std::invalid_argument("duplicate solver entry");
        if (problem.solver() != solvers[index].name) continue;
        if (found) throw std::invalid_argument("duplicate solver entry");
        found = true;
    }
    if (!found) throw std::invalid_argument("unknown BabelSim solver: " + problem.solver());
    for (std::size_t index = 0; index < count; ++index)
        if (problem.solver() == solvers[index].name) { status = solvers[index].run(problem); break; }
    // 用户函数的任意非零值均是失败；不能让负返回码在全局 maximum 中被 0 掩盖。
    status = ParallelContext::world().maximum(status < 0 ? 1 : status);
    if (status == 0) problem.finish();
    return status;
}

}  // 匿名命名空间
}  // babelsim 命名空间

int babelsim::runApplication(int argc, char* argv[], const SolverEntry* solvers, std::size_t count) {
    int initialized = 0;
    int finalized = 0;
    if (MPI_Initialized(&initialized) != MPI_SUCCESS ||
        MPI_Finalized(&finalized) != MPI_SUCCESS || finalized) {
        std::cerr << "invalid MPI application lifecycle\n";
        return 1;
    }
    const bool owns_mpi = initialized == 0;
    const int init_status = owns_mpi ? MPI_Init(&argc, &argv) : MPI_SUCCESS;
    if (init_status != MPI_SUCCESS) {
        std::cerr << "MPI_Init failed with code " << init_status << '\n';
        return 1;
    }
    int status = 1;
    try {
        status = babelsim::run(babelsim::parseArguments(argc, argv), solvers, count);
    } catch (const std::exception& error) {
        int rank = 0;
        const int rank_status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank_status != MPI_SUCCESS) rank = -1;
        std::cerr << "babelsim-solve rank " << rank << ": " << error.what() << '\n';
        // 单个 rank 的 I/O 失败不能让其他 rank 阻塞在后续 halo 交换或集体通信；
        // 正常的不收敛通过状态码 2 返回。
        const int abort_status = MPI_Abort(MPI_COMM_WORLD, 1);
        if (abort_status != MPI_SUCCESS) {
            std::cerr << "MPI_Abort failed with code " << abort_status << '\n';
        }
        return 1;
    }
    // Finalize 不再放在可能抛异常的 try 块内；避免 finalize 失败后异常路径
    // 再次调用 MPI_Comm_rank/MPI_Abort，违反 MPI 生命周期。
    const int finalize_status = owns_mpi ? MPI_Finalize() : MPI_SUCCESS;
    if (finalize_status != MPI_SUCCESS) {
        std::cerr << "MPI_Finalize failed with code " << finalize_status << '\n';
        return 1;
    }
    return status;
}
