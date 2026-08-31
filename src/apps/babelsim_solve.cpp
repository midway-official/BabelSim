#include "babelsim/case.h"
#include "babelsim/parallel.h"
#include "solver_selection.h"

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

int run(const Arguments& arguments) {
    Case problem(arguments.case_directory, arguments.time_name);
    const int status = runSolver(problem);
    if (status == 0) problem.finish();
    return status;
}

}  // 匿名命名空间
}  // babelsim 命名空间

int main(int argc, char* argv[]) {
    const int init_status = MPI_Init(&argc, &argv);
    if (init_status != MPI_SUCCESS) {
        std::cerr << "MPI_Init failed with code " << init_status << '\n';
        return 1;
    }
    int status = 1;
    try {
        status = babelsim::run(babelsim::parseArguments(argc, argv));
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
    const int finalize_status = MPI_Finalize();
    if (finalize_status != MPI_SUCCESS) {
        std::cerr << "MPI_Finalize failed with code " << finalize_status << '\n';
        return 1;
    }
    return status;
}
