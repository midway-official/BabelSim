#include "babelsim/case.h"
#include "babelsim/parallel.h"
#include "babelsim/application.h"

#include <mpi.h>

#include <cstring>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace babelsim {

const SolverRegistration*& SolverRegistration::first() noexcept {
    // 函数内的零初始化头指针不依赖跨翻译单元的构造顺序。
    static const SolverRegistration* head = nullptr;
    return head;
}

SolverRegistration::SolverRegistration(const char* name, int (*run)(Case&)) noexcept
    : m_name(name), m_run(run), m_next(first())
{
    first() = this;
}

SolverRegistration::~SolverRegistration() noexcept {
    // 既支持正常静态析构，也避免临时注册对象离开作用域后留下悬空描述项。
    const SolverRegistration** link = &first();
    while (*link != nullptr && *link != this) link = &(*link)->m_next;
    if (*link != nullptr) *link = m_next;
}

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

}  // 匿名命名空间
}  // babelsim 命名空间

int babelsim::runApplication(int argc, char* argv[]) {
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
        const Arguments arguments = parseArguments(argc, argv);
        const SolverRegistration* solvers = SolverRegistration::first();
        if (solvers == nullptr) throw std::invalid_argument("no registered solvers");
        for (const auto* entry = solvers; entry != nullptr; entry = entry->m_next) {
            if (entry->m_name == nullptr || entry->m_run == nullptr)
                throw std::invalid_argument("invalid solver registration");
            if (entry->m_name[0] == '\0') throw std::invalid_argument("empty solver name");
            for (const auto* previous = solvers; previous != entry; previous = previous->m_next)
                if (std::strcmp(previous->m_name, entry->m_name) == 0)
                    throw std::invalid_argument("duplicate solver registration: " + std::string(entry->m_name));
        }
        Case problem(arguments.case_directory, arguments.time_name);
        const SolverRegistration* selected = solvers;
        while (selected != nullptr && problem.solver() != selected->m_name) selected = selected->m_next;
        if (selected == nullptr) throw std::invalid_argument("unknown BabelSim solver: " + problem.solver());
        status = selected->m_run(problem);
        // 任意非零值均是失败；负返回码不能在全局 maximum 中被 0 掩盖。
        status = ParallelContext::world().maximum(status < 0 ? 1 : status);
        if (status == 0) problem.finish();
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
