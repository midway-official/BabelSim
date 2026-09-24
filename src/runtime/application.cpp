#include "babelsim/case.h"
#include "babelsim/parallel.h"
#include "babelsim/application.h"
#include "internal/petsc_session.h"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace babelsim {

const SolverRegistration*& SolverRegistration::first() noexcept {
    // 函数内的零初始化头指针不依赖跨翻译单元的构造顺序。
    static const SolverRegistration* head = nullptr;
    return head;
}

SolverRegistration::SolverRegistration(const char* name, SolverResult (*run)(Case&)) noexcept
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
    std::filesystem::path performance_directory;
};

Arguments parseArguments(int argc, char* argv[]) {
    Arguments result;
    for (int index = 1; index < argc;) {
        const std::string option = argv[index++];
        if ((option != "-case" && option != "-time" && option != "-performance") || index == argc) {
            throw std::invalid_argument(
                "usage: babelsim-solve -case <case-directory> [-time <name>] "
                "[-performance <directory>]");
        }
        if (option == "-case") result.case_directory = argv[index];
        else if (option == "-time") result.time_name = argv[index];
        else result.performance_directory = argv[index];
        ++index;
    }
    if (result.case_directory.empty()) {
        throw std::invalid_argument("babelsim-solve needs -case <case-directory>");
    }
    return result;
}

const char* statusName(SolveStatus status) {
    switch (status) {
    case SolveStatus::Converged: return "converged";
    case SolveStatus::MaxIterations: return "maxIterations";
    case SolveStatus::NumericalFailure: return "numericalFailure";
    }
    return "unknown";
}

void writePerformance(
    const std::filesystem::path& directory,
    const Case& problem,
    const SolverResult& result,
    double mpi_init_seconds,
    double case_setup_seconds,
    double solver_seconds,
    double application_seconds)
{
    const ParallelContext parallel = ParallelContext::world();
    if (parallel.rank == 0) std::filesystem::create_directories(directory);
    parallel.barrier();
    const PerformanceCounters counters = problem.performance();
    const PerformanceCounters local = problem.localPerformance();
    const MeshPartitionInfo partition = problem.mesh().partitionInfo();
    const double solver_compute_seconds = std::max(
        0.0, solver_seconds - local.output_seconds);
    std::ostringstream name;
    name << "rank-" << std::setw(4) << std::setfill('0') << parallel.rank << ".json";
    std::ofstream output(directory / name.str());
    if (!output) throw std::runtime_error("cannot write performance report");
    output << std::setprecision(17)
           << "{\n  \"rank\": " << parallel.rank
           << ",\n  \"ranks\": " << parallel.size
           << ",\n  \"status\": \"" << statusName(result.status) << "\""
           << ",\n  \"mpiInitSeconds\": " << mpi_init_seconds
           << ",\n  \"caseSetupSeconds\": " << case_setup_seconds
           << ",\n  \"solverSeconds\": " << solver_seconds
           << ",\n  \"solverComputeSeconds\": " << solver_compute_seconds
           << ",\n  \"applicationSeconds\": " << application_seconds
           << ",\n  \"partition\": {"
           << "\n    \"globalCells\": " << partition.global_cells
           << ",\n    \"localCells\": " << partition.local_cells
           << ",\n    \"ownedCells\": " << partition.owned_cells
           << ",\n    \"ghostCells\": " << partition.ghost_cells
           << ",\n    \"localFaces\": " << partition.local_faces
           << ",\n    \"ownedFaces\": " << partition.owned_faces
           << ",\n    \"communicationFaces\": " << partition.communication_faces
           << ",\n    \"neighbourRanks\": " << partition.neighbour_ranks
           << "\n  }"
           << ",\n  \"linearSolves\": " << counters.linear_solves
           << ",\n  \"krylovIterations\": " << counters.krylov_iterations
           << ",\n  \"sparseMatvecs\": " << counters.sparse_matvecs
           << ",\n  \"haloExchanges\": " << counters.halo_exchanges
           << ",\n  \"haloBytes\": " << counters.halo_bytes
           << ",\n  \"globalReductions\": " << counters.global_reductions
           << ",\n  \"equationAssemblies\": " << counters.equation_assemblies
           << ",\n  \"matrixPatternBuilds\": " << counters.matrix_pattern_builds
           << ",\n  \"matrixValueUpdates\": " << counters.matrix_value_updates
           << ",\n  \"rhsOnlySolves\": " << counters.rhs_only_solves
           << ",\n  \"trueResidualChecks\": " << counters.true_residual_checks
           << ",\n  \"preconditionerSetups\": " << counters.preconditioner_setups
           << ",\n  \"preconditionerApplications\": " << counters.preconditioner_applications
           << ",\n  \"outputWrites\": " << counters.output_writes
           << ",\n  \"elapsedSeconds\": " << counters.elapsed_seconds
           << ",\n  \"assemblySeconds\": " << counters.assembly_seconds
           << ",\n  \"patternBuildSeconds\": " << counters.pattern_build_seconds
           << ",\n  \"matrixUpdateSeconds\": " << counters.matrix_update_seconds
           << ",\n  \"rhsSeconds\": " << counters.rhs_seconds
           << ",\n  \"residualCheckSeconds\": " << counters.residual_check_seconds
           << ",\n  \"preconditionerSeconds\": " << counters.preconditioner_seconds
           << ",\n  \"preconditionerApplySeconds\": " << counters.preconditioner_apply_seconds
           << ",\n  \"linearSolveSeconds\": " << counters.linear_solve_seconds
           << ",\n  \"sparseMatvecSeconds\": " << counters.sparse_matvec_seconds
           << ",\n  \"haloSeconds\": " << counters.halo_seconds
           << ",\n  \"globalReductionSeconds\": " << counters.global_reduction_seconds
           << ",\n  \"outputSeconds\": " << counters.output_seconds
           << ",\n  \"local\": {"
           << "\n    \"linearSolves\": " << local.linear_solves
           << ",\n    \"krylovIterations\": " << local.krylov_iterations
           << ",\n    \"sparseMatvecs\": " << local.sparse_matvecs
           << ",\n    \"haloExchanges\": " << local.halo_exchanges
           << ",\n    \"haloBytes\": " << local.halo_bytes
           << ",\n    \"globalReductions\": " << local.global_reductions
           << ",\n    \"equationAssemblies\": " << local.equation_assemblies
           << ",\n    \"matrixPatternBuilds\": " << local.matrix_pattern_builds
           << ",\n    \"matrixValueUpdates\": " << local.matrix_value_updates
           << ",\n    \"rhsOnlySolves\": " << local.rhs_only_solves
           << ",\n    \"trueResidualChecks\": " << local.true_residual_checks
           << ",\n    \"preconditionerSetups\": " << local.preconditioner_setups
           << ",\n    \"preconditionerApplications\": " << local.preconditioner_applications
           << ",\n    \"outputWrites\": " << local.output_writes
           << ",\n    \"elapsedSeconds\": " << local.elapsed_seconds
           << ",\n    \"assemblySeconds\": " << local.assembly_seconds
           << ",\n    \"patternBuildSeconds\": " << local.pattern_build_seconds
           << ",\n    \"matrixUpdateSeconds\": " << local.matrix_update_seconds
           << ",\n    \"rhsSeconds\": " << local.rhs_seconds
           << ",\n    \"residualCheckSeconds\": " << local.residual_check_seconds
           << ",\n    \"preconditionerSeconds\": " << local.preconditioner_seconds
           << ",\n    \"preconditionerApplySeconds\": " << local.preconditioner_apply_seconds
           << ",\n    \"linearSolveSeconds\": " << local.linear_solve_seconds
           << ",\n    \"sparseMatvecSeconds\": " << local.sparse_matvec_seconds
           << ",\n    \"haloSeconds\": " << local.halo_seconds
           << ",\n    \"globalReductionSeconds\": " << local.global_reduction_seconds
           << ",\n    \"outputSeconds\": " << local.output_seconds
           << ",\n    \"equationSystems\": [";
    for (std::size_t index = 0; index < local.equation_systems.size(); ++index) {
        const auto& equation = local.equation_systems[index];
        if (index != 0) output << ',';
        output << "\n      {\"identity\": \"" << equation.identity
               << "\", \"kspType\": \"" << equation.ksp_type
               << "\", \"pcType\": \"" << equation.pc_type
               << "\", \"lastStatus\": \"" << equation.last_status
               << "\", \"matrixPatternBuilds\": " << equation.matrix_pattern_builds
               << ", \"matrixValueUpdates\": " << equation.matrix_value_updates
               << ", \"rhsOnlySolves\": " << equation.rhs_only_solves
               << ", \"linearSolves\": " << equation.linear_solves
               << ", \"krylovIterations\": " << equation.krylov_iterations
               << ", \"preconditionerSetups\": " << equation.preconditioner_setups
               << ", \"trueResidualChecks\": " << equation.true_residual_checks
               << ", \"lastInitialResidual\": " << equation.last_initial_residual
               << ", \"lastFinalResidual\": " << equation.last_final_residual
               << ", \"lastRelativeResidual\": " << equation.last_relative_residual
               << ", \"patternBuildSeconds\": " << equation.pattern_build_seconds
               << ", \"matrixUpdateSeconds\": " << equation.matrix_update_seconds
               << ", \"preconditionerSeconds\": " << equation.preconditioner_seconds
               << ", \"rhsSeconds\": " << equation.rhs_seconds
               << ", \"solveSeconds\": " << equation.solve_seconds
               << ", \"residualCheckSeconds\": " << equation.residual_check_seconds << '}';
    }
    output << "\n    ]\n  }\n}\n";
    if (!output) throw std::runtime_error("cannot finish performance report");
}

}  // 匿名命名空间
}  // babelsim 命名空间

int babelsim::runApplication(int argc, char* argv[], ApplicationErrorHandler onError) {
    int initialized = 0;
    int finalized = 0;
    if (MPI_Initialized(&initialized) != MPI_SUCCESS ||
        MPI_Finalized(&finalized) != MPI_SUCCESS || finalized) {
        if (onError) onError("invalid MPI application lifecycle");
        return 1;
    }
    const bool owns_mpi = initialized == 0;
    const auto mpi_start = std::chrono::steady_clock::now();
    const int init_status = owns_mpi ? MPI_Init(&argc, &argv) : MPI_SUCCESS;
    if (init_status != MPI_SUCCESS) {
        if (onError) onError("MPI_Init failed");
        return 1;
    }
    const auto application_start = std::chrono::steady_clock::now();
    const double mpi_init_seconds = std::chrono::duration<double>(
        application_start - mpi_start).count();
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
        const auto case_start = std::chrono::steady_clock::now();
        Case problem(arguments.case_directory, arguments.time_name);
        const double case_setup_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - case_start).count();
        const SolverRegistration* selected = solvers;
        while (selected != nullptr && problem.solver() != selected->m_name) selected = selected->m_next;
        if (selected == nullptr) throw std::invalid_argument("unknown BabelSim solver: " + problem.solver());
        const auto solver_start = std::chrono::steady_clock::now();
        const SolverResult result = selected->m_run(problem);
        const double solver_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solver_start).count();
        if (!arguments.performance_directory.empty()) {
            const double application_seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - application_start).count();
            writePerformance(
                arguments.performance_directory, problem, result,
                mpi_init_seconds, case_setup_seconds, solver_seconds, application_seconds);
        }
        // Keep the established CLI contract: both nonconvergence and numerical
        // failure are exit 2; configuration/application errors remain exit 1.
        status = ParallelContext::world().maximum(
            result.status == SolveStatus::Converged ? 0 : 2);
    } catch (const std::exception& error) {
        if (onError) onError(error.what());
        // 单个 rank 的 I/O 失败不能让其他 rank 阻塞在后续 halo 交换或集体通信；
        // 正常的不收敛通过状态码 2 返回。
        const int abort_status = MPI_Abort(MPI_COMM_WORLD, 1);
        if (abort_status != MPI_SUCCESS) {
            if (onError) onError("MPI_Abort failed");
        }
        return 1;
    }
    // Finalize 不再放在可能抛异常的 try 块内；避免 finalize 失败后异常路径
    // 再次调用 MPI_Comm_rank/MPI_Abort，违反 MPI 生命周期。
    int petsc_finalize_status = 0;
    try {
        detail::finalizePetscSession();
    } catch (const std::exception& error) {
        petsc_finalize_status = 1;
        if (onError) onError(error.what());
    }
    if (petsc_finalize_status != 0) return 1;
    const int finalize_status = owns_mpi ? MPI_Finalize() : MPI_SUCCESS;
    if (finalize_status != MPI_SUCCESS) {
        if (onError) onError("MPI_Finalize failed");
        return 1;
    }
    return status;
}
