#include "babelsim/case.h"
#include "babelsim/field_io.h"
#include "babelsim/incompressible_io.h"
#include "babelsim/mesh_io.h"
#include "babelsim/parallel.h"
#include "babelsim/parallel_writer.h"
#include "babelsim/thermal.h"
#include "babelsim/thermal_io.h"

#include <mpi.h>

#include <filesystem>
#include <iomanip>
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

int runSimpleFoam(
    const CaseDefinition& definition,
    const ParallelContext& parallel,
    const std::string& requested_time)
{
    // 根 rank 读取并分发局部几何；每个 rank 只构造自己的 owned+ghost Mesh。
    // Field 同样直接绑定局部 Mesh，避免先创建全局 Field 再复制内部值。
    Mesh local_mesh = readDistributedMesh(definition.mesh_file, parallel);
    IncompressibleFields fields(local_mesh);
    readFieldFile(definition.fields_directory / "U.field", fields.velocity);
    readFieldFile(definition.fields_directory / "p.field", fields.pressure);
    fields.face_flux.fill(0.0);

    const IncompressibleCaseControl controls = readIncompressibleCase(definition);
    RunTime run_time = RunTime::forMesh(
        local_mesh, simpleRunTimeControl(controls.methods, controls.simple));
    SimpleSolver solver(run_time, fields, controls.fluid, controls.simple);
    SimpleIterationResult result;
    int completed = 0;
    for (int iteration = 1; iteration <= controls.simple.max_iterations; ++iteration) {
        result = solver.iterate();
        completed = iteration;
        if (parallel.rank == 0 &&
            (iteration == 1 || iteration % 100 == 0 || result.converged || !result.healthy)) {
            std::cout << "SIMPLE " << std::setw(5) << iteration << std::scientific
                      << std::setprecision(4) << " mass=" << result.continuity.relative
                      << " dU=" << result.relative_velocity_change
                      << " linP=" << result.pressure.relative_residual
                      << " linear=" << (result.linear_converged ? "ok" : "inexact")
                      << '\n';
        }
        if (!result.healthy || result.converged) break;
    }
    if (!result.healthy) {
        throw std::runtime_error("SIMPLE encountered a numerical failure");
    }

    OutputControl output = readOutputControl(definition);
    if (!requested_time.empty()) output.time_name = requested_time;
    const auto time_directory = definition.root / output.directory / output.time_name;
    writeOwnedFieldCsv(time_directory, fields.velocity, parallel);
    writeOwnedFieldCsv(time_directory, fields.pressure, parallel);
    writeOwnedResultMetadata(
        time_directory, local_mesh, parallel, output.time_name,
        {{fields.velocity.name(), "vector", fields.velocity.location()},
         {fields.pressure.name(), "scalar", fields.pressure.location()}});

    if (parallel.rank == 0) {
        std::cout << "SIMPLE completed=" << completed
                  << " converged=" << (result.converged ? "true" : "false")
                  << " linear=" << (result.linear_converged ? "ok" : "inexact")
                  << " mass=" << std::scientific << result.continuity.relative
                  << " dU=" << result.relative_velocity_change << '\n';
    }
    return result.converged ? 0 : 2;
}

int runHeatFoam(
    const CaseDefinition& definition,
    const ParallelContext& parallel,
    const std::string& requested_time)
{
    Mesh local_mesh = readDistributedMesh(definition.mesh_file, parallel);
    ScalarField temperature(local_mesh, FieldLocation::Cell, "T");
    readFieldFile(definition.fields_directory / "T.field", temperature);

    const ThermalCaseControl controls = readThermalCase(definition);
    RunTime run_time = RunTime::forMesh(local_mesh, controls.runtime);
    const HeatResult result = solveTransientHeat(
        run_time, temperature, controls.material, controls.volumetric_source);
    if (!result.converged) {
        throw std::runtime_error("heat solver did not converge");
    }

    OutputControl output = readOutputControl(definition);
    if (!requested_time.empty()) output.time_name = requested_time;
    const auto time_directory = definition.root / output.directory / output.time_name;
    writeOwnedFieldCsv(time_directory, temperature, parallel);
    writeOwnedResultMetadata(
        time_directory, local_mesh, parallel, output.time_name,
        {{temperature.name(), "scalar", temperature.location()}});

    if (parallel.rank == 0) {
        std::cout << "Heat completed steps=" << result.steps << std::scientific
                  << " residual=" << result.linear.relative_residual << '\n';
    }
    return 0;
}

int run(const Arguments& arguments, const ParallelContext& parallel) {
    const CaseDefinition definition = readCase(arguments.case_directory);
    if (definition.solver == "heatFoam") {
        return runHeatFoam(definition, parallel, arguments.time_name);
    }
    if (definition.solver == "simpleFoam") {
        return runSimpleFoam(definition, parallel, arguments.time_name);
    }
    throw std::runtime_error("unknown BabelSim solver: " + definition.solver);
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
        const babelsim::ParallelContext parallel = babelsim::ParallelContext::world();
        status = babelsim::run(babelsim::parseArguments(argc, argv), parallel);
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
