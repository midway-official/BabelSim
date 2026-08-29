#include "babelsim/case.h"
#include "babelsim/field_io.h"
#include "babelsim/incompressible_io.h"
#include "babelsim/mesh_io.h"
#include "babelsim/parallel.h"
#include "babelsim/parallel_writer.h"

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

template <typename T>
void localize(const Field<T>& global, Field<T>& local) {
    copyBoundaryConditions(global, local);
    for (Index cell = 0; cell < local.mesh().cellCount(); ++cell) {
        local[cell] = global[local.mesh().globalCellId(cell)];
    }
}

int runSimpleFoam(
    const CaseDefinition& definition,
    const ParallelContext& parallel,
    const std::string& requested_time)
{
    const Mesh global_mesh = readMeshFile(definition.mesh_file);
    IncompressibleFields global_fields(global_mesh);
    readFieldFile(definition.fields_directory / "U.field", global_fields.velocity);
    readFieldFile(definition.fields_directory / "p.field", global_fields.pressure);

    Mesh local_mesh = decompose(global_mesh, parallel);
    IncompressibleFields fields(local_mesh);
    localize(global_fields.velocity, fields.velocity);
    localize(global_fields.pressure, fields.pressure);
    fields.face_flux.fill(0.0);

    const IncompressibleCaseControl controls = readIncompressibleCase(definition);
    SimpleSolver solver(fields, controls.fluid, controls.methods, controls.simple, parallel);
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
                      << " linP=" << result.pressure.relative_residual << '\n';
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
                  << " mass=" << std::scientific << result.continuity.relative
                  << " dU=" << result.relative_velocity_change << '\n';
    }
    return result.converged ? 0 : 2;
}

int run(const Arguments& arguments, const ParallelContext& parallel) {
    const CaseDefinition definition = readCase(arguments.case_directory);
    if (definition.solver == "simpleFoam") {
        return runSimpleFoam(definition, parallel, arguments.time_name);
    }
    throw std::runtime_error("unknown BabelSim solver: " + definition.solver);
}

}  // 匿名命名空间
}  // babelsim 命名空间

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    try {
        const babelsim::ParallelContext parallel = babelsim::ParallelContext::world();
        const int status = babelsim::run(babelsim::parseArguments(argc, argv), parallel);
        MPI_Finalize();
        return status;
    } catch (const std::exception& error) {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cerr << "babelsim-solve rank " << rank << ": " << error.what() << '\n';
        // 单个 rank 的 I/O 失败不能让其他 rank 阻塞在后续 halo 交换或集体通信；
        // 正常的不收敛通过状态码 2 返回。
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return 1;
}
