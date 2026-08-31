#include "babelsim/case.h"
#include "babelsim/mpi_support.h"
#include "test_util.h"
#include <iostream>

namespace babelsim { int runCoupledScalar(Case&); }

// 测试启动器，不是用户编写的 Solver。实际应用由 babelsim-solve 承担这一职责。
int main(int argc, char** argv) {
    babelsim::detail::checkMpi(MPI_Init(&argc, &argv), "MPI_Init");
    int result = 1;
    try {
        if (argc != 3) throw std::invalid_argument("expected case directory and run name");
        babelsim::Case problem(argv[1], argv[2]);
        babelsim::ScalarField& first = problem.scalarField("T");
        for (int i = 0; i < 24; ++i)
            problem.vectorField("scratch" + std::to_string(i), babelsim::Vec3{});
        require(&first == &problem.scalarField("T"), "Case invalidated a Field reference");
        const babelsim::TensorField& stress = problem.tensorField("stress");
        require(near(stress[0][2][2], 9.0), "tensor input order is incorrect");
        bool rejected = false;
        try { problem.vectorField("T", babelsim::Vec3{}); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Case accepted one name for different Field types");
        result = babelsim::runCoupledScalar(problem);
        if (result == 0) problem.finish();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        babelsim::detail::checkMpi(MPI_Abort(MPI_COMM_WORLD, 1), "MPI_Abort");
    }
    babelsim::detail::checkMpi(MPI_Finalize(), "MPI_Finalize");
    return result;
}
