#include "babelsim/incompressible.h"
#include "babelsim/mesh_io.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace babelsim;

int main(int argc, char* argv[]) {
    try {
        if (argc < 2 || argc > 5) {
            throw std::invalid_argument(
                "usage: cavity <BabelSim-mesh> [max-iterations] [mu] [csv]");
        }
        const int maximum_iterations = argc > 2 ? std::stoi(argv[2]) : 5000;
        const double viscosity = argc > 3 ? std::stod(argv[3]) : 0.01;
        Mesh mesh = readMeshFile(argv[1]);
        if (maximum_iterations <= 0 || !(viscosity > 0.0)) {
            throw std::invalid_argument("maximum iterations and viscosity must be positive");
        }
        IncompressibleFields fields(mesh);
        for (Side side : {Side::XMin, Side::XMax, Side::YMin}) {
            fields.velocity.setBoundary(
                static_cast<Index>(side),
                BoundaryCondition<Vec3>::fixedValue({}));
        }
        fields.velocity.setBoundary(
            static_cast<Index>(Side::YMax),
            BoundaryCondition<Vec3>::fixedValue({1.0, 0.0, 0.0}));
        for (Side side : {Side::ZMin, Side::ZMax}) {
            fields.velocity.setBoundary(
                static_cast<Index>(side), BoundaryCondition<Vec3>::symmetry());
            fields.pressure.setBoundary(
                static_cast<Index>(side), BoundaryCondition<double>::symmetry());
        }

        Methods methods;
        methods.gradient = GradientMethod::GreenGauss;
        methods.convection = ConvectionMethod::Upwind;
        methods.diffusion = DiffusionMethod::Orthogonal;
        SimpleControl control;
        control.max_iterations = maximum_iterations;
        control.velocity_relaxation = 0.5;
        control.pressure_relaxation = 0.3;
        control.continuity_tolerance = 1e-8;
        control.velocity_tolerance = 1e-6;
        control.velocity_solver.absolute_tolerance = 1e-14;
        control.velocity_solver.relative_tolerance = 1e-8;
        control.pressure_solver.absolute_tolerance = 1e-14;
        control.pressure_solver.relative_tolerance = 1e-8;

        SimpleSolver solver(fields, {1.0, viscosity}, methods, control);
        SimpleIterationResult result;
        int completed = 0;
        for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
            result = solver.iterate();
            completed = iteration;
            if (iteration == 1 || iteration % 100 == 0 ||
                result.converged || !result.healthy) {
                std::cout << "SIMPLE " << std::setw(5) << iteration
                          << std::scientific << std::setprecision(4)
                          << " mass=" << result.continuity.relative
                          << " dU=" << result.relative_velocity_change
                          << " linP=" << result.pressure.relative_residual << '\n';
            }
            if (!result.healthy) {
                throw std::runtime_error("SIMPLE numerical failure");
            }
            if (result.converged) {
                break;
            }
        }

        const Index centre = mesh.cellId(
            mesh.dimensions[0] / 2 - 1, mesh.dimensions[1] / 2 - 1, 0);
        std::cout << "completed=" << completed
                  << " converged=" << std::boolalpha << result.converged
                  << " centre=" << fields.velocity[centre]
                  << " mass=" << result.continuity.relative
                  << " dU=" << result.relative_velocity_change << '\n';
        if (argc > 4) {
            std::ofstream output(argv[4]);
            if (!output) {
                throw std::runtime_error("cannot create output CSV");
            }
            output << "id,x,y,z,u,v,w,p\n" << std::setprecision(17);
            for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
                const Vec3& point = mesh.cell_centres[static_cast<std::size_t>(cell)];
                const Vec3& velocity = fields.velocity[cell];
                output << cell << ',' << point.x << ',' << point.y << ',' << point.z
                       << ',' << velocity.x << ',' << velocity.y << ',' << velocity.z
                       << ',' << fields.pressure[cell] << '\n';
            }
        }
        return result.converged ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
