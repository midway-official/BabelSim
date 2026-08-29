#include "babelsim/incompressible.h"
#include "babelsim/legacy_taiho.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace babelsim;

int main(int argc, char* argv[]) {
    try {
        if (argc < 2 || argc > 6) {
            throw std::invalid_argument(
                "usage: imported_simple <Taiho-mesh> [max-iterations] [rho] [mu] [csv]");
        }
        const int maximum_iterations = argc > 2 ? std::stoi(argv[2]) : 5000;
        const double density = argc > 3 ? std::stod(argv[3]) : 1.0;
        const double viscosity = argc > 4 ? std::stod(argv[4]) : 0.01;
        ImportedTaihoMesh imported = readTaihoMesh(argv[1]);
        IncompressibleFields fields(imported.mesh);
        applyImportedBoundaryConditions(imported, fields.velocity, fields.pressure);

        Methods methods;
        methods.gradient = GradientMethod::GreenGauss;
        methods.convection = ConvectionMethod::Upwind;
        methods.diffusion = DiffusionMethod::Orthogonal;
        SimpleControl control;
        control.max_iterations = maximum_iterations;
        control.velocity_relaxation = 0.5;
        control.pressure_relaxation = 0.3;
        control.continuity_tolerance = 1e-7;
        control.velocity_tolerance = 1e-6;
        control.velocity_solver.absolute_tolerance = 1e-14;
        control.velocity_solver.relative_tolerance = 1e-7;
        control.velocity_solver.max_iterations = 200;
        control.pressure_solver.absolute_tolerance = 1e-14;
        control.pressure_solver.relative_tolerance = 1e-7;
        control.pressure_solver.max_iterations = 1000;

        SimpleSolver solver(fields, {density, viscosity}, methods, control);
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

        Vec3 sum_velocity{};
        double sum_pressure = 0.0;
        double maximum_velocity = 0.0;
        for (Index cell = 0; cell < imported.mesh.cellCount(); ++cell) {
            sum_velocity += fields.velocity[cell];
            sum_pressure += fields.pressure[cell];
            maximum_velocity = std::max(maximum_velocity, norm(fields.velocity[cell]));
        }
        std::cout << "completed=" << completed
                  << " converged=" << std::boolalpha << result.converged
                  << " cells=" << imported.mesh.cellCount()
                  << " sumU=" << sum_velocity
                  << " sumP=" << sum_pressure
                  << " max|U|=" << maximum_velocity
                  << " mass=" << result.continuity.relative
                  << " dU=" << result.relative_velocity_change << '\n';

        if (argc == 6) {
            std::ofstream output(argv[5]);
            if (!output) {
                throw std::runtime_error("cannot create output CSV");
            }
            output << "id,x,y,z,u,v,w,p\n" << std::setprecision(17);
            for (Index cell = 0; cell < imported.mesh.cellCount(); ++cell) {
                const Vec3& centre = imported.mesh.cell_centres[static_cast<std::size_t>(cell)];
                const Vec3& velocity = fields.velocity[cell];
                output << cell << ',' << centre.x << ',' << centre.y << ',' << centre.z
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

