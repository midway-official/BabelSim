#include "babelsim/case.h"
#include "babelsim/numerics_io.h"
#include "babelsim/field_io.h"
#include "babelsim/mesh_io.h"

#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    const CaseDefinition cavity = readCase("cases/cavity");
    require(cavity.solver == "simple", "cavity solver selection is incorrect");
    const Parameters properties(cavity.physics_file);
    const Parameters solution(cavity.solution_file);
    require(near(properties.number("density"), 1.0), "cavity density is incorrect");
    require(near(properties.number("dynamicViscosity"), 0.01), "cavity viscosity is incorrect");
    require(solution.integer("maxIterations", 0) == 5000, "cavity iteration limit is incorrect");

    const Mesh mesh = readMeshFile(cavity.mesh_file);
    ScalarField pressure(mesh, FieldLocation::Cell, "p");
    VectorField velocity(mesh, FieldLocation::Cell, "U");
    readFieldFile(cavity.fields_directory / "U.field", velocity);
    readFieldFile(cavity.fields_directory / "p.field", pressure);
    require(velocity.boundary(static_cast<Index>(Side::YMax)).type == BoundaryType::FixedValue,
            "lid velocity condition was not read");
    require(near(velocity.boundary(static_cast<Index>(Side::YMax)).value, {1.0, 0.0, 0.0}),
            "lid velocity value was not read");
    require(pressure.boundary(static_cast<Index>(Side::ZMin)).type == BoundaryType::Symmetry,
            "scalar symmetry condition was not read");

    const CaseDefinition channel = readCase("cases/poiseuille");
    const Mesh channel_mesh = readMeshFile(channel.mesh_file);
    VectorField channel_velocity(channel_mesh, FieldLocation::Cell, "U");
    readFieldFile(channel.fields_directory / "U.field", channel_velocity);
    require(channel_velocity.boundary(static_cast<Index>(Side::XMin)).type == BoundaryType::FixedValue,
            "inlet condition was not read");
    require(near(channel_velocity.boundary(static_cast<Index>(Side::XMin)).value, {1.0, 0.0, 0.0}),
            "inlet velocity value was not read");

    const CaseDefinition transport = readCase("cases/transport");
    const Parameters transport_properties(transport.physics_file);
    const Methods methods = readMethodsFile(transport.methods_file);
    require(transport.solver == "transport", "transport solver selection is incorrect");
    require(near(transport_properties.number("diffusivity"), 0.01), "transport diffusivity is incorrect");
    require(methods.convectionFor("C") == ConvectionMethod::Upwind,
            "Field-specific convection method was not read");
    std::cout << "case_io_test: SIMPLE, heat-compatible and transport dictionaries passed\n";
}
