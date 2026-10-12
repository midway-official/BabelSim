#include "babelsim/case.h"
#include "internal/field_access.h"
#include "support/simple_reference.h"
#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    {
        Case problem("cases/heat", "lifecycle");
        const Parameters& physics = problem.read.physics();
        require(&physics == &problem.read.physics(), "Case recreated its physics dictionary");
        require(near(physics.positive("density"), 1.0) &&
                near(physics.positive("heatCapacity"), 1.0) &&
                near(physics.nonnegative("conductivity"), 0.1) &&
                near(physics.number("source"), 0.0),
                "Case physics values do not match the selected dictionary");
        ScalarField& first = problem.load.scalarField("T");
        readEquationControl(problem, "temperature", first);
        for (int i = 0; i < 64; ++i)
            problem.create.vectorField("scratch" + std::to_string(i), Vec3{});
        require(&first == &problem.load.scalarField("T"), "Case invalidated a Field reference");
        bool rejected = false;
        try { problem.create.vectorField("T", Vec3{}); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Case accepted one name for different Field types");
        TensorField& tensor = problem.create.tensorField("gradient", Tensor3{});
        ScalarField& initialized = problem.create.scalarField("initialized", 1.0);
        bool initialization_rejected = false;
        try { problem.create.scalarField("initialized", 2.0); }
        catch (const std::logic_error&) { initialization_rejected = true; }
        require(initialization_rejected,
                "a repeated create call silently ignored its initialization value");
        bool load_create_rejected = false;
        try { problem.load.scalarField("initialized"); }
        catch (const std::logic_error&) { load_create_rejected = true; }
        require(load_create_rejected,
                "a created field was silently treated as a file-loaded field");
        require(near(detail::fieldData(initialized)[0], 1.0),
                "created field initialization was not applied");
        problem.output(tensor);
        problem.output(first, false);
        auto& faceVector = problem.create.faceVectorField("faceVector", {1.0, 2.0, 3.0});
        auto& faceTensor = problem.create.faceTensorField("faceTensor");
        auto& faceScalar = problem.create.faceScalarField("faceScalar", 2.5);
        require(near(detail::fieldData(faceVector)[0], {1.0, 2.0, 3.0}) &&
                near(detail::fieldData(faceScalar)[0], 2.5), "face initial values were ignored");
        require(&problem.fields.faceVectorField("faceVector") == &faceVector &&
                &problem.fields.faceTensorField("faceTensor") == &faceTensor &&
                &problem.fields.faceScalarField("faceScalar") == &faceScalar,
                "face lookup did not return the declared field");
        rejected = false;
        try { problem.output(faceScalar); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "unsupported face output was silently accepted");
        ScalarField copy(first);
        rejected = false;
        try { problem.output(copy); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "non-owned Field was registered for Case output");

        // A correction is a normal named field: creation, equation binding,
        // lookup and output all refer to the same object and preserve the input.
        const ScalarField inputBefore = first;
        auto& correction = problem.create.homogeneousCorrection(first, "TPrime");
        for (int i = 0; i < 64; ++i)
            problem.create.scalarField("coefficient" + std::to_string(i), 1.0);
        require(&correction == &problem.fields.scalarField("TPrime"),
                "registering more fields invalidated the correction reference");
        require(near(math::normL2(correction), 0.0), "correction copied input values");
        for (Index patch = 0; patch < problem.mesh().patchCount(); ++patch) {
            const auto kind = first.boundary(patch).type;
            const auto expected = kind == BoundaryType::FixedValue ? BoundaryType::FixedValue
                : kind == BoundaryType::Symmetry ? BoundaryType::Symmetry : BoundaryType::ZeroGradient;
            require(correction.boundary(patch).type == expected &&
                    near(correction.boundary(patch).value, 0.0), "correction boundary mapping changed");
        }
        auto equation = testEquation(correction);
        equ::laplacian(equation, 1.0, -1);
        correction.fill(3.0);
        require(equ::solve(equation).converged() &&
                math::normL2(problem.fields.scalarField("TPrime")) < 1e-10,
                "equation solve did not update the registered correction");
        require(near(math::normL2(first - inputBefore), 0.0), "correction changed its input field");
        require(near(first.boundary(0).value, inputBefore.boundary(0).value),
                "correction changed input boundary values");
        bool duplicate_rejected = false;
        try { problem.create.homogeneousCorrection(first, "TPrime"); }
        catch (const std::logic_error&) { duplicate_rejected = true; }
        require(duplicate_rejected, "correction creation silently reset an existing field");
        bool face_rejected = false;
        try { problem.create.homogeneousCorrection(faceScalar, "invalidCorrection"); }
        catch (const std::invalid_argument&) { face_rejected = true; }
        require(face_rejected, "correction accepted a face field");
        // Failed creation must not reserve a name or leave a partially declared field.
        auto& afterFailure = problem.create.scalarField("invalidCorrection", 4.0);
        require(near(detail::fieldData(afterFailure)[0], 4.0), "failed creation reserved its name");
        const auto otherMesh = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
        ScalarField foreign(otherMesh, FieldLocation::Cell);
        bool mesh_rejected = false;
        try { problem.create.homogeneousCorrection(foreign, "foreignCorrection"); }
        catch (const std::invalid_argument&) { mesh_rejected = true; }
        require(mesh_rejected, "correction accepted a foreign mesh");
        problem.output(correction);
        problem.finalizeInitialization();
        problem.write();
        const auto written = problem.definition().root / problem.read.outputControl().directory
            / "lifecycle" / "rank-0000" / "TPrime.csv";
        require(std::filesystem::is_regular_file(written), "registered correction was not written");
    }
    {
        Case problem("cases/cavity", "lifecycle-simple");
        problem.load.vectorField("U"); problem.load.scalarField("p"); problem.create.faceScalarField("phi");
        problem.read.physics().positive("density");
        problem.read.physics().positive("dynamicViscosity");
        readSimpleControl(problem.read.solution());
        readEquationControl(problem, "momentum", problem.load.vectorField("U"));
        readEquationControl(problem, "pressureCorrection", problem.load.scalarField("p"));
        problem.create.scalarField("couplingState", 0.0);
        problem.finalizeInitialization();
        problem.finalizeInitialization(); // idempotent, and no time advancement
        bool rejected = false;
        try { problem.create.scalarField("tooLate", 0.0); }
        catch (const std::logic_error&) { rejected = true; }
        require(rejected, "Case start did not close declarations");
    }
    // 前一个 Case 必须先完整释放 RunTime，才能在同一线程创建下一次计算。
    Case second("cases/heat", "lifecycle2");
    require(second.load.scalarField("T").size() == static_cast<std::size_t>(second.mesh().cellCount()),
            "Case recreated a Field with the wrong storage");
    std::cout << "case_lifecycle_test: stable references and deterministic destruction passed\n";
}
