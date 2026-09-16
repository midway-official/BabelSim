#include "babelsim/case.h"
#include "internal/field_access.h"
#include "support/simple_reference.h"
#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    {
        Case problem("cases/heat", "lifecycle");
        const Parameters& physics = problem.physics();
        require(&physics == &problem.physics(), "Case recreated its physics dictionary");
        require(near(physics.positive("density"), 1.0) &&
                near(physics.positive("heatCapacity"), 1.0) &&
                near(physics.nonnegative("conductivity"), 0.1) &&
                near(physics.number("source"), 0.0),
                "Case physics values do not match the selected dictionary");
        problem.validate();
        ScalarField& first = problem.scalarField("T");
        for (int i = 0; i < 64; ++i)
            problem.createVectorField("scratch" + std::to_string(i), Vec3{});
        require(&first == &problem.scalarField("T"), "Case invalidated a Field reference");
        bool rejected = false;
        try { problem.createVectorField("T", Vec3{}); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Case accepted one name for different Field types");
        TensorField& tensor = problem.createTensorField("gradient", Tensor3{});
        ScalarField& initialized = problem.createScalarField("initialized", 1.0);
        bool initialization_rejected = false;
        try { problem.createScalarField("initialized", 2.0); }
        catch (const std::logic_error&) { initialization_rejected = true; }
        require(initialization_rejected,
                "a repeated create call silently ignored its initialization value");
        bool load_create_rejected = false;
        try { problem.scalarField("initialized"); }
        catch (const std::logic_error&) { load_create_rejected = true; }
        require(load_create_rejected,
                "a created field was silently treated as a file-loaded field");
        require(near(detail::fieldData(initialized)[0], 1.0),
                "created field initialization was not applied");
        problem.output(tensor);
        problem.output(first, false);
        problem.createFaceVectorField("faceVector");
        problem.createFaceTensorField("faceTensor");
        rejected = false;
        try { problem.output(problem.createFaceScalarField("faceScalar")); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "unsupported face output was silently accepted");
        ScalarField copy(first);
        rejected = false;
        try { problem.output(copy); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "non-owned Field was registered for Case output");
    }
    {
        Case problem("cases/cavity", "lifecycle-simple");
        problem.vectorField("U"); problem.scalarField("p"); problem.createFaceScalarField("phi");
        problem.physics().positive("density");
        problem.physics().positive("dynamicViscosity");
        readSimpleControl(problem.solution());
        problem.createScalarField("couplingState", 0.0);
        problem.validate();
        problem.createScalarField("afterValidation", 0.0);
        problem.start();
        bool rejected = false;
        try { problem.createScalarField("tooLate", 0.0); }
        catch (const std::logic_error&) { rejected = true; }
        require(rejected, "Case start did not close declarations");
    }
    // 前一个 Case 必须先完整释放 RunTime，才能在同一线程创建下一次计算。
    Case second("cases/heat", "lifecycle2");
    require(second.scalarField("T").size() == static_cast<std::size_t>(second.mesh().cellCount()),
            "Case recreated a Field with the wrong storage");
    std::cout << "case_lifecycle_test: stable references and deterministic destruction passed\n";
}
