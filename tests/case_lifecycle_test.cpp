#include "babelsim/case.h"
#include "babelsim/simple.h"
#include "test_util.h"

#include <iostream>

using namespace babelsim;

int main() {
    {
        Case problem("cases/heat", "lifecycle");
        ScalarField& first = problem.scalarField("T");
        for (int i = 0; i < 64; ++i)
            problem.vectorField("scratch" + std::to_string(i), Vec3{});
        require(&first == &problem.scalarField("T"), "Case invalidated a Field reference");
        bool rejected = false;
        try { problem.vectorField("T", Vec3{}); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Case accepted one name for different Field types");
        TensorField& tensor = problem.tensorField("gradient", Tensor3{});
        problem.output(tensor);
        problem.output(first, false);
        problem.faceVectorField("faceVector");
        problem.faceTensorField("faceTensor");
        rejected = false;
        try { problem.output(problem.faceField("faceScalar")); }
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
        SimpleSolver simple(problem);
        problem.scalarField("couplingState", 0.0);
        problem.validate();
        problem.scalarField("afterValidation", 0.0);
        problem.start();
        bool rejected = false;
        try { problem.scalarField("tooLate", 0.0); }
        catch (const std::logic_error&) { rejected = true; }
        require(rejected, "Case start did not close declarations");
    }
    // 前一个 Case 必须先完整释放 RunTime，才能在同一线程创建下一次计算。
    Case second("cases/heat", "lifecycle2");
    require(second.scalarField("T").size() == static_cast<std::size_t>(second.mesh().cellCount()),
            "Case recreated a Field with the wrong storage");
    std::cout << "case_lifecycle_test: stable references and deterministic destruction passed\n";
}
