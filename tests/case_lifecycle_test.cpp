#include "babelsim/case.h"
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
    }
    // 前一个 Case 必须先完整释放 RunTime，才能在同一线程创建下一次计算。
    Case second("cases/heat", "lifecycle2");
    require(second.scalarField("T").size() == static_cast<std::size_t>(second.mesh().cellCount()),
            "Case recreated a Field with the wrong storage");
    std::cout << "case_lifecycle_test: stable references and deterministic destruction passed\n";
}
