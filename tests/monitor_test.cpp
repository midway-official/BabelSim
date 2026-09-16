#include "babelsim/monitor.h"
#include "babelsim/runtime.h"
#include "test_util.h"
#include <iostream>
#include <sstream>

int main() {
    using namespace babelsim;
    auto mesh = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
    auto runtime = RunTime::forMesh(mesh);
    const monitor::Reporter reporter("generic", 10);
    std::ostringstream output;
    auto* previous = std::cout.rdbuf(output.rdbuf());
    reporter.record({{"energy", 1.5}, {"samples", 2}, {"ready", true}, {"phase", "initial"}});
    for (int i = 1; i <= 12; ++i)
        reporter.iteration(i, 12, i == 3, {{"error", 0.25}, {"accepted", false}});
    std::cout.rdbuf(previous);
    require(output.str() ==
        "generic energy=1.5 samples=2 ready=true phase=initial\n"
        "generic 1 error=0.25 accepted=false\n"
        "generic 3 error=0.25 accepted=false\n"
        "generic 10 error=0.25 accepted=false\n"
        "generic 12 error=0.25 accepted=false\n",
        "generic reporter metrics/cadence or caller-requested event changed");
    bool rejected = false;
    try { monitor::Reporter invalid("generic", 0); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "zero report interval accepted");
    std::cout << "monitor_test: generic observations, cadence and urgent event passed\n";
}
