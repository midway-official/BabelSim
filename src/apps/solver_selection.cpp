#include "babelsim/case.h"
#include "solver_selection.h"

#include <stdexcept>

namespace babelsim {

// 新求解器只在这里增加声明和一行选择；不需要继承、注册宏或修改通用启动器。
int runHeat(Case& problem);
int runSimple(Case& problem);
int runTransport(Case& problem);

int runSolver(Case& problem) {
    if (problem.solver() == "heat") return runHeat(problem);
    if (problem.solver() == "simple") return runSimple(problem);
    if (problem.solver() == "transport") return runTransport(problem);
    throw std::runtime_error("unknown BabelSim solver: " + problem.solver());
}

}  // babelsim 命名空间
