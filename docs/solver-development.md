# BabelSim 新物理求解器开发

本指南给出从一个普通 C++ 函数开始的完整流程。先读 [DSL 语义参考](procedural-dsl.md)、
[案例配置](case-structure.md) 和 [架构边界](architecture.md)。

## 1. 单文件入口

~~~cpp
#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/solver.h"
#include "babelsim/monitor.h"

namespace babelsim {
SolverResult runMySolver(Case& problem) {
    // load -> assemble -> solve -> observe
}
const SolverRegistration mySolver("mySolver", runMySolver);
}
~~~

不需要 Solver 类、继承、Manager 或 Factory。Case 拥有命名场；算法状态留在函数中；math
临时量用值语义。application 负责 MPI 生命周期、Case 构造和错误映射，Physics 不接触 MPI。

## 2. 配置和声明阶段

case.bs 只声明 solver、mesh、fields、physics、methods、solution、control、output。
physics 放物性、模型和模型常数；methods 放 interpolation、gradient、convection、diffusion、
time；solution 放 scalarSolver/vectorSolver、松弛、迭代上限和容差；control 放
startTime、endTime、deltaT；output 放 directory、timeName、writeInterval、writeFields、
excludeFields。

~~~cpp
auto& T = problem.scalarField("T");
auto& phi = problem.createFaceField("phi");
const auto& physics = problem.physics();
const double rho = physics.positive("density");
const double D = physics.nonnegative("diffusivity");
const int limit = problem.solution().integer("maxIterations", 1000, 1, 1000000);
const auto linear = readLinearControl(problem, T);
problem.validate();
~~~

加载场和创建场必须区分，复用场使用 existing*Field。validate 只校验配置，不推进时间、
不写文件；完成它之后再开始算法。

## 3. 让代码读成 PDE

~~~cpp
auto time = time::start(problem);
auto history = time::history(T);
auto equation = equ::createEquation(T);
while (time.value() < time.end()) {
    time.advance();
    history.save(T, time.dt());

    equation.reset();
    equ::ddt(equation, rho, history);
    equ::div(equation, phi);
    equ::laplacian(equation, D, -1);
    equ::source(equation, source);

    const auto solved = equ::solve(equation, T, linear);
}
~~~

Equation 明确绑定未知量；reset、ddt、div、laplacian、source 的顺序就是离散算法。math
只表示显式场数学，equ 只表示加入方程的离散项。不要让运行时猜测 PDE，也不要访问 CSR、
Eigen 或原始 Field 数组。

## 4. 循环、收敛和监视

稳态使用一个主要外循环；瞬态使用显式时间循环和自己的内层迭代。允许 pressure 方程的
non-orthogonal 修正子循环。残差应在组装后、relax/solve 前检查，避免为诊断再组装：

~~~cpp
const double r = diagnostics::relativeResidual(equation, T);
const double d = diagnostics::relativeChange(T, previous);
const bool converged = solved.converged() && std::isfinite(r)
    && r <= tolerance && d <= changeTolerance;
report.iteration(iter + 1, limit, converged, {{"r", r}, {"change", d}});
~~~

Reporter 是通用报告器，不知道 SIMPLE、RANS 或收敛标准；Solver 自己决定继续并返回
SolverResult。不要在 Physics 直接 return 2。

## 5. SIMPLE 和 RANS

稳态 SIMPLE 与瞬态 SIMPLE 是两个独立 main.cpp。每个 main 显式写动量组装、组装后残差、
松弛、速度求解、Rhie–Chow、压力修正、速度/通量修正、湍流输运和收敛。Rhie–Chow 只在
对应 main 中用 math::grad、math::interpolate、math::flux、math::add/subtract 组合，
不能放入通用算子。

RANS 的 SA、k–omega、k–epsilon 各自拥有输运方程和闭合；SIMPLE 只使用
rans::load、effectiveViscosity、deviatoricStressRemainder 和 solveTransport。模型常数放
physics，模型数值控制放 solution。

## 6. 输出、注册和验证

文件场默认输出，程序场在 output.bs 的 writeFields 或 problem.output(field) 中选择，Face
场当前不直接写出。每个 Solver 文件末尾放唯一 SolverRegistration；Makefile 加入自己的
源文件，不修改集中式 Solver 表。

最低验证命令：

~~~bash
make -j4 all
make -j4 test
python3 tests/solver_workflow_test.py
python3 tests/external_solver_test.py
~~~

涉及 MPI、非正交或 RANS 时再运行对应目标；同时覆盖配置错误、load/create/existing 顺序、
解析或制造解、1/2/4 rank 和结构化失败返回。
