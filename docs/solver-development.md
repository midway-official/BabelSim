# 求解器开发：从一个函数到一个可运行算例

本文面向掌握 PDE/FVM 和基础 C++ 的 Solver 作者。若要修改离散、通信、存储或线性求解，请读
[架构中的维护者视角](architecture.md)。两类工作不应该从同一堆内部类开始学习。

## 1. 最少需要知道什么

普通 Solver 只包含：

```cpp
#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/solver.h"
```

认识以下概念即可：

| 概念 | 作者负责的内容 |
| --- | --- |
| Case | 从算例获取命名场、物性和算法参数；时间循环 |
| ScalarField / VectorField | 温度、浓度、速度等数学场 |
| Boundary | 按 patch 名称定义数学边界 |
| eqn / math | 隐式方程项 / 显式场运算 |
| solve | 求一个离散后的方程并检查线性收敛 |
| diagnostics | 全局数学范数、场变化、守恒误差 |
| 算法循环 | 组织预测、求解、修正和收敛；普通函数和循环即可 |

不需要先学习 RunTime、FieldHistory、EquationSolver、LDU/CSR、MPI、智能指针或工作区。
`solver.h` 的包含关系不再带入 `runtime.h` 或 MPI/Eigen 头文件；自动测试检查这条边界。

必要 C++ 只有普通函数、变量、引用、运算符、条件和循环。`ScalarField& T` 中的
`&` 表示使用 Case 已拥有的场；不要写成 `auto T = ...`，那会复制 Field，
后续修改不再作用于 Case 自动输出的对象。可用 `auto&`，不必使用模板或自行管理所有权。

## 2. 仓库外新增 Solver：一个源文件，不修改框架

例如开发扩散带体源的浓度方程。在自己的项目里创建 `species.cpp`：

```cpp
#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/solver.h"

namespace babelsim {
int runSpecies(Case& problem) {
    ScalarField& C = problem.scalarField("C");
    const double D = problem.physics().nonnegative("diffusivity");
    const double S = problem.physics().number("source");

    while (problem.loop()) {
        if (!solve(eqn::ddt(C) ==
                   eqn::laplacian(D, C) + eqn::source(S)).converged()) return 2;
    }
    return 0;
}
const SolverRegistration species("species", runSpecies);
}
```

这就是完整 Solver，不需要再写 Case reader、执行入口类、析构函数或输出代码。
也不需要同时提供一个专用头文件和 solveTransientXXX 库函数；原来的 Heat/Transport
重复库式入口已经删除。测试需要内存输入时直接调用同一套通用数学 API。

不要为具体 Solver 向 `include/babelsim/` 增加头文件。公共头表达稳定的框架概念，
具体 Solver 只在自己的实现文件注册运行入口；算法较复杂时，其辅助头也必须留在所属
`src/physics/<solver>/` 目录，并且不能成为其他 Solver 的依赖。
在同一文件加上启动入口：

```cpp
int main(int argc, char* argv[]) {
    return babelsim::runApplication(argc, argv);
}
```

编译自己的文件，再链接已构建的框架库：

```bash
g++ -std=c++17 -I/home/midway/BabelSim/include -c species.cpp -o species.o
mpic++ species.o /home/midway/BabelSim/build/libbabelsim.a -o species
./species -case ./case -time serial
TMPDIR=/tmp mpirun -np 4 ./species -case ./case -time mpi4
```

Solver 编译不需要 MPI/Eigen 头或 `-Isrc`；最终用 MPI 链接器解决框架的并行依赖。
这不是动态插件：每个 Solver 在自己的源文件以一行 `SolverRegistration` 注册自己，
同一个可执行程序可以链接多个 Solver；通用 main 不需要知道它们的函数名。
不修改 Framework、Case reader、Runtime、线性代数或内置应用，也不创建注册宏和基类。

如果希望加入仓库内置命令，将函数和注册行放入 `src/physics/species/main.cpp`，
不包含独立程序的 `int main`。Makefile 自动收集这个文件及同目录的其他 `.cpp`，
不需要修改 `src/apps/babelsim_solve.cpp` 或 Makefile 的 Solver 名单。
这是选择内置分发方式的可选步骤，不是独立二次开发的前提。

注册行放在函数外，名称使用字符串字面量。名称必须唯一，注册顺序不影响选择。
含注册行的目标文件应直接参与链接，如上面的 `species.o`；不要仅把它藏进按需抽取的
静态库，否则链接器可能不载入它。仓库 Makefile 已直接链接各 Solver 的 `main.o`。
原 `SolverEntry` 表和带表参数的 `runApplication` 已移除；迁移为同文件注册后使用两参数入口。

`tests/external/solver.cpp` 给出方程、双场耦合及矢量响应三个实际例子；
`make test-external` 会将文件复制到临时目录，仅使用公开头和预编译库完成构建与 1/2/4 进程运行。

`return 0` 表示计算成功，启动器保证最终输出；`return 2` 表示没有收敛。
不得忽略求解失败后继续推进时间。配置和文件错误由框架给出路径/行号，并停止整个并行作业。

## 3. 创建 Case 并运行

复制相近的案例作为起点，例如 `cases/transport`。修改：

- `case.bs`：选择 `solver species`；
- `fields/initial/C.field`：初值、边界；
- 物性文件：`diffusivity`、`source`；
- `control.bs`：起止时间、步长；
- `methods.bs`、`solution.bs`：离散格式、线性求解器；
- `output.bs`：结果目录、写出间隔。

`solution.bs` 必须同时包含 `scalarSolver` 和 `vectorSolver` 的完整配置，即使新 PDE
只使用其中一类。通用 Case 读取两项，不需要 Solver 自己处理；漏项会直接报错。

线性后端还支持 `gmres amg`、`cg amg`、`bicgstab amg` 和独立的 `amg none`。例如：

```text
scalarSolver gmres amg 1e-14 1e-9 800 gmresRestart=30 amgMaxLevels=12 amgCoarseSize=48 amgSmoothingSteps=2
vectorSolver gmres amg 1e-12 1e-8 800 gmresRestart=30 amgMaxLevels=12 amgCoarseSize=48 amgSmoothingSteps=2
```

这些是 Case 的计算后端选择，不是 Solver 代码中的对象。AMG/GMRES 的层级、Krylov 基和
全局归约由默认后端管理；新增 PDE 只需要组合 `eqn`/`math`/`solve`，不需要修改或调用
线性代数 API。MPI 下独立 AMG 使用局部子域 V-cycle 加全局残差判据；若要开发全局粗网格
或 GPU AMG，应替换 ComputeBackend，而不是在 Physics 中加入通信代码。

移除新 Solver 不读取的旧物性/算法参数。`Parameters` 会检查重复键、非有限数字、
缺少项和未消费项，避免拼错参数后悄悄使用默认值。它不是新的物性模型，只是命名配置的读取器。
读取物理参数使用 `problem.physics()`，与 `case.bs` 中的 `physics` 条目对应；
算法/求解控制使用 `problem.solution()`，与 `solution` 条目对应。

```bash
make -j4
build/babelsim-solve -case cases/species -time serial
mpirun -np 4 build/babelsim-solve -case cases/species -time mpi4
build/babelsim-post -case cases/species -time mpi4/all -format vtk tecplot
```

打开 `post/mpi4/series.pvd` 即可查看时间序列。Solver 不关心进程数。
自定义 main 只调用 runApplication；不要自行处理 MPI 初始化、异常退出或最终输出。

## 4. 添加对流、变系数和边界

添加对流只需读入速度、形成面通量，再增加一个数学项：

```cpp
VectorField& U = problem.vectorField("U");
ScalarField& phi = problem.faceFlux("phi", U);
// 在时间循环内：
solve(eqn::ddt(C) + eqn::div(phi, C) ==
      eqn::laplacian(D, C) + eqn::source(S));
```

这里只省略了示例中的收敛检查，实际代码仍应检查返回值。若速度随时间更新，在更新后用
`math::evaluate(math::flux(U), phi)` 更新通量；不要把构造时的通量误当成自动跟随 U 的表达式。

标量与矢量方程的公开 `solve` 都返回一个 `SolveResult`，统一用 `.converged()` 检查。
矢量方程必须三个分量全部收敛才成功；公开相对残差是最差分量，不能只检查 x 分量。

扩散系数可换为普通命名场：`ScalarField& k = problem.scalarField("k");`。
均匀、分片或温度相关物性属于物理层。当前文件读取支持均匀初值；任意非均匀文件和材料模型库
还未实现，不应宣称已有通用热物性系统。

边界通常来自 Case，也可按数学语义设置：

```cpp
T.boundary("hot") = fixedValue(500.0);
T.boundary("cold") = fixedValue(300.0);
T.boundary("side") = zeroGradient();
```

框架把它们带入离散。不要修改矩阵系数或 halo 边界数组。

## 5. 新增算法驱动 Solver：先用普通函数，不先建类

算法驱动只是在同一个函数里组织多个方程、修正和收敛。实际可编译、串并行运行的最小耦合例子在
`tests/examples/coupled_scalar.cpp`：

\[
\partial_t T=D\nabla^2T+aC,\qquad
\partial_t C=D\nabla^2C+aT.
\]

两个未知量仍由同一套 `solve/eqn` 求解；不创建耦合执行器、解析器或另一套矩阵。

算法需要保存某个数学状态时可写：

```cpp
ScalarField& previous = problem.scalarField("previous", 0.0);
previous.assign(T);
const double change = diagnostics::relativeChange(T, previous);
```

带初值的重载创建不读取文件、也不自动输出的中间场；Case 拥有它，算法不管理分配或析构。
这是数学上的赋值和变化范数，不是存储访问。初始边界是零梯度；需要别的数学条件时显式设置。

所有场在 start/首次 loop 之前创建。validate() 只做校验，算法私有对象的构造不会
抢先结束声明阶段；算法驱动主程序完成声明后调用 problem.start()。
已有命名场在开始计算后仍能查找，但新名称会报错。Case 返回的引用必须保留 `&`，
也不能超出 Case 生命周期。

按坐标定义非均匀源或物性，不需要本地索引：

```cpp
VectorField& force = problem.vectorField("force", Vec3{});
force.evaluate([](Vec3 x) { return Vec3{1.0 + x.y, 0.0, 0.0}; });
problem.output(force);        // 派生 cell 场加入自动输出
problem.output(force, false); // 取消输出
```

位置函数必须只依赖坐标和捕获的物理参数，不依赖进程数或调用次数。
已知场的点值物性可写 k.evaluate(T, [](double temperature) { return 1.0 + 0.01*temperature; })；
跨值类型也可用 energy.evaluate(U, [](Vec3 velocity) { return 0.5*squaredNorm(velocity); })。
同网格同位置检查由框架完成，不向回调暴露存储指针。
张量中间场用 tensorField(name, Tensor3{})；faceField/faceVectorField/faceTensorField
创建面场。当前输出格式只保存 cell 场，选择面场或非 Case 自有场输出会报错。

有限次耦合修正使用普通 `for` 循环是合理的，它表达数值算法；禁止的是手写
cell/face 索引循环来实现已有离散。停止判据必须覆盖**全部耦合未知量**，不能只看一个场。
`diagnostics::relativeChange` 已经是全局值，不能再写 rank 本地范数决定是否继续。

时间步内可重复求同一场。历史场在进入下一个物理时间步时推进一次，不随 solve 调用次数推进。
BDF2 首步自动用 Euler；当前 BDF2 只支持均匀步长，非整数步数的终点会被拒绝。

只有算法确实跨多个步骤共享状态时，才在自己的 `src/physics/<solver>/` 内建立私有算法对象。
这不是公共 Solver API：所有具体求解器都只通过同目录的 `SolverRegistration` 被选择，
`include/babelsim/` 不发布 Heat、Transport、稳态 SIMPLE 或瞬态 SIMPLE 的类。

稳态 SIMPLE 的全部代码在 `src/physics/simple/`；瞬态版本在
`src/physics/transient_simple/`。两个目录分别拥有 `main.cpp`、私有 `algorithm.h/.cpp`、
`momentum.cpp`、`pressure.cpp` 和 `state.h`。瞬态版本复用稳态方程写法和校正顺序，
仅在自己的动量方程加入 `eqn::ddt(rho,U)` 并增加物理时间循环；没有给稳态类增加模式开关，
也没有修改 Runtime、MPI 或离散后端。

普通 Case 用户只写 `solver simple` 或 `solver transientSimple`。开发新的 Solver 时也不要
包含任何其他求解器目录的 algorithm.h/state.h；应组合公共 Case/Field/eqn/math/solve。
`simple_common.h` 只是两个内置 SIMPLE 的 Physics 私有数据类型，不是二次开发入口。

算法需要了解非正交方法时，使用只读的 `numericalMethods()`；需要输出诊断时，
使用 `diagnostics::report("说明本次迭代的数学诊断")`，不要自己查询主进程或包含 runtime.h。
这不要求普通方程驱动 Solver 创建更多对象，默认 solve/math 的方法仍来自 Case。

`make test-external` 会把稳态 SIMPLE 私有模块作为维护对象复制到仓库外重新编译，
同时验证真正的外部 Solver 只使用公共 Framework 头。新增自己的算法可以保留少量本模块
私有文件，但不能包含框架 internal/ 或其他求解器的私有头。

## 方程控制和动量响应

普通标量/矢量源统一使用 `eqn::source(F)` 或 `eqn::source(a,F)`。
已知矢量场源不再需要转换成手写逐分量循环。

```cpp
SolveResult result = solveWithResponse(
    eqn::div(phi, U) == -math::grad(p) + eqn::laplacian(nu, U),
    rAU, relaxed(0.7));
```

rAU 是对角体积响应，不是矩阵引用。必须检查 result.converged()；
它的缩放约定与框架 SIMPLE 一致，详见 [eqn/math](eqn-math.md)。

对于相容的全 Neumann 标量 Poisson 方程，使用 `solve(equation, referenceValue(0.0))`
指定零空间规范。这不是任意方程的通用点约束，不能用来代替固定值边界。
通常 Solver 只用默认 solve；这些控制用于有明确数值需求的耦合算法。

## 6. 两类开发者的交接点

| 想做的事 | 普通 Solver 作者 | 框架维护者 |
| --- | --- | --- |
| 新 PDE / 组合已有算子 | 新增一个函数和 Case | 无需修改底层 |
| 新耦合策略 / 停止准则 | 方程顺序、修正、全局 diagnostics | 提供已有场运算和归约 |
| 新物性关系 | 更新物性数学场 | 有通用需求才补场运算核 |
| 新梯度或对流格式 | 在 Case 选方法 | 扩展 Method/离散核并验证 |
| 新线性求解器 / MPI 后端 | 不参与 | 实现/替换 ComputeBackend；Solver 与 FVM 表达式不变 |
| 新数据布局 / 存储优化 | Solver 不变 | 保持所有权、halo、数值与性能契约 |

不要通过创建新 Context/Manager/Wrapper 解决缺少数学运算的问题。先判断它是一个可复用的
数学操作，还是只有本算法需要的步骤；后者留在算法局部。

## 7. 学习与验收顺序

新作者推荐顺序：读一个 Case → 读 heat/main.cpp → 修改一个源项 → 加 div →
读耦合双场例子 → 再读 SIMPLE 的 main/momentum/pressure。无需先通读框架实现。

维护者则从 architecture.md 的所有权与执行契约开始，读 Runtime、离散、代数和并行测试。
`make test-architecture` 独立检查头文件依赖与层次约束，且已纳入下面两个测试目标。
每次改变接口都跑 `make test test-workflow test-external`，涉及执行层还需
`make test-mpi test-mpi-poiseuille`。

当前仍有边界：显式操作有时需要一个命名结果场；没有一般表达式模板、自动单位检查、
通用非线性求解器、restart 或动态插件。不要用“代码短”掩盖这些限制，也不要为未来可能的需求
提前引入复杂框架。
