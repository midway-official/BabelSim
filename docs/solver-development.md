# 求解器开发：从一个函数到一个可运行算例

本文面向掌握 PDE/FVM 和基础 C++ 的 Solver 作者。若要修改离散、通信、存储或线性求解，请读
[架构中的维护者视角](architecture.md)。两类工作不应该从同一堆内部类开始学习。

## 1. 最少需要知道什么

普通 Solver 只包含：

```cpp
#include "babelsim/case.h"
#include "babelsim/solver.h"
```

认识以下概念即可：

| 概念 | 作者负责的内容 |
| --- | --- |
| Case | 从算例获取命名场、物性和算法参数；时间循环 |
| ScalarField / VectorField | 温度、浓度、速度等数学场 |
| Boundary | 按 patch 名称定义数学边界 |
| fvm / fvc | 隐式方程项 / 显式场运算 |
| solve | 求一个离散后的方程并检查线性收敛 |
| diagnostics | 全局数学范数、场变化、守恒误差 |
| 算法循环 | 组织预测、求解、修正和收敛；普通函数和循环即可 |

不需要先学习 RunTime、FieldHistory、EquationSolver、LDU/CSR、MPI、智能指针或工作区。
`solver.h` 的包含关系不再带入 `runtime.h` 或 MPI/Eigen 头文件；自动测试检查这条边界。

必要 C++ 只有普通函数、变量、引用、运算符、条件和循环。`ScalarField& T` 中的
`&` 表示使用 Case 已拥有的场；不要写成 `auto T = ...`，那会复制 Field，
后续修改不再作用于 Case 自动输出的对象。可用 `auto&`，不必使用模板或自行管理所有权。

## 2. 新增方程驱动 Solver：只新增一个源文件

例如开发扩散带体源的浓度方程。在 `src/physics/species/main.cpp` 写：

```cpp
#include "babelsim/case.h"
#include "babelsim/solver.h"

namespace babelsim {
int runSpecies(Case& problem) {
    ScalarField& C = problem.scalarField("C");
    const double D = problem.properties().nonnegative("diffusivity");
    const double S = problem.properties().number("source");

    while (problem.loop()) {
        if (!solve(fvm::ddt(C) ==
                   fvm::laplacian(D, C) + fvm::source(S)).converged()) return 2;
    }
    return 0;
}
}
```

这就是完整 Solver，不需要再写 Case reader、执行入口类、析构函数或输出代码。
也不需要同时提供一个专用头文件和 solveTransientXXX 库函数；原来的 Heat/Transport
重复库式入口已经删除。测试需要内存输入时直接调用同一套通用数学 API。
`Makefile` 自动收集 `src/physics/*/main.cpp`；只在
`src/apps/solver_selection.cpp` 增加声明和一行显式选择：

```cpp
int runSpecies(Case& problem);
// 在 runSolver 内添加：
if (problem.solver() == "species") return runSpecies(problem);
```

因此，最小开发量是**新增一个 Solver 文件，修改一个选择表文件**。不需要新增头文件、
修改通用启动器、编写注册宏、继承基类或修改 MPI/代数层。复杂算法可按数学步骤拆文件，
但不要把 SIMPLE 的多文件规模当作所有 PDE 的必需模板。

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

移除新 Solver 不读取的旧物性/算法参数。`Parameters` 会检查重复键、非有限数字、
缺少项和未消费项，避免拼错参数后悄悄使用默认值。它不是新的物性模型，只是命名配置的读取器。

```bash
make -j4
build/babelsim-solve -case cases/species -time serial
mpirun -np 4 build/babelsim-solve -case cases/species -time mpi4
build/babelsim-post -case cases/species -time mpi4/all -format vtk tecplot
```

打开 `post/mpi4/series.pvd` 即可查看时间序列。Solver 不关心进程数。
不要自行写一个没有初始化 MPI 的 main，再用 mpirun 复制启动；应使用通用启动器。

## 4. 添加对流、变系数和边界

添加对流只需读入速度、形成面通量，再增加一个数学项：

```cpp
VectorField& U = problem.vectorField("U");
ScalarField& phi = problem.faceFlux("phi", U);
// 在时间循环内：
solve(fvm::ddt(C) + fvm::div(phi, C) ==
      fvm::laplacian(D, C) + fvm::source(S));
```

这里只省略了示例中的收敛检查，实际代码仍应检查返回值。若速度随时间更新，在更新后用
`fvc::evaluate(fvc::flux(U), phi)` 更新通量；不要把构造时的通量误当成自动跟随 U 的表达式。

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

两个未知量仍由同一套 `solve/fvm` 求解；不创建耦合执行器、解析器或另一套矩阵。

算法需要保存某个数学状态时可写：

```cpp
ScalarField& previous = problem.scalarField("previous", 0.0);
previous.assign(T);
const double change = diagnostics::relativeChange(T, previous);
```

带初值的重载创建不读取文件、也不自动输出的中间场；Case 拥有它，算法不管理分配或析构。
这是数学上的赋值和变化范数，不是存储访问。初始边界是零梯度；需要别的数学条件时显式设置。

有限次耦合修正使用普通 `for` 循环是合理的，它表达数值算法；禁止的是手写
cell/face 索引循环来实现已有离散。停止判据必须覆盖**全部耦合未知量**，不能只看一个场。
`diagnostics::relativeChange` 已经是全局值，不能再写 rank 本地范数决定是否继续。

时间步内可重复求同一场。历史场在进入下一个物理时间步时推进一次，不随 solve 调用次数推进。
BDF2 首步自动用 Euler；当前 BDF2 只支持均匀步长，非整数步数的终点会被拒绝。

只有算法确实跨多个模块复用时才建立类似 `SimpleSolver` 的算法对象。
SIMPLE 的公开 API 是 loop、五个步骤和 converged；内部状态不属于普通调用者要学习的内容。

## 6. 两类开发者的交接点

| 想做的事 | 普通 Solver 作者 | 框架维护者 |
| --- | --- | --- |
| 新 PDE / 组合已有算子 | 新增一个函数和 Case | 无需修改底层 |
| 新耦合策略 / 停止准则 | 方程顺序、修正、全局 diagnostics | 提供已有场运算和归约 |
| 新物性关系 | 更新物性数学场 | 有通用需求才补场运算核 |
| 新梯度或对流格式 | 在 Case 选方法 | 扩展 Method/离散核并验证 |
| 新线性求解器 / MPI 后端 | 不参与 | 修改代数或执行层 |
| 新数据布局 / 存储优化 | Solver 不变 | 保持所有权、halo、数值与性能契约 |

不要通过创建新 Context/Manager/Wrapper 解决缺少数学运算的问题。先判断它是一个可复用的
数学操作，还是只有本算法需要的步骤；后者留在算法局部。

## 7. 学习与验收顺序

新作者推荐顺序：读一个 Case → 读 heat/main.cpp → 修改一个源项 → 加 div →
读耦合双场例子 → 再读 SIMPLE 的 main/momentum/pressure。无需先通读框架实现。

维护者则从 architecture.md 的所有权与执行契约开始，读 Runtime、离散、代数和并行测试。
`make test-architecture` 独立检查头文件依赖与层次约束，且已纳入下面两个测试目标。
每次改变接口都跑 `make test test-workflow`，涉及执行层还需
`make test-mpi test-mpi-poiseuille`。

当前仍有边界：显式操作有时需要一个命名结果场；没有一般表达式模板、自动单位检查、
通用非线性求解器、restart 或动态插件。不要用“代码短”掩盖这些限制，也不要为未来可能的需求
提前引入复杂框架。
