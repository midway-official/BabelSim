# BabelSim 求解器开发指南

本文面向熟悉 PDE/FVM、掌握基础 C++、但不需要了解 MPI 或稀疏矩阵实现的开发者。

## 先选对层次

新增 Solver 时只应依赖以下概念：

```text
Mesh、ScalarField/VectorField、边界条件、Material/物性、fvm、fvc、RunTime、solve
```

不要在 Physics 源码中使用：

```text
MPI_*、ParallelContext、HaloExchange、rank、ghost、buffer
Equation/LDU、SparseAssembly、Eigen、PreparedLinearSolver、DistributedLinearSolver
Field::data()、mutableData()、values()、手写单元循环来组装 PDE
```

这些并非“不能存在”，而是 Framework Internal API。Runtime 已经处理它们。

## 第一个标量 Solver：热/扩散

1. 在 Case 中放置初值和边界条件；
2. 定义物性；
3. 创建 `RunTime`；
4. 用 `fvm` 写方程；
5. 在时间循环内调用 `solve`。

```cpp
#include "babelsim/thermal.h"
#include "babelsim/fvm.h"

void solveMyScalar(
    RunTime& run_time,
    ScalarField& T,
    double rho_cp,
    double k,
    double Q)
{
    while (run_time.loop()) {
        const SolveResult result = solve(
            run_time,
            fvm::ddt(rho_cp, T)
                == fvm::laplacian(k, T) + fvm::source(Q));
        if (!result.converged()) {
            throw std::runtime_error("temperature linear solve did not converge");
        }
    }
}
```

这段源码已经是完整的物理表达；不会因 `mpirun -np 4` 而改变。若系数随空间变化，把
`rho_cp` 或 `k` 换成预先建立的 cell scalar Field 即可。

## 边界条件

Field 创建后按 patch 名称指定数学边界：

```cpp
T.boundary("hot") = fixedValue(500.0);
T.boundary("cold") = fixedValue(300.0);
T.boundary("wall") = zeroGradient();
U.boundary("wall") = fixedValue(Vec3{});
U.boundary("symmetry") = symmetry();
```

边界条件会在 fvm 装配和 fvc 面重构时自动参与。不应通过修改矩阵系数或面数据实现边界。

## 使用显式量

显式量由 `fvc` 描述，写入开发者提供的工作 Field：

```cpp
VectorField grad_p(mesh, FieldLocation::Cell, "grad(p)");
ScalarField phi(mesh, FieldLocation::Face, "phi");
ScalarField div_phi(mesh, FieldLocation::Cell, "div(phi)");

run_time.evaluate(fvc::grad(p), grad_p);
run_time.evaluate(fvc::flux(U), phi);
run_time.evaluate(fvc::div(phi), div_phi);
```

`evaluate` 会验证 Mesh/位置、同步必要的输入和结果 ghost，并按 Case 指定的方法处理非正交
和偏斜；它不会生成隐式方程。

## 对流扩散方程模板

\[
\frac{\partial C}{\partial t}+\nabla\cdot(\phi C)
=\nabla\cdot(D\nabla C)+S
\]

对应：

```cpp
solve(
    run_time,
    fvm::ddt(C)
      + fvm::div(phi, C)
        == fvm::laplacian(D, C) + fvm::source(S));
```

`phi` 是 face scalar Field，`C` 是 cell scalar Field，`D` 可以是常数、cell Field 或
face Field。对流格式、梯度和扩散的非正交方式由案例 `numerics` 配置，不应硬编码到 Solver。

## 编写矢量方程

矢量未知量使用相同表达式规则。压力梯度是显式 `fvc::grad(p)`，因此动量方程可读为：

```cpp
solve(
    run_time,
    fvm::div(phi, U)
        == -fvc::grad(p) + fvm::laplacian(mu, U));
```

若算法需要从动量对角得到 mobility（例如 SIMPLE 的 `rAU`），传入 `VectorEquationControl`；
这仍然不暴露矩阵或行号：

```cpp
VectorEquationControl control;
control.relaxation = 0.7;
control.mobility = &rAU;
solve(run_time, momentum_equation, control);
```

除压力速度耦合等确有领域意义的算法外，普通 PDE 不应创建专用 `TemperatureMatrix`、
`EquationManager` 或 `SolverManager`。

## 新 Solver 的文件与 Case

建议将一个 Solver 放在 `src/physics/<领域>/`，并保持一个短的核心源文件。Case 负责
Mesh、初值、物性、方法、线性控制和输出；Solver 负责方程与算法。完整组织见
[case-structure.md](case-structure.md)。

新增过程：

1. 新建物理/算法源文件，只包含 Public Solver API；
2. 为物性和 Case 字典增加小型值对象/读取器；
3. 在 `babelsim-solve` 中按 `case.bs` 的 `solver` 选择它；
4. 增加至少一个解析解、守恒、对称或基准回归；
5. 为串行和 MPI 运行比较设置合理的绝对/相对容差。

只有当某段算法跨多个 Solver 可复用，并有明确数学/数值语义时，才应下沉为通用 `fvm/fvc`
或新的 Method；只有当它是稳定的领域算法（如 `MomentumInterpolation`）时才应作为专用组件。

## 不要误用 Runtime

`RunTime` 是执行语义而不是通用 Manager。它只负责时间、方程离散、solve、Field 同步和
全局诊断。Solver 不保存 MPI 对象，也不创建线性求解器。对稳态 SIMPLE 可以由应用控制外
迭代；对瞬态标量输运使用 `while (run_time.loop())`。

如果新物理确实要逐单元计算经验源项，应先考虑把它定义为可复用的 Field 计算；仅在无法用
`fvc` 表达时，在物理模块中实现一个短且有明确名称的计算核，并避免接触 ghost、halo 和 MPI。
