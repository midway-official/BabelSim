# BabelSim DSL 与运行时用户手册

本文描述当前公开的 C++17 求解器接口和唯一支持的 Case 文件格式。实现以
`include/babelsim/*.h` 为准；本文与头文件冲突时以头文件为准。求解器由普通 C++ 循环驱动，
框架提供 Case、Field、`math`、`equ`、时间历史、并行运行域和结果写出。

## 最小求解器

下面的入口求解 `ddt(T) = div(k grad(T)) + Q`。方程名称、空间配置和线性求解配置都来自
Case 的具名字典：

```cpp
#include "babelsim/application.h"
#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/time.h"

using namespace babelsim;

SolverResult solveHeat(Case& problem) {
    auto& T = problem.scalarField("T");
    const double rhoCp = problem.physics().positive("density") *
                         problem.physics().positive("heatCapacity");
    const double conductivity = problem.physics().nonnegative("conductivity");
    const double source = problem.physics().number("source");
    const auto writeInterval = readWriteInterval(problem);
    const auto control = readEquationControl(problem, "temperature", T, {"diffusion"});
    auto equation = equ::createEquation(T, control);
    auto stepper = time::start(problem);
    auto history = time::history(T);

    while (!stepper.finished()) {
        stepper.advance();
        history.save(T, stepper.dt());
        equation.reset();
        equ::ddt(equation, rhoCp, history);
        equ::laplacian(equation, conductivity, -1.0, "diffusion");
        equ::source(equation, source);
        const auto result = equ::solve(equation);
        if (!result.converged()) return SolverResult::notConverged();
        if (stepper.step() % writeInterval == 0 || stepper.finished())
            write(problem, stepper);
    }
    return SolverResult::completed();
}

const SolverRegistration registration("myHeat", solveHeat);

int main(int argc, char** argv) { return runApplication(argc, argv); }
```

`equ::createEquation(problem, name, field, terms)` 读取 Case 中的完整配置，并绑定方程身份。
如果方程需要在循环外创建一次，也可以先调用 `readEquationControl`，再使用
`equ::createEquation(field, control)`。后者仍然必须传入已解析的 `EquationControl`。

## Case 文件

`case.bs` 只负责选择求解器和文件：

```text
solver simple
mesh mesh/cavity.mesh
fields fields/initial
physics physics/simple.bs
methods numerics/methods.bs
solution numerics/solution.bs
control control.bs
output output.bs
ghostLayers 3
```

推荐文件路径相对于 Case 根目录；重复键、未知键和缺少必需文件会被拒绝。reader
按根目录拼接输入路径，并未统一禁止绝对输入路径。输出路径有单独限制。目录
通常包含：

```text
case.bs
mesh/*.mesh
fields/initial/*.field
physics/*.bs
numerics/methods.bs
numerics/solution.bs
control.bs
output.bs
results/
```

`physics/*.bs` 和 `solution.bs` 返回只读的 `Parameters`。求解器使用
`word`、`number`、`positive`、`nonnegative`、`integer`、`boolean` 等方法读取值；
`problem.validate()` 会拒绝已经配置但没有被读取的键。

## 当前数值格式

`numerics/methods.bs` 只允许 `time`、`equation.*` 和 `operation.*`：

```text
time steady
equation.momentum.interpolation linear
equation.momentum.gradient leastSquares
equation.momentum.convection linearUpwind
equation.momentum.diffusion corrected
equation.momentum.term.diffusion.diffusion limitedCorrected
operation.pressureGradient.interpolation linear
operation.pressureGradient.gradient leastSquares
operation.pressureGradient.convection upwind
operation.pressureGradient.diffusion orthogonal
```

每个具名方程和算子都必须给出完整的 `interpolation`、`gradient`、`convection`、
`diffusion`。项覆盖只写差异项，最终配置由方程基础配置和项覆盖合并得到。系数插值和
系数梯度可用 `coefficientInterpolation`、`coefficientGradient` 单独设置。

`numerics/solution.bs` 同时包含求解器自己的算法键和按方程命名的线性键：

```text
maxIterations 5000
velocityRelaxation 0.5
pressureRelaxation 0.3
continuityTolerance 1e-8
velocityTolerance 1e-6

equation.momentum.kspType bcgs
equation.momentum.pcType bjacobi
equation.momentum.absoluteTolerance 1e-12
equation.momentum.relativeTolerance 1e-8
equation.momentum.maxIterations 1000
```

`cg`、`bcgs`、`gmres`、`fgmres` 是 PETSc KSP 方法；`none`、`icc`、`hypre`、`gamg`、
`bjacobi`、`asm`、`jacobi` 是 PETSc PC 类型。每个实际使用的方程必须提供 `kspType`、
`pcType`、`absoluteTolerance`、`relativeTolerance`、`maxIterations` 五个核心键。`icc` 只用于
串行求解；MPI 可选择 `hypre`、`gamg` 或 `bjacobi`。GAMG 层级选项和 warm-start 键也使用
相同的 `equation.<name>.*` 前缀。
通用后端默认以 8 个矩阵版本为上限复用预条件器；累积 RHS 工作量达到同一上限后，
在下一次矩阵变化时提前重建，矩阵未变化的 RHS-only 调用仍复用原预条件器。
迭代增长也会提前触发重建；失败的旧
预条件器会重建后重试一次。矩阵更新和真实残差标准保持独立。
`preconditionerRebuildInterval 1` 可恢复逐矩阵重建，
`preconditionerReuseIterationGrowth` 设置提前重建的增长倍数（默认 2）。

完整的配置示例、内置方程名称和迁移规则见
[方程、项与算子的数值配置](numerical-configuration.md)。

## Field 与 Case 生命周期

`scalarField`、`vectorField`、`tensorField` 加载 `fields/initial` 中的单元场；
`createScalarField`、`createVectorField`、`createTensorField` 创建程序拥有的场；
`existing*Field` 只取得已经声明的对象。面场使用对应的 `createFace*Field` 接口。
文件加载的单元场默认进入输出列表，程序创建的场需要在 `output.bs` 中列出或通过
`problem.output(field)` 选择。

Case 先完成场和方程声明，再调用 `problem.validate()` 或 `time::start`。时间推进由
`TimeStepper` 显式执行：`advance()` 更新时刻和步长，`history.save()` 保存历史，
`write(problem, stepper)` 写出结果。`Case::write()` 和 `Case::finish()` 只负责输出，
不会替求解器判断收敛。

## `math`：整场显式计算

`math` 操作立即执行并返回独立 Field，不形成延迟表达式。需要空间重构的调用必须带有
完整 `OperatorOptions`，通常直接使用方程项配置：

```cpp
const auto options = equation.options("diffusion");
auto gradient = math::grad(T, options);
auto faceCoefficient = math::interpolate(k, options.coefficientOptions());
auto diffusiveFlux = math::flux(k, T, options);
auto divergence = math::div(diffusiveFlux);
```

也可以通过 `problem.methods().operationOptions("pressureGradient")` 取得独立算子配置。
`math::div(faceFlux)` 是面通量守恒求和，`math::reconstruct(field, gradient)` 使用调用者
提供的梯度；这两个操作没有额外的格式选择。梯度、插值、对流、扩散、拉普拉斯和面通量
重构都不会读取隐式的全局或最近一次方程配置。

## `equ`：过程式方程

`Equation<T>` 绑定一个未知 Field 和一个 `EquationControl`。典型装配顺序是：

```cpp
equation.reset();
equ::ddt(equation, capacity, history);
equ::div(equation, phi, 1.0, "convection");
equ::laplacian(equation, diffusivity, -1.0, "diffusion");
equ::source(equation, forcing);
equ::relax(equation, previous, relaxation);
const auto residual = diagnostics::relativeResidual(equation, unknown);
const auto result = equ::solve(equation);
```

`reset()` 清除上一次装配；`copy()` 复制方程身份和配置但拥有独立的系数。方程项必须在
创建时声明，调用未声明的项会报错。`equ::solve(equation)` 使用绑定的线性配置；需要
临时覆盖时显式传入 `LinearSolverConfig`。

## 运行时和并行边界

普通 Solver 只操作 `Case`、Field、`math`、`equ` 和时间 API。`RunTime` 由 Case 内部创建，
负责当前 Mesh、MPI 全局归约、halo、矩阵装配和线性后端生命周期；Physics 不构造它，
也不接触 Eigen、CSR/LDU 或 rank 本地索引。`runApplication` 负责 MPI 生命周期、Solver
注册和退出码映射。串行和 MPI 使用同一套方程与配置格式。

结果目录、原生 mesh/field 与重启格式见 [网格、场与结果](mesh-and-results.md)。
构建和外部链接见 [构建说明](building.md)。

结果由 Case 写成每个 rank 的拥有单元 CSV 和元数据；`babelsim-post` 再按 global cell ID
合并为 VTK 或 Tecplot。不同 MPI 规模应使用不同的 `-time` 运行名称，避免混合结果目录。

## 内置求解器的稳定方程名

| 求解器 | 方程名 |
| --- | --- |
| `heat` | `temperature` |
| `transport` | `transport` |
| SIMPLE、transient SIMPLE、PISO | `momentum`、`pressureCorrection` |
| k-omega | `kTransport`、`omegaTransport` |
| k-epsilon | `kTransport`、`epsilonTransport` |
| Spalart--Allmaras | `nuTildaTransport` |

方程名同时用于方法配置、线性配置和方程身份检查。新增 Solver 时，应在自己的 Case 中为
每个实际使用的方程写出对应的 `methods.bs` 和 `solution.bs` 块，并在算法启动前调用
`problem.validate()`。

## 验证入口

```bash
make test
make test-full
make test-architecture
make test-external
make test-workflow
make test-rans
make test-scalar-time
make test-simple-parallel
make test-piso-coupling
make test-mpi
```

这些入口分别检查公开 API、配置拒绝、分层约束、外部 Solver、内置工作流、RANS 方程和 MPI
数值一致性。它们不把短时回归结果当成长时间物理精度结论。

### 时间层与面通量的冻结契约

`equ::laplacian` 装配时冻结扩散系数、边界约束及显式非正交面通量，
包括 `limitedCorrected` 当时计算的限幅结果。`equ::faceFlux(A, x)` 仅用当前
`x` 更新隐式正交通量，因此恢复的扩散通量对应已装配的矩阵；方程复制、缩放、
相加也保留这些冻结贡献。

`math::flux(U)` 使用 `U.setBoundaryFlux(phi)` 已保存的边界流向。首次没有
上下文时，用 owner 速度与外法向面积矢量的点积初始化流向；输出缓冲区的旧值
不参与入口/出口分支判定。求值后的通量成为后续边界求值的上下文。

`heat` 和 `transport` 在 BDF2 的第二步起，用
`q* = (1 + r) q_old - r q_older`（`r = dt / dt_previous`）作为空间显式修正
的装配状态，时间导数仍读取真实物理历史。只更新场值，保留物理边界约束。
首步仍用 Euler，缩短的末步使用实际步长比。这保证显式修正的二阶时间一致性；
显式处理仍有稳定性限制，不代表任意时间步下无条件稳定。

## 线性结果、诊断与性能

`equ::solve` 用绑定方程的线性配置，返回结构化状态。最终检查未预条件的全局 L2 真实残差：
`||b-Ax|| <= max(atol, rtol * max(||r0||, ||b||))`。向量三个分量分别求解，全部通过才
收敛；迭代数是实际分量总数。欠松弛前的方程诊断与 KSP 线性残差是不同指标。
具体公式和复用控制见 [数值配置](numerical-configuration.md)。

`monitor::Reporter` 只观测和打印，Solver 判断收敛/时间步接受。应用把完成映射为退出码 0，
迭代耗尽/数值失败映射为 2，配置或应用错误为 1。实际求解失败不能因已有结果文件而判成功。
`-performance` 的循环窗口从首次迭代观测后开始；排除整个首迭代和此前初始化。有效持续
时间、完整工作量及各 rank 关键路径的解释见 [性能测量](performance/README.md)。

新增 Solver 优先参照已测试的 [双场示例](../tests/examples/coupled_scalar.cpp)，为其 Case
提供完整具名方法/线性配置，再执行 `make test-external`。本页是开发入口；全部重载与
参数契约查 [math.h](../include/babelsim/math.h)、[equ.h](../include/babelsim/equ.h)、
[case.h](../include/babelsim/case.h)，不把概述当作全部 API 清单。
