# 方程、项与算子的数值配置

BabelSim 只接受当前的具名配置格式。`numerics/methods.bs` 描述时间格式、方程级空间离散和
可选的方程项/独立算子覆盖；`numerics/solution.bs` 描述算法控制和每个方程的线性求解器。
数值配置按名称绑定到方程，因此同一个 Case 可以为动量、压力修正、温度和湍流输运选择不同
的离散与线性算法。

## `methods.bs`

文件只允许三类键：一个 `time`、具名方程键和具名算子键。

```text
time euler

equation.temperature.interpolation linear
equation.temperature.gradient leastSquares
equation.temperature.convection linearUpwind
equation.temperature.diffusion corrected

equation.temperature.term.diffusion.diffusion limitedCorrected
equation.temperature.term.diffusion.coefficientInterpolation linear
equation.temperature.term.diffusion.coefficientGradient greenGauss

operation.pressureGradient.interpolation linear
operation.pressureGradient.gradient leastSquares
operation.pressureGradient.convection upwind
operation.pressureGradient.diffusion orthogonal
```

每个 `equation.<name>` 和 `operation.<name>` 都必须给出完整的四项空间配置：
`interpolation`、`gradient`、`convection`、`diffusion`。系数重构的
`coefficientInterpolation` 和 `coefficientGradient` 是可选项，省略时分别使用该配置的
插值和梯度方法。方程项覆盖只写需要改变的键，其余键从该方程的基础配置继承。

支持的取值如下：

| 配置项 | 取值 |
| --- | --- |
| `time` | `steady`、`euler`、`bdf2` |
| `interpolation`、`coefficientInterpolation` | `linear`、`corrected` |
| `gradient`、`coefficientGradient` | `greenGauss`、`leastSquares` |
| `convection` | `upwind`、`linearUpwind`、`central` |
| `diffusion` | `orthogonal`、`corrected`、`limitedCorrected` |

这些值按表中的字符串精确匹配；下划线拼写、旧别名和类型级默认值都不再兼容。
重复键、未知键、缺少 `time` 或不完整的方程/算子配置都会在读取或使用配置时失败。框架不再
从全局方法、场名或字段类型推断空间格式，也不再保存一个可供其它方程继承的全局离散对象。

## `solution.bs`

算法级键仍由具体 Solver 读取，例如 SIMPLE 的外迭代、欠松弛和收敛阈值；线性求解参数
必须按方程完整写出：

```text
maxIterations 1200
velocityRelaxation 0.5
pressureRelaxation 0.3
continuityTolerance 1e-7
velocityTolerance 1e-6

equation.momentum.kspType bcgs
equation.momentum.pcType bjacobi
equation.momentum.absoluteTolerance 1e-12
equation.momentum.relativeTolerance 1e-8
equation.momentum.maxIterations 800

equation.pressureCorrection.kspType cg
equation.pressureCorrection.pcType hypre
equation.pressureCorrection.absoluteTolerance 1e-12
equation.pressureCorrection.relativeTolerance 1e-8
equation.pressureCorrection.maxIterations 800
```

`kspType` 支持 `cg`、`bcgs`、`gmres`、`fgmres`；`pcType` 支持 `none`、`icc`、`hypre`、
`gamg`、`bjacobi`、`asm`、`jacobi`。`icc` 是串行预条件器；MPI 算例应选用 `bjacobi`、
`hypre` 或 `gamg`。`warmStart` 和 `amgMaxLevels`、`amgCoarseSize`、`amgSmoothingSteps`
（仅用于 `gamg`）是可选的方程级参数。原 Eigen ILUT 的 drop tolerance/fill factor 在 PETSc
中没有等价配置；`bjacobi` 使用每个 rank 的本地子问题预条件器，数值效果不等同于 ILUT。
每个实际创建的方程都必须有五个核心线性参数，缺少任意一个都会报错。算法控制键和方程线性键
必须使用各自的完整名字；框架不再提供标量/矢量类型默认值或按场覆盖。

## C++ 使用方式

内置 Solver 或外部 Solver 通过稳定的方程名字取得已经解析的配置：

```cpp
auto temperature = equ::createEquation(
    problem, "temperature", T, {"diffusion"});

temperature.reset();
equ::ddt(temperature, rhoCp, history);
equ::laplacian(temperature, conductivity, -1.0, "diffusion");
equ::source(temperature, source);
const auto result = equ::solve(temperature);
```

如果方程在循环中反复创建，可在初始化阶段读取一次：

```cpp
const auto control = readEquationControl(
    problem, "temperature", T, {"diffusion"});
auto temperature = equ::createEquation(T, control);
```

`EquationControl` 同时持有方程级空间配置、项级解析配置和线性求解配置。调用
`equation.options("diffusion")` 可取得项的最终配置；调用 `math` 时应显式传入
`OperatorOptions`：

```cpp
const auto options = temperature.options("diffusion");
auto gradT = math::grad(T, options);
auto faceK = math::interpolate(conductivity, options.coefficientOptions());
```

`math::div(faceFlux)`、`math::reconstruct(field, gradient)` 等本身不需要选择空间格式的
守恒操作仍直接使用输入。需要梯度、插值、对流或扩散重构的操作必须收到完整的显式选项。

## 内置方程名称

| 求解器/模型 | 方程名称 | 未知场 | 常用项 |
| --- | --- | --- | --- |
| `heat` | `temperature` | `T` | `diffusion` |
| `transport` | `transport` | `C` | `convection`、`diffusion` |
| SIMPLE、transient SIMPLE、PISO | `momentum` | `U` | `convection`、`diffusion` |
| SIMPLE、transient SIMPLE、PISO | `pressureCorrection` | `pPrime` | `diffusion` |
| k-omega、k-epsilon | `kTransport` | `k` | `convection`、`diffusion` |
| k-omega | `omegaTransport` | `omega` | `convection`、`diffusion` |
| k-epsilon | `epsilonTransport` | `epsilon` | `convection`、`diffusion` |
| Spalart--Allmaras | `nuTildaTransport` | `nuTilda` | `convection`、`diffusion` |

方程名称是配置和方程身份的一部分。同一个 Case 中，同名方程不能绑定到不同未知场。每个
方程项必须在创建方程时声明；未声明的项名会立即报错。

## 迁移规则

1. 把 `methods.bs` 中的全局或按场离散行改成 `equation.<name>.<option>`；按项差异改成
   `equation.<name>.term.<term>.<option>`。
2. 把 `solution.bs` 中的线性求解行改成完整的 `equation.<name>.kspType`、
   `pcType`、`absoluteTolerance`、`relativeTolerance` 和 `maxIterations`。
3. 用 `equ::createEquation(problem, name, field, terms)` 或
   `readEquationControl(problem, name, field, terms)` 创建方程；不能从字段类型推断配置。
4. 独立 `math` 调用保存并传递对应的 `OperatorOptions`，不要读取隐含的最近方程或全局默认。
5. 为新方程同时加入 `methods.bs` 和 `solution.bs` 的具名块，并在 `problem.validate()` 前
   完成声明和配置检查。

旧的全局/场级离散键、`scalarSolver`/`vectorSolver`、按场线性覆盖、无配置的
`createEquation(field)`、`readLinearControl` 以及无配置求解后端均已删除；它们不会再被解析、
导出或作为兼容入口保留。

## 验证

```bash
make test
make test-architecture
make test-external
make test-workflow
make test-rans
make test-mpi
```

这些测试覆盖配置拒绝、项级继承、公开方程 API、外部 Solver、MPI 后端和内置算例。它们验证
的是配置和数值执行契约，不替代长时间的物理精度验证。
