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

每个方程还支持通用的预条件器复用控制，不依赖 SIMPLE/PISO 或算例名称：

```text
equation.momentum.preconditionerRebuildInterval 8
equation.momentum.preconditionerReuseIterationGrowth 2
```

默认在至多 8 个不同矩阵版本间复用预条件器；当前矩阵 A 仍在每次系数改变时更新。
相同矩阵只更新 RHS，不触发周期重建。每次有实际 Krylov 迭代的逻辑 RHS 求解还累积
工作预算；达到同一间隔上限后，在下一次矩阵变化时提前重建。这使一版矩阵承载很多
RHS 时及时重新准备预条件器，避免旧层级的少量迭代增长在整批求解中累积。
计数在各 rank 上一致，不依赖计时或额外归约。某个 RHS 的迭代数超过新预条件器参考次数的
2 倍时，下次求解提前重建；矢量方程按分量保存参考次数。配置改变也重新准备预条件器。
旧预条件器的求解未满足当前矩阵的真实残差要求时，重建并从原初猜重试一次；第二次失败
仍按原收敛/失败合同返回。重试使用相同容差，并计入实际迭代数和耗时。
`preconditionerRebuildInterval 1` 恢复每次矩阵数值更新都准备预条件器的策略。
间隔必须为正整数，增长倍数必须有限且不小于 1。`pcType none` 不使用这项复用策略。
复用会改变 Krylov 路径和迭代次数；用同一精度要求的总计算时间和结果验证评价收益。

`pcType hypre` 的 CPU 后端统一使用 HMIS 粗化、ext+i 插值、每行最多 4 个插值项，
CG、GMRES、FGMRES 使用 1 层 aggressive coarsening；BiCGStab 使用常规粗化和 CF 松弛顺序，
避免在已验证的非奇异扩散算子上发生 shadow-residual breakdown。
其他 Krylov 方法在串行保留 CF 松弛顺序，分布式 MPI 使用自然顺序；
两者均保留对称 SOR/Jacobi 平滑、一次 V-cycle 和零预条件器内部容差。该设置在后端
按 Krylov 类型和执行方式选择，适用于所有选择 Hypre 的方程，不读取算例名称或几何类别。
它遵循 [Hypre 的 CPU BoomerAMG 配置建议](https://hypre.readthedocs.io/en/stable/solvers-boomeramg.html)，
并以本机持续运行和真实残差验证选择松弛顺序。
配置改变 Krylov 路径，但不改变方程线性容差、失败判断或最终的 `b-A*x` 检查。
后端通过临时私有 PETSc options 数据库解析高级选项，不修改全局 options；
`amgMaxLevels`、`amgCoarseSize`、`amgSmoothingSteps` 继续仅用于 GAMG。

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

反复求解的方程在初始化时创建一次，循环中 reset/装配/solve，以保留后端生命周期：

```cpp
const auto control = readEquationControl(
    problem, "temperature", T, {"diffusion"});
auto temperature = equ::createEquation(T, control);
```

`EquationControl` 同时持有方程级空间配置、项级解析配置和线性求解配置。调用
`equation.options("diffusion")` 可取得项的最终配置；调用 `math` 时应显式传入
`OperatorOptions`（下例 k 是已加载/创建的 ScalarField 系数；常数系数直接传给 laplacian）：

```cpp
const auto options = temperature.options("diffusion");
auto gradT = math::grad(T, options);
auto faceK = math::interpolate(k, options.coefficientOptions());
```

`math::div(faceFlux)`、`math::reconstruct(field, gradient)` 等本身不需要选择空间格式的
守恒操作仍直接使用输入。需要梯度、插值、对流或扩散重构的操作必须收到完整的显式选项。

## 内置方程名称

| 求解器/模型 | 方程名称 | 未知场 | 常用项 |
| --- | --- | --- | --- |
| `heat` | `temperature` | `T` | `diffusion` |
| `transport` | `transport` | `C` | `convection`、`diffusion` |
| SIMPLE、transient SIMPLE、PISO | `momentum` | `U` | `convection`、`diffusion` |
| SIMPLE、transient SIMPLE | `pressureCorrection` | `pPrime` | `diffusion` |
| PISO | `pressureCorrection` | `p`（总压力），潜势初始化另用 `pPrime` | `diffusion` |
| k-omega、k-epsilon | `kTransport` | `k` | `convection`、`diffusion` |
| k-omega | `omegaTransport` | `omega` | `convection`、`diffusion` |
| k-epsilon | `epsilonTransport` | `epsilon` | `convection`、`diffusion` |
| Spalart--Allmaras | `nuTildaTransport` | `nuTilda` | `convection`、`diffusion` |

方程名称是配置和方程身份的一部分。Case 读取配置时，同名方程不能绑定到不同未知场；PISO 先读取绑定 p 的
控制，再显式传给私有潜势方程。后端身份还含标量/向量和未知场名，不能仅按配置名推断系统。每个
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

## 容差、初猜与残差口径

五个核心线性键必须显式提供。`warmStart` 可选、默认 false；`preconditionerRebuildInterval`
默认 8，`preconditionerReuseIterationGrowth` 默认 2。GAMG 可选键默认分别为
`amgMaxLevels=12`、`amgCoarseSize=48`、`amgSmoothingSteps=2`；这些键不配置 Hypre。
容差需有限且为正，迭代上限需为正，详细检查见 [solver_control.h](../include/babelsim/solver_control.h)。

最终收敛按当前矩阵的未预条件全局 L2 残差判定：

```text
r0 = b - A*x0
scale = max(||r0||_2, ||b||_2)
||b-A*x||_2 <= max(absoluteTolerance, relativeTolerance * scale)
```

scale 为零时相对残差显示用 1，绝对判据仍生效。重新准备 PC 后的重试保持原初猜及同一
当前矩阵和目标，不因第一次尝试改变尺度。向量各分量均需收敛，实际迭代数求和。
外层 `diagnostics::relativeResidual` 观察原方程，使用 `||Ax||+||b||` 尺度，与线性
KSP 判据不同；不能把两类指标互换或用一个放宽另一个。

使用 CG 需确认矩阵和 PC 的适用性质，设置名本身不能证明对称正定。更换 Krylov/PC 或
复用策略可能改变迭代路径；用同精度结果、真实残差及总计算时间评估。配置收益不等于
旧 Eigen 与当前 PETSc 的直接对照，历史报告只表示当时实验。
