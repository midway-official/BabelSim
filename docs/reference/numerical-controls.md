# 数值配置参考

[配置参考](README.md) · [数值使用说明](../user-guide/numerics.md)

## methods 语法与优先顺序

```text
time euler
equation.transport.interpolation linear
equation.transport.gradient leastSquares
equation.transport.convection upwind
equation.transport.diffusion corrected
equation.transport.coefficientInterpolation corrected
equation.transport.coefficientGradient leastSquares
equation.transport.term.convection.convection linearUpwind
```

基础方程配置必须完整。term 配置只填写要覆盖的项；装配时优先顺序为基础方程 → 命名 term 覆盖 → 显式调用参数。调用参数由 C++ 求解器提供，用户不能通过未声明的名字改变任意运算。

| option | 值 | 必需/默认 |
| --- | --- | --- |
| `interpolation` | linear / corrected | 基础配置必需 |
| `gradient` | greenGauss / leastSquares | 基础配置必需 |
| `convection` | upwind / linearUpwind / central | 基础配置必需 |
| `diffusion` | orthogonal / corrected / limitedCorrected | 基础配置必需 |
| `coefficientInterpolation` | linear / corrected | 缺省继承基础 interpolation |
| `coefficientGradient` | greenGauss / leastSquares | 缺省继承基础 gradient |

`time` 必需，可选 `steady/euler/bdf2`，只能定义一次。基础与覆盖键均严格校验，未知方法名、重复键、废弃全局写法与未使用配置会报错。

独立已知场运算使用完整的 `operation.<name>.<option>` 配置，C++ 通过 `problem.methods().operationOptions(name)` 读取。例如由代码明确读取 `velocityGradient` 后，可配置该 operation 的四个基础选项。未使用的 operation 不应放入字典。

## 正式方程名称

| 求解器/模型 | 方程名 | 主要未知量 | 支持的命名项 |
| --- | --- | --- | --- |
| heat | temperature | T | diffusion |
| transport | transport | C | convection、diffusion |
| flow | momentum | U | convection、diffusion |
| SIMPLE 系列 | pressureCorrection | 压力校正量 | diffusion |
| PISO | pressureCorrection | 压力求解过程 | diffusion |
| kOmega | kTransport、omegaTransport | k、omega | convection、diffusion |
| kEpsilon | kTransport、epsilonTransport | k、epsilon | convection、diffusion |
| SA | nuTildaTransport | nuTilda | convection、diffusion |

基础键使用方程名而非场名。自定义求解器的名称由 `equ::createEquation` 的参数定义。

## solution：逐方程线性配置

为每个实际创建的方程设置 `equation.<name>.<key>`：

| key | 必需/默认 | 范围与含义 |
| --- | --- | --- |
| kspType | 必需 | cg / bcgs / gmres / fgmres |
| pcType | 必需 | none / icc / hypre / gamg / bjacobi / asm / jacobi |
| absoluteTolerance | 必需 | 有限正数 |
| relativeTolerance | 必需 | 有限正数 |
| maxIterations | 必需 | 正整数 |
| warmStart | false | 布尔，使用当前未知量作为初始猜测 |
| preconditionerRebuildInterval | 8 | 正整数，PC 周期重建控制 |
| preconditionerReuseIterationGrowth | 2 | 有限数且至少 1，迭代数增长触发重建 |
| amgMaxLevels | 12 | 正整数，仅 GAMG |
| amgCoarseSize | 48 | 正整数，仅 GAMG |
| amgSmoothingSteps | 2 | 正整数，仅 GAMG |

GAMG 三项不用于修改 Hypre 的层级。`icc` 在 MPI 多进程下不支持。历史 ILUT、drop/fill 等键不是现行后端配置。

```text
equation.temperature.kspType bcgs
equation.temperature.pcType hypre
equation.temperature.absoluteTolerance 1e-14
equation.temperature.relativeTolerance 1e-10
equation.temperature.maxIterations 1000
equation.temperature.warmStart false
equation.temperature.preconditionerRebuildInterval 8
equation.temperature.preconditionerReuseIterationGrowth 2
```

PC 复用只影响预条件器生命周期，矩阵贡献仍按装配更新。复用 PC 的求解失败时，后端可使用当前矩阵重建并从原始猜测重试一次；这不是继续迭代旧矩阵。

线性接受阈值按真实残差 `||b-Ax|| <= max(atol, rtol*max(||r0||, ||b||))` 判断，零缩放时报告使用安全缩放。向量方程要求全部分量成功；报告迭代数与残差须按该接口的合并规则理解。

## solution：SIMPLE / transientSimple

这些键是 solution 的顶层项，不带 `equation.` 前缀：

| key | 默认 | 约束/用途 |
| --- | --- | --- |
| maxIterations | 1000 | 正整数，SIMPLE 外层或每个物理步的内层上限 |
| nonOrthogonalCorrections | 1 | 整数 0..20，修正扩散的追加压力校正 |
| velocityRelaxation | 0.7 | `(0,1]`，速度预测方程松弛 |
| pressureRelaxation | 0.3 | `(0,1]`，压力更新松弛 |
| continuityTolerance | 1e-8 | 正数，质量指标阈值 |
| velocityTolerance | 1e-7 | 正数，速度变化阈值 |
| momentumTolerance | 1e-6 | 正数，动量方程指标阈值 |
| pressureCorrectionTolerance | 1e-6 | 正数，压力校正指标阈值 |

正交扩散实际只需单次对应压力求解；非正交配置根据算法执行基础求解加追加校正。稳态 SIMPLE 要求 methods 为 steady；transientSimple 要求 euler/bdf2。达到上限却未满足耦合条件返回未收敛，不应通过最后一条日志伪装成功。

## solution：PISO

| key | 默认 | 约束/用途 |
| --- | --- | --- |
| maxIterations | 1 | 正整数，每步耦合流程次数 |
| nCorrectors | 2 | 整数 1..100，最低压力校正次数 |
| maxCorrectors | nCorrectors | 整数 nCorrectors..100 |
| couplingTolerance | 最大 double | 正数，默认近似不限制相应耦合误差 |
| nonOrthogonalCorrections | 1 | 整数 0..20 |
| initialNonOrthogonalCorrections | nonOrthogonalCorrections | 整数 0..20，初始投影使用 |
| velocityRelaxation | 1 | `(0,1]`，1 关闭动量预测松弛 |
| continuityTolerance | 1e-8 | 正数 |
| velocityTolerance | 1e-7 | 正数 |
| momentumTolerance | 1e-6 | 正数 |
| pressureCorrectionTolerance | 1e-6 | 正数 |
| initializePotentialFlow | 0 | 整数 0 或 1 |
| turbulenceCoupling | iterated | iterated / segregated |
| turbulenceMaxIterations | 1000 | 正整数 |
| turbulenceCorrectors | 按模式选择 | 正整数，不能超过 turbulenceMaxIterations |

启用模型、`maxIterations=1` 且 `iterated` 时，turbulenceCorrectors 默认取 turbulenceMaxIterations；其他情况下默认 1。segregated 模式要求 maxIterations=1、turbulenceCorrectors=1，启用模型时还要求 turbulenceRelaxation=1。

接受步至少需要线性成功、质量与配置的耦合条件通过。单次瞬态流程不要求速度已经稳态，但 iterated 湍流仍需收敛；多次流程还使用相应外层变化/残差条件。Courant 数只记录，不自动控制步长。该求解器不读取 pressureRelaxation，将 SIMPLE 参数块直接复制过来会触发未使用项错误。

## 自定义教程

coupledScalar 的 solution 顶层项为 `couplingIterations`（默认 100，1..10000）和 `couplingTolerance`（默认 1e-10，正数），以及完整的 `equation.temperature`、`equation.concentration`。它要求 time euler。不要将这些教程键加入内置热/流动求解器。
