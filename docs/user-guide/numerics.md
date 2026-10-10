# 07 数值设置与收敛

[用户手册](README.md) · 上一章：[选择物理模型](models.md) · 下一章：[运行与并行计算](running-and-parallel.md)

## 三类控制分别放置

`control.bs` 管理物理时间；`methods.bs` 管理离散格式；`solution.bs` 管理线性求解和算法迭代。增加线性迭代次数不会自动减小时间步，收紧线性容差也不会自动增加 PISO 校正次数。

## 完整的方程离散设置

```text
time euler
equation.temperature.interpolation linear
equation.temperature.gradient greenGauss
equation.temperature.convection upwind
equation.temperature.diffusion orthogonal
```

每个被求解器声明的方程都需要四个基础空间选项，即使该方程本次没有对流项。逐项覆盖采用：

```text
equation.transport.term.convection.convection linearUpwind
equation.transport.term.diffusion.diffusion limitedCorrected
```

覆盖值只应用于调用时使用同名 term 的装配。独立显式运算采用 `operation.<name>.<option>`，必须由 C++ 求解器主动读取同名配置。旧式全局 `gradient`、按字段写 `T.gradient` 等配置不属于当前接口。

## 格式的选择

| 设置 | 取值与使用 |
| --- | --- |
| 插值 | `linear`；`corrected` 加入梯度修正 |
| 梯度 | `greenGauss`；`leastSquares` |
| 对流 | `upwind` 稳健但较耗散；`linearUpwind` 增加重构修正；`central` 中心格式 |
| 扩散 | `orthogonal`；`corrected` 非正交修正；`limitedCorrected` 限制修正 |
| 时间 | `steady` 无时间项；`euler`；`bdf2`，第一步以 Euler 启动 |

格式名称不能替代精度验证。强对流、非正交网格和边界附近的重构都需要检查局部过冲与关注量误差。空间格式改变时，使用新的运行名，并保持物理输入、网格和时刻一致。

## 时间步

```text
startTime 0
endTime 0.05
deltaT 0.01
```

`deltaT` 为正，结束时间不得早于起始时间。Euler 可以在最后一步截短到 `endTime`；当前配置驱动的 BDF2 要求时间区间与固定时间步匹配。完整规则见[数值配置参考](../reference/numerical-controls.md)。

PISO 日志提供 `CoMax`、`CoMean`。它们是观测值，软件没有据此自动减小 `deltaT` 或以固定 CFL 上限拒绝时间步。至少比较两组较小时间步的关键结果，尤其是瞬态载荷、涡量和周期。

## 线性求解

```text
equation.temperature.kspType bcgs
equation.temperature.pcType hypre
equation.temperature.absoluteTolerance 1e-14
equation.temperature.relativeTolerance 1e-10
equation.temperature.maxIterations 1000
```

核心五项必须逐方程提供。可用 KSP 为 `cg`、`bcgs`、`gmres`、`fgmres`；可用 PC 为 `none`、`icc`、`hypre`、`gamg`、`bjacobi`、`asm`、`jacobi`。`icc` 只适用于串行；CG 应用于满足相应对称性与正定性要求的系统，不能因为它迭代少就用于任意对流矩阵。

后端检查真实的未预条件残差 `r = b - A*x`。停止阈值为 `max(atol, rtol * max(||r0||, ||b||))`。日志相对残差的缩放与外层方程残差不同，应结合对应算法定义阅读。矩阵更新与预条件器复用是不同动作；复用 PC 不表示复用旧矩阵。

可选 `warmStart`、`preconditionerRebuildInterval`、`preconditionerReuseIterationGrowth` 和 GAMG 设置的默认值及范围见[完整参考](../reference/numerical-controls.md)。先得到稳定基准，再调整复用策略并比较迭代数与时间。

## 外层迭代与时间步接受

SIMPLE 系列同时检查线性成功、速度变化、动量/压力校正残差、质量误差及模型耦合。PISO 的一次时间推进则需要区分线性成功、质量误差、压力速度耦合和湍流耦合；`maxIterations 1` 不要求速度达到稳态。

PISO 的 `nCorrectors` 为最低压力校正次数，`maxCorrectors` 允许追加校正；只有明确设置有效的 `couplingTolerance` 才对相应耦合指标施加实际限制，默认极大值近似不限制该指标。`accepted=true` 因此只说明满足了本次设置的接受条件。

## 一个可操作的调参顺序

先核对网格与边界，再用基础离散和合理时间步得到可运行基准；随后查看线性残差、外层变化和质量误差，判断问题处于哪个层次。只在线性阶段受限时调整 KSP/PC；只在耦合阶段受限时调整对应校正次数、外层阈值或松弛。最后完成网格、时间步和格式敏感性比较，并保存每组输入。
