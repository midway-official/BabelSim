# 06 选择物理模型

[用户手册](README.md) · 上一章：[网格与边界条件](mesh-and-boundaries.md) · 下一章：[数值设置与收敛](numerics.md)

## 求解器选择

| `case.bs` 的 solver | 必需输入场 | 主要用途 | 时间配置 |
| --- | --- | --- | --- |
| `heat` | `T` | 常物性导热、体积热源 | `steady` / `euler` / `bdf2` |
| `transport` | `C`、给定 `U` | 指定速度的标量对流扩散 | `steady` / `euler` / `bdf2` |
| `simple` | `U`、`p` | 不可压缩稳态流动 | `steady` |
| `transientSimple` | `U`、`p` | 每个时间步内迭代 SIMPLE | `euler` / `bdf2` |
| `piso` | `U`、`p` | 瞬态压力速度校正 | `euler` / `bdf2` |

`heat` 和 `transport` 的驱动循环仍按 control 时间推进；使用 `steady` 时应提供非空时间区间，时间项不参与装配。若 `endTime == startTime`，它们不会执行一次方程求解。稳态 `simple` 在自己的迭代循环收敛后写出起始时刻的结果。

## 热传导

```text
density 1
heatCapacity 1
conductivity 0.1
source 0
```

求解 `rho * cp * dT/dt - div(k * grad(T)) = Q`。密度和热容为正，导热率可为零，源项可正可负。`T` 的固定值和固定梯度分别控制边界温度和温度法向导数。使用 `equation.temperature` 设置温度方程。

## 给定流场的标量输运

```text
storage 1
diffusivity 0.01
source 0
```

求解 `storage * dC/dt + div(U * C) - div(diffusivity * grad(C)) = source`。`U.field` 是给定流场，求解器不会同时求解动量与压力；对流项不会自动再乘 `storage`。使用 `equation.transport` 设置浓度方程。需要双向耦合或温度反馈时，应使用相应的自定义求解器。

## 不可压缩流动

```text
density 1
dynamicViscosity 0.01
turbulenceModel none
```

密度与动力黏度为正。动量方程包含惯性、对流、压力梯度和黏性扩散，压力校正用于控制离散质量不平衡。`simple` 使用稳态外层迭代；`transientSimple` 在每个时间步内完成耦合迭代；`piso` 按预测与多次压力校正推进。

流动方程使用 `equation.momentum` 与 `equation.pressureCorrection`。虽然配置名相同，SIMPLE 的校正未知量与 PISO 的压力求解过程并不完全相同，不应直接依据名字推断算法。

稳态 SIMPLE 的 `velocityRelaxation`、`pressureRelaxation` 调整外层更新；瞬态 SIMPLE 的压力松弛还影响该算法的速度/通量校正。PISO 的 `velocityRelaxation` 作用于动量预测方程，取 1 时关闭这种松弛。

## RANS 模型

RANS 是流动求解器中的模型选择，不能写成 `solver kOmega`。

| `turbulenceModel` | 追加输入场 | 模型说明 |
| --- | --- | --- |
| `none` | 无 | 层流，不创建湍流输运方程 |
| `SA` | `nuTilda`、`wallDistance` | Spalart–Allmaras 单方程模型；壁距必须提供 |
| `kOmega` | `k`、`omega` | 当前实现为 Wilcox1988m，非 SST |
| `kEpsilon` | `k`、`epsilon` | 标准 k-epsilon 形式 |

模型创建 `mut` 和有效黏度，湍流动力黏度 `mut` 默认加入输出选择。新增输入场的每个 patch 都必须配置边界；新增输运方程的 methods/solution 必须按名称完整配置。模型方程名依次为 `nuTildaTransport`、`kTransport`、`omegaTransport`、`epsilonTransport`。

例如选择 kOmega，至少需要：

```text
# physics
density 1
dynamicViscosity 0.01
turbulenceModel kOmega
```

同时提供 `k.field`、`omega.field`，并配置 `equation.kTransport` 和 `equation.omegaTransport`。系数、下限和黏度比限制见[物理参数参考](../reference/physical-models.md)。当前模型不会根据入口速度自动推算湍流初值，也不提供通用壁面函数配置。

## 选择建议

首次检查配置与边界时，优先使用规模小的层流/标量算例，确认单位、几何和数据路径；再逐步引入复杂模型。非定常问题的 PISO 时间步需要通过 Courant 数、压力速度耦合误差和关注量的时间步敏感性共同判断。稳态求解器不能用于证明非定常涡结构。

关闭湍流时使用 `turbulenceModel none`，并移除对应的湍流输入、方程设置和不再存在的输出场。产品目前没有提供自动选择湍流模型或物理有效性认证；已有研究算例的结果范围由各自 README 和验证报告说明。
