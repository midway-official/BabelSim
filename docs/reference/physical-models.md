# 物理模型参数参考

[配置参考](README.md) · [模型使用说明](../user-guide/models.md)

下表给出当前源码读取的键。标记必需的参数没有通用默认值；正数均要求有限且大于零，非负数要求有限且不小于零。单位由算例统一约定，软件不自动换算。

## 热传导 heat

| physics 键 | 必需/默认 | 含义 |
| --- | --- | --- |
| `density` | 必需，正数 | 密度 rho |
| `heatCapacity` | 必需，正数 | 比热容 cp |
| `conductivity` | 必需，非负 | 导热率 k |
| `source` | 必需，有限 | 体积源 Q，可正可负 |

输入 `T`，方程为 `rho*cp*dT/dt - div(k*grad(T)) = Q`，配置名 `temperature`。当前内置热求解器读取常数物性；变物性可通过公开 Field 系数接口开发。

## 标量输运 transport

| physics 键 | 必需/默认 | 含义 |
| --- | --- | --- |
| `storage` | 必需，正数 | 时间项系数 |
| `diffusivity` | 必需，非负 | 扩散系数 |
| `source` | 必需，有限 | 标量体源 |

输入 `C` 与 `U`，配置名 `transport`。给定速度只用来产生对流通量，当前求解器不会随标量反过来更新速度。

## 流动 simple / transientSimple / piso

| physics 键 | 必需/默认 | 含义 |
| --- | --- | --- |
| `density` | 必需，正数 | 流体密度 |
| `dynamicViscosity` | 必需，正数 | 分子动力黏度 |
| `turbulenceModel` | `none` | none / SA / kOmega / kEpsilon |

输入 `U`、`p`，方程配置名 `momentum`、`pressureCorrection`。压力使用物理压力含义。模型名称经过字母数字规范化处理：`laminar/off` 对应 none，`spalartallmaras` 对应 SA，`wilcox1988` 对应当前 kOmega 实现，`standardkepsilon` 对应 kEpsilon；新算例建议使用表中标准写法。

## RANS 共用控制

| 配置位置与键 | 默认 | 约束 |
| --- | --- | --- |
| solution `turbulenceRelaxation` | `0.7` | `(0, 1]` |
| solution `turbulenceTolerance` | `1e-6` | 正数 |
| physics `maxTurbulentViscosityRatio` | `1e5` | 正数 |

这些模型控制仅在相应模型启用并读取时使用。模型会限制输运量与黏度，工程验证应记录限幅影响，而不能将限幅后数值有限视为物理充分验证。

## kOmega：Wilcox1988m

需要 `k.field`、`omega.field`；方程名 `kTransport`、`omegaTransport`。湍流黏度以 `rho*k/omega` 为基础并施加模型限制。

| physics 键 | 默认 |
| --- | --- |
| `kMin` | `1e-12` |
| `omegaMin` | `1e-12` |
| `kOmegaBetaStar` | `0.09` |
| `kOmegaBeta` | `3/40`，数值 `0.075` |
| `kOmegaGamma` | `5/9`，数值约 `0.5555555555555556` |
| `kOmegaSigmaK` | `0.5` |
| `kOmegaSigmaOmega` | `0.5` |

所有表中系数与下限均要求正数。输入字典应写数值，不能写 `3/40` 或 `5/9` 表达式。该模型不是 k-omega SST，没有 SST 混合函数参数。

## kEpsilon

需要 `k.field`、`epsilon.field`；方程名 `kTransport`、`epsilonTransport`。湍流黏度以 `rho*Cmu*k*k/epsilon` 为基础并施加模型限制。

| physics 键 | 默认 |
| --- | --- |
| `kMin` | `1e-12` |
| `epsilonMin` | `1e-12` |
| `kEpsilonCmu` | `0.09` |
| `kEpsilonC1` | `1.44` |
| `kEpsilonC2` | `1.92` |
| `kEpsilonSigmaK` | `1` |
| `kEpsilonSigmaEpsilon` | `1.3` |

以上均为正数。近壁条件必须与所研究问题及当前模型实现匹配，软件不会从网格 patch 的 `wall` 标签自动生成完整湍流壁面条件。

## SA

需要 `nuTilda.field` 和 `wallDistance.field`；方程名 `nuTildaTransport`。壁距作为几何输入提供，不会由本求解器自动计算。模型由 `nuTilda` 和黏性修正产生湍流黏度。

| physics 键 | 默认 |
| --- | --- |
| `saNuTildaMin` | `1e-14` |
| `saWallDistanceMin` | `1e-12` |
| `saCb1` | `0.1355` |
| `saCb2` | `0.622` |
| `saSigma` | `2/3` 的数值 |
| `saKappa` | `0.41` |
| `saCw2` | `0.3` |
| `saCw3` | `2` |
| `saCv1` | `7.1` |
| `saCt3` | `1.2` |
| `saCt4` | `0.5` |
| `saCw1` | `saCb1/saKappa^2 + (1+saCb2)/saSigma` 按已读取系数计算 |

以上均要求正数；默认公式在源码中计算，用户自定义时写最终数值。SA 与两方程模型的输入量、单位和边界并不互换。

## 自定义教程的参数

[耦合标量示例](../examples/coupled_scalar.cpp)使用自定义 `solver coupledScalar`，physics 必需键为非负 `diffusivity` 和非负 `coupling`。这两个键不因此成为内置 flow/heat 求解器的通用参数。
