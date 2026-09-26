# Physics 与 DSL 数学物理审查（2026-09-26）

审查基线：`005c8181222c451af68b94da277b162a9aea6a35`。先查阅下列公开原始/官方资料，再从连续方程、求解器算法、DSL 装配、离散算子和数值反例逐层核对。生产源码未修改；本次新增内容仅为审查报告及证据程序/日志。

**结论：不能认定当前全部求解器与 DSL 算子均正确。确认 5 项问题，其中 4 项 P1、1 项 P2。三种 RANS 模型的主要闭合公式基本与其声明变体一致，严重问题集中在面通量、边界上下文和瞬态算法组合。**

P1 表示会破坏守恒、边界条件、物理时间演化或接受明显失败的线性解，应优先修复；P2 表示在已支持的数值组合下损失预期精度。严重性依赖下文列出的触发条件，不意味着所有现有算例均失败。

## 参考依据

- [TMR：Spalart–Allmaras 模型及变体](https://tmbwg.github.io/turbmodels/spalart.html)：标准 SA、ft2、壁面距离与壁面条件。该地址为 [NASA 所说明的 TMR 新地址](https://www.nasa.gov/nasa-turbulence-modeling-resource/)。
- [TMR：Wilcox k–ω](https://tmbwg.github.io/turbmodels/wilcox.html)：明确区分 1988、1998、2006 和 `m` 变体。
- [CFD Direct：标准 k–ε](https://doc.cfd.direct/notes/cfd-general-principles/the-k-epsilon-turbulence-model)：输运项、产生/耗散项与有效扩散。
- [CFD Direct：压力—速度耦合](https://doc.cfd.direct/notes/cfd-general-principles/pressure-velocity-coupling)、[SIMPLE](https://doc.cfd.direct/notes/cfd-general-principles/steady-state-solution)、[瞬态/PISO](https://doc.cfd.direct/notes/cfd-general-principles/transient-solution)：动量对角响应、压力方程和面通量校正。
- [CFD Direct：法向梯度与非正交修正](https://doc.cfd.direct/notes/cfd-general-principles/surface-normal-gradient)、[二阶时间格式](https://doc.cfd.direct/notes/cfd-general-principles/second-order-time-schemes)、[欠松弛](https://doc.cfd.direct/notes/cfd-general-principles/under-relaxation)、[入口/出口混合边界](https://doc.cfd.direct/notes/cfd-general-principles/mixed-inlet-outlet-condition)。

文献用于确定数学契约；以下错误结论来自当前源码和实际运行的反例，不是仅以“与 OpenFOAM 写法不同”判错。

## F1 · P1：limitedCorrected 恢复的面通量不属于刚求解的矩阵

位置：[procedural_equation.cpp:313](/home/midway/BabelSim/src/discretization/procedural_equation.cpp:313)、[operators.cpp:441](/home/midway/BabelSim/src/discretization/operators.cpp:441)、[operators.cpp:1258](/home/midway/BabelSim/src/discretization/operators.cpp:1258)。

触发：扩散方法为 `limitedCorrected`，非正交修正被限幅，并且解在装配之后发生变化。影响 SIMPLE、transient SIMPLE、PISO 的压力面通量修正，也影响任何使用 `equ::faceFlux` 的标量扩散方程。

设一个面的正交通量为 `O(x)`，非正交通量为 `C(g)`。装配使用旧迭代值冻结：

```
C_old = clamp(C(grad(x_old)), -abs(O(x_old)), abs(O(x_old)))
```

`DiffusionSnapshot` 保存了系数、边界和旧梯度；但是 `faceFlux` 把 snapshot 中的单元值换成新解，再调用 `diffusionFlux`。限幅上界于是改成 `abs(O(x_new))`。旧梯度并不足以冻结这个非线性修正。

对于只有扩散项的积分方程，应满足：

```
V * div(equ::faceFlux(A, x)) = A*x - b_boundary
```

在 3×3×1 剪切网格、零 Dirichlet 边界上，严格线性容差的实测结果：

| 方法 | 求解后矩阵残差 L2 | 求解后面通量散度积分 L2 | 两者恒等式误差 L2 |
|---|---:|---:|---:|
| corrected | 1.999e-15 | 6.512e-15 | 6.079e-15 |
| limitedCorrected | 2.688e-15 | 8.726647 | 8.726647 |

因此这不是线性容差不足；压力方程可以已求准，恢复出的压力通量却不能消除对应连续性缺陷。[非正交修正资料](https://doc.cfd.direct/notes/cfd-general-principles/surface-normal-gradient)解释了显式修正与隐式部分的分拆；这里违反的是冻结装配与通量恢复之间的离散一致性。

建议：装配时保存实际的显式面修正通量（已经限幅），恢复时只以新解更新隐式正交通量；或者保存足以精确重建旧限幅的完整状态。回归测试必须在未知量改变/求解之后验证恒等式，而不只在装配状态验证。

证据：[operator_probes.cpp](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/operator_probes.cpp)、[运行日志](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/operator_probes.log)。

## F2 · P1：math::flux(U) 用零初始化输出覆盖 inletOutlet 的已有流向上下文

位置：[math.h:150](/home/midway/BabelSim/include/babelsim/math.h:150)、[operators.cpp:919](/home/midway/BabelSim/src/discretization/operators.cpp:919)、[fvm_execution.cpp:246](/home/midway/BabelSim/src/discretization/fvm_execution.cpp:246)。

触发：`U` 有 `inletOutlet` 边界，并且计算面通量前已经用 `setBoundaryFlux(phi)` 提供回流信息。影响所有通过返回值接口 `math::flux(U)` 计算通量的调用，尤其是 SIMPLE/transient SIMPLE，以及 PISO 的初始化和旧时间层通量重建。PISO 内部 calculated predictor 的边界已是数值迹，不能把这一特定分支也一概说成必然受影响。

`math::flux` 创建初值为零的新面场；底层 `flux()` 将该输出场写进输入 `U` 的边界上下文，又以输出原值 `previous_flux=0` 判定所有边界流向。已有的负通量因而丢失，混合边界使用零梯度分支。执行层随后还会把新结果再次保存到输入的上下文中。

两单元、面积为 1 的 `plus_x` 面：

```
单元 Ux = +1
inletOutlet 入流指定 Ux = -2
已有外向通量 phi = -1
math::interpolate(U) 调用前边界值：-2
math::flux(U) 返回：+1      （应为 -2）
math::interpolate(U) 调用后边界值：+1
```

该调用既改变流向也改变随后表达式的边界取值。[混合入口/出口边界的定义](https://doc.cfd.direct/notes/cfd-general-principles/mixed-inlet-outlet-condition)要求按已有面通量符号选择入流固定值或出流零梯度，不能让输出缓冲区的初始化值替代它。

建议：将流向上下文作为输入读取/显式传入；不要在求值前以输出场覆盖它。无既有上下文的初始化应有独立且明确的策略。返回值接口与写入已有输出场的接口应遵守相同契约。

证据同 F1。

## F3 · P1：默认单次 PISO 与 RANS 的 0.7 欠松弛组合改变物理时间演化

位置：[piso/main.cpp:80](/home/midway/BabelSim/src/physics/piso/main.cpp:80)、[piso/main.cpp:198](/home/midway/BabelSim/src/physics/piso/main.cpp:198)、[k_omega.cpp:32](/home/midway/BabelSim/src/physics/RANS/k_omega.cpp:32)、[k_omega.cpp:84](/home/midway/BabelSim/src/physics/RANS/k_omega.cpp:84)。k–ε 和 SA 有相同的默认松弛及每次调用仅解一次的结构。

触发：PISO `maxIterations=1`（默认），启用任一 RANS 模型，`turbulenceRelaxation<1`（默认 0.7）。即使每个线性方程求解到很严的容差，问题仍然存在。

方程欠松弛相当于向原系统增加：

```
E * (q_new - q_iter),  E = (1/alpha - 1) diag(A)
```

在单次瞬态调用中 `q_iter=q_old`，而 `diag(A)` 含 `rho*V/dt`，附加项不是随 `dt→0` 消失的误差。固定网格且质量项主导时，它使有效时间容量趋向 `rho/alpha`。空间扩散对角也被松弛，会使原本均匀的衰减出现与单元邻接有关的差异。`equ::relax` 本身实现符合[标准方程欠松弛公式](https://doc.cfd.direct/notes/cfd-general-principles/under-relaxation)；错误发生在将一次尚未收敛的松弛迭代直接当作物理时间步的组合方式。[瞬态算法参考](https://doc.cfd.direct/notes/cfd-general-principles/transient-solution)也讨论了时间对角增强和不必欠松弛的常规瞬态路径。

生产 `babelsim-solve` 实测无流动均匀衰减，`k(0)=omega(0)=1`、`t=0.1`。解析解为：

```
omega(t) = 1/(1 + 0.075*t)
k(t) = (1 + 0.075*t)^(-0.09/0.075)
```

所以 `k_exact=0.991073660642`。取 global_id 排序后首单元的结果：

| alpha | dt | k(t=0.1) | k 绝对误差 |
|---|---:|---:|---:|
| 0.7 | 0.01 | 0.998271441583 | 7.198e-3 |
| 0.7 | 0.0025 | 0.996297735858 | 5.224e-3 |
| 0.7 | 0.0001 | 0.993909197607 | 2.836e-3 |
| 1.0 | 0.01 | 0.991074324405 | 6.638e-7 |
| 1.0 | 0.0025 | 0.991073826626 | 1.660e-7 |

完整证据还包含全场最小/最大值，及 k–ε、SA 的对应偏差。不能仅由表中有限个 dt 宣称证明了极限；不一致性的极限结论来自上述附加时间项推导，数值试验与其相符。

建议：单次 PISO 的湍流输运采用 `alpha=1`；若用户要求小于 1，则在同一物理时间层进行足够内迭代，或明确拒绝这一不一致的配置。关闭欠松弛仅解决本项，不自动保证冻结非线性源项的 BDF2 总体二阶精度。

证据：[piso_decay_probe.py](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/piso_decay_probe.py)、[完整数据](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/piso_decay.json)。

## F4 · P1：单次 PISO 接受湍流线性方程 MaxIterations

位置：[piso/main.cpp:201](/home/midway/BabelSim/src/physics/piso/main.cpp:201)、[piso/main.cpp:212](/home/midway/BabelSim/src/physics/piso/main.cpp:212)。

触发：`maxIterations=1`，湍流输运的任意线性解达到最大迭代数但未发生 NaN/Inf 或其他数值故障。

`TransportResult::healthy()` 允许 `MaxIterations`；`result.linearConverged()` 只放进 `turbulenceConverged`；最终 `stepAccepted` 的 `linearConverged` 却只包含动量和压力结果。`maxIterations==1` 又跳过含湍流结果的 `outerSettled`，使失败的湍流线性解不影响整步接受。

生产求解器反例：k–ω 单步，`kTransport.pcType=jacobi`，`kTransport.maxIterations=1`，其他线性求解器保持严格容差：

```
kTransport.lastStatus            = maxIterations
kTransport.lastRelativeResidual  = 0.19483548594795416
总体 status                      = converged
控制台                           = linear=ok converged=true
进程退出码                       = 0
```

PISO 可不要求每个物理时间步达到稳态外迭代收敛，但这不等于可以忽略输运方程的线性失败。建议将湍流 `linearConverged()` 无条件加入时间步接受条件，并与外层变化量/非线性收敛条件分开。

证据：[复现程序](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/piso_failure_probe.py)、[控制台日志](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/piso_failure.log)、[逐方程状态与真残差](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/piso-failure-performance/rank-0000.json)。

## F5 · P2：heat/transport 的 BDF2 在显式空间修正下退化为一阶

位置：[heat/main.cpp:25](/home/midway/BabelSim/src/physics/heat/main.cpp:25)、[transport/main.cpp:27](/home/midway/BabelSim/src/physics/transport/main.cpp:27)、[operators.cpp:516](/home/midway/BabelSim/src/discretization/operators.cpp:516)、[operators.cpp:1261](/home/midway/BabelSim/src/discretization/operators.cpp:1261)。

触发：BDF2，且非正交 `corrected` 扩散或 `linearUpwind` 对流的显式修正非零。该问题不等于 BDF2 时间权重写错。

两个求解器每步保存旧值后只装配、求解一次。设空间算子分为隐式部分 `L_implicit` 和显式修正 `C`，实际执行的是：

```
D_BDF2(q_new) + L_implicit(q_new) + C(q_old) = source
```

`C(q_old)` 与应在新时间层评估的 `C(q_new)` 相差一般为 `O(dt)`，因此 BDF2 的二阶时间导数无法保住全方程二阶。DSL 按契约冻结输入是正确的；求解器必须理解这种延迟修正语义。[二阶时间格式](https://doc.cfd.direct/notes/cfd-general-principles/second-order-time-schemes)与[非正交修正](https://doc.cfd.direct/notes/cfd-general-principles/surface-normal-gradient)分别给出了新时间层离散和显式修正需要更新的背景。

证据程序采用这两个求解器相同的时间推进/装配次序，并显式传入相同的 Euler/BDF2 权重，通过公共 DSL 构造无源热扩散及常速输运，未调用生产入口；这点与 F3/F4 的生产可执行文件试验区别明确。固定空间网格，终止时间 0.2，以 20/40/80/160 步比较相邻时间分辨率的最终场 L2 差值。二阶时差值比趋向 4，一阶时趋向 2：

| BDF2 空间组合 | 前一组差值比 | 后一组差值比 | 观察 |
|---|---:|---:|---|
| orthogonal 扩散 | 4.1180 | 4.0486 | 二阶对照 |
| corrected 扩散 | 1.7591 | 1.8918 | 趋于一阶 |
| upwind 对流 | 4.0320 | 4.0188 | 二阶对照 |
| linearUpwind 对流 | 2.0493 | 2.0228 | 一阶 |

扩散两组使用同一个剪切网格；orthogonal 在该网格上的空间误差不被当作正确空间解，仅作为时间阶对照。比较是在各自固定半离散算子下进行，不涉及网格误差向时间误差的混入。

建议：增加步内修正迭代至适当容差，或者对显式修正采用时间一致的二阶外推/分裂方法，并处理可变步长与启动。未经处理的组合应明确标注整体时间阶。对单次 PISO 的冻结动量/RANS 系数也应补做独立时间阶验证，不能从 `equ::ddt` 的单元测试推断全算法二阶。

证据：[time_order_probe.cpp](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/time_order_probe.cpp)、[时间阶日志](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/time_order.log)。

## 连续模型与 DSL 对照

以下“符合”仅指列出的公式、单位和已检查路径，不覆盖前述缺陷，也不构成全部几何/参数下的验证声明。

| 对象 | 应实现的内容 | 审查结论 |
|---|---|---|
| heat | `rho*cp*dT/dt - div(k grad T) = Q` | `ddt(rho*cp)`、`laplacian(k,-1)`、`source(Q)` 的符号、体积和单位对应正确；F5 影响时间组合 |
| transport | `storage*dC/dt + div(U C) - div(D grad C)=S` | 与当前参数定义一致，`phi=U_f·Sf` 为体积通量；不是自动乘 storage 的对流方程；F2/F5 需修复 |
| 三种不可压流求解器 | `rho*dU/dt + rho*div(U U) = -grad(p) + div(tau)`，`div(phi)=0` | rho、动态黏度、压力梯度与积分矩阵单位匹配；`rAU=V/aP` 正确，不应直接写 `1/aP` |
| 稳态 SIMPLE | 动量预测、Rhie–Chow、压力增量、速度/通量修正 | 压力源 `-div(phiH)` 与负 Laplacian、`phi += faceFlux` 的符号相容；内部面区域包含 MPI 分区面；受 F1/F2 影响 |
| transient SIMPLE | 每物理时间步冻结历史，步内重装配收敛 | U/RANS 历史未在内迭代重复推进；压力/速度/通量使用 alphaP 的变体使每次内迭代不必严格守恒，但最终有质量残差门槛，不应仅因中途不守恒判错 |
| PISO | 每个 corrector 重建 H(U)，避免重复投影已经守恒的通量 | 使用 `rhs + V grad(p)` 去除压力源、`apply(A,U)` 重建非对角作用的代数正确；时间通量历史的 Euler/可变步 BDF2 权重与 ddt 相符；受 F1–F4 影响 |
| RANS 应力 | 隐式 `div(muEff grad U)` 加显式 `div(muEff*(grad U^T - 2/3 I div U))` | 没有把完整应力重复计入；张量按 `gradU[i][j]=dUi/dxj` 和 `div(T)i=dTij/dxj` 一致使用；变量黏度应力测试通过 |

SIMPLE 中 `relax` 后再整体 `scale(alphaU)` 并不改变该动量预测线性方程的解，但会使提取的对角响应不同于直接使用松弛对角的写法。源码明确采用这一归一化变体；本次未将其仅凭写法差异认定为错误。transient SIMPLE 没有 PISO 那样的旧面通量历史校正，应另做非定常压力/速度网格与时间收敛验证；本次没有足够独立证据将它列成第六项已确认问题。

### 三种湍流模型

**SA**：[spalart_allmaras.cpp](/home/midway/BabelSim/src/physics/RANS/spalart_allmaras.cpp) 对应保留 ft2、去掉 primary trip term 的标准正变量 SA；不是 SA-noft2，也不是 SA-neg。`mu_t=rho*nuTilda*fv1`、扩散 `(mu+rho*nuTilda)/sigma`、`rho*cb2/sigma*|grad nuTilda|²`、反应项正负拆分及默认常数逐项相符。`S_tilde` 正值下限属于参考资料列出的允许正则化类别，宜在用户文档中注明选择。[SA 定义](https://tmbwg.github.io/turbmodels/spalart.html)。

**k–ω**：[k_omega.cpp](/home/midway/BabelSim/src/physics/RANS/k_omega.cpp) 对应 Wilcox1988m：`mu_t=rho*k/omega`；betaStar=0.09、beta=3/40、gamma=5/9；两扩散系数为 `mu+0.5*mu_t`。k 耗散以 `rho*betaStar*omega*k` 隐式线性化，omega 耗散以 `rho*beta*omega*omega_new` 处理，产生项比例 `gamma*P*omega/k` 正确。缺少 2006 的交叉扩散/Pope 修正/黏度限制器不是这个声明版本的 bug。[Wilcox1988m 定义](https://tmbwg.github.io/turbmodels/wilcox.html)。

**k–ε**：[k_epsilon.cpp](/home/midway/BabelSim/src/physics/RANS/k_epsilon.cpp) 的 `mu_t=Cmu*rho*k²/epsilon`、`mu+mu_t/sigma_k`、`mu+mu_t/sigma_epsilon`、`P`、`C1*P*epsilon/k`、`rho*epsilon`、`C2*rho*epsilon²/k` 及默认常数相符。扩散中的 sigma 是除数，与 k–ω 使用乘数不同，当前实现没有混淆。[标准 k–ε 参考](https://doc.cfd.direct/notes/cfd-general-principles/the-k-epsilon-turbulence-model)。

两个两方程模型使用 `2*dev(symm(grad U)):dev(symm(grad U))` 计算产生率；在目标不可压连续方程 `div U=0` 下对应标准应变产生率。离散 cell-gradient 的散度未必与守恒面通量散度相同，不应把本实现推广为可压缩 RANS。

`k` 的各向同性 Reynolds 应力项未显式加入动量。Wilcox 的 `m` 命名允许此选择；对于常密度不可压形式，也可解释为吸收到修正压力。需要在压力输出/边界的物理意义中统一说明，不能不加说明地把所有 `p` 都视为已恢复的静压。本次不将这一声明变体本身列作公式遗漏错误。

### DSL 算子

- `equ::source`、`reaction`、`ddt` 乘一次单元体积；`rhs`、`diagonal`、`apply` 保持积分量；`math::div(faceFlux)` 除以体积。非单位体积测试通过。
- 可变步 BDF2 令 `r=dt/dt_previous`，LHS 权重为 `(1+2r)/((1+r)dt)`，RHS 旧层权重为 `(1+r)/dt`、`-r²/((1+r)dt)`。当前实现正确，首步以 Euler 启动；二次多项式测试通过。
- `equ::laplacian(...,-1)` 表示 LHS 的负散度扩散，正交部分形成正对角/负邻接；Dirichlet 和固定外法向梯度源项符号已核对。普通 corrected 的冻结通量一致性额外反例通过；limitedCorrected 见 F1。
- Upwind 的 owner/neighbour 配对和正负通量矩阵系数、central 权重、linearUpwind 延迟修正 RHS 符号已检查。高阶修正的时间层组合见 F5；linearUpwind 未限幅，不能宣称一般保正/TVD。
- 插值权重、LS 梯度、Green–Gauss、法向梯度、tensor divergence 与逐点转置/trace/isotropic 路径已检查并运行对应已有测试。几何采用 `c=(Sf·d)/|d|²`、`T=Sf-c*d` 的 minimum-correction 分拆；这是可用分拆，不能直接套用 over-relaxed 分拆的稳定角度结论。
- 普通 Field 赋值保留目标物理边界；派生表达式使用 calculated traces。RANS 的内部场下限赋值不会把 SA 的零 Dirichlet 壁面值整体抬高；现有测试验证了壁面有效黏度返回分子黏度。
- MPI 的体积分/范数按 owned 实体归约，Interior 面区域不等于“本 rank 内部面”。本次 RANS 在 1/2/4 rank 测试通过，但新增小反例仅串行运行，未宣称已穷尽 MPI 路径。

## 验证结果及边界

实际运行并通过的现有测试：`procedural_equation_test`、`numerical_contract_test`、`math_runtime_test`、`operators_test`、`time_history_test`、`field_boundary_test`。其中 numerical contract 包含 30 个求解器/预条件器/尺度组合。完整日志：[unit_tests.log](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/unit_tests.log)。

运行 `tests/rans_validation_test.py`，三个模型的发布公式项、均匀衰减时间阶、SIMPLE 集成、1/2/4 进程一致性以及裁剪后的拒绝收敛测试均通过。Euler 误差比约 2、BDF2 误差比约 4。该脚本的瞬态生产路径使用 transient SIMPLE 并进行步内迭代，不能覆盖默认单次 PISO 的 F3/F4。[RANS 日志](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26-evidence/rans_validation.log)。

本次没有跑全量 Makefile test/validate、完整 Ghia/NACA/圆柱高成本基准，也没有建立所有边界类型和网格畸变下的完备证明。已确认问题均有小规模可复现证据；其余路径仅给出公式和现有测试支持范围。

SA 读取用户提供的 `wallDistance`，本次未证明任一实际复杂几何算例的该字段是真正最近壁面距离。k–ε 为高 Reynolds 数形式，源码未提供通用壁函数/低 Reynolds 数壁面阻尼流程；k–ω 的 omega 壁面约束也依赖输入。现有制造/均匀衰减验证不能替代壁面摩阻、对数律、分离或自由剪切流的物理验证。此类内容列为适用范围和后续验证需求，不伪装成已复现源码缺陷。

建议修复次序：先 F1/F2 的 DSL 契约，再修复 F3/F4 的 PISO 时间步策略，随后处理 F5 的时间一致性；每项加入对应反例作为回归验收。完成这些后再做网格/时间步独立性及壁面流动验证。

## 后续修复

F1–F5 已完成最小修复，具体实现、语义约束及修复后测试见
[修复与验证报告](physics-fixes-2026-09-26.md)。本文数值和原始日志保留为修复前证据。
