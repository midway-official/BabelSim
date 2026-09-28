# NACA0012 PISO / k–ω 长时失稳审计（2026-09-28）

## 文档导航

- [结论与正式运行验收](#当前结论与验证范围)
- [原始故障证据](#1-原始故障可以确定性复现)
- [同起点对照试验](#2-单变量定位试算)
- [PISO、非正交通量、DSL 与 RANS 审计](#3-实现审计问题在哪一层)
- [回归检查与参考资料](#4-已完成的检查)

## 当前结论与验证范围

已找到可复现的主要失稳路径：**单步内 PISO 压力—速度校正不足，误差先在压力校正中增长，随后污染湍流生产项和 k/ω**。不能把问题简单归因于 AMG、k–ω 初始场或会话中断。

本次修复包含总压力 PISO 方程、可配置耦合残差验收/追加校正、真实 Courant 数，以及湍流限幅前最小值诊断。NACA 配置保留 Δt=0.01 s、BDF2、动量 linearUpwind、k/ω upwind、Hypre AMG、4 MPI ranks、U 松弛 1、k–ω 模型，至少 8 次 PISO 校正，`couplingTolerance=1e-5`，最多 50 次。

**同一 t=8 s 保存场的 8 次校正对照已正常退出至 t=9.5 s；从原始初始场到 30 s 的正式重跑现已正常完成。** 正式进程退出码为 0，状态记录为 `dt=0.01 s`、`endTime=30 s`、4 ranks、耗时 4,352 s。最终元数据和保存场已逐 rank 核验：时间均为 30 s，覆盖原网格全部 127,013 个全局单元，字段完整且有限。

## 正式运行最终结果（2026-09-28）

从原始初始场完成了 `t=0–30 s` 正式计算。`status.json` 记录 `exit_code=0`、`wall_seconds=4352`、`deltaT=0.01`、`endTime=30`、`ranks=4`。日志共 3,000 个时间步，全部 `accepted=true`、`linear=ok`、`converged=true`；没有拒绝步或线性求解失败。全程最大 `rCoupling=9.1766e-6`，低于配置的 `couplingTolerance=1e-5`。按实际面通量定义的 CoMax 启动峰值为 `239.899`，终点 CoMax=`46.4556`、CoMean=`0.316034`。

最终 `final` 快照四个 rank 的元数据均为 `time=30`、`ranks=4`、`global_cell_count=127013`。逐字段合并后，每个字段都恰有 127,013 个唯一全局单元 ID，所有值均有限。终点范围：速度模 `6.41e-6–2.88319 m/s`，`p=−4.42244–0.511471`，`k=4.46907e-8–0.101278`，`omega=0.636935–124803 s^-1`，`mut=5.27407e-13–9.00291e-4`。最终步限幅前 `kMinRaw=4.46907e-8`、`omegaMinRaw=0.636935`，`omegaBoundedCells=0`。因此这次正式重跑通过了本审计覆盖的 PISO/k–ω 长时运行和终点保存场检查；它不替代升阻力、网格收敛或实验数据验证。

复核入口：正式日志为 [`piso_audit_20260928/formal/solver.log`](piso_audit_20260928/formal/solver.log)，
退出状态为 [`piso_audit_20260928/formal/status.json`](piso_audit_20260928/formal/status.json)。
完整按时间输出的 CSV 场保存在本地结果目录 `../results/piso_komega_bdf2_dt0p01_T30_mpi4_coupled/`；
该大体积结果目录不纳入版本库，本文记录了终点元数据和逐字段完整性检查结果。

## 1. 原始故障可以确定性复现

原始两次 dt=0.01 s、3 次 PISO 校正的日志：

- `../results/piso_komega_bdf2_dt0p01_T30_mpi4/solver.log`
- `formal_dt0p01_T30_mpi4_retry1/solver.log`

它们在抽查步上数值一致；失稳在人工中断前已经发生。

| 时间 / s | 压力增量指标 dP | 最终耦合残差 | 湍流初始残差 rTurb |
|---:|---:|---:|---:|
| 5.0 | 0.0128877 | 3.49735e-7 | 0.00181952 |
| 7.0 | 0.0315369 | 8.35041e-7 | 0.00114897 |
| 7.5 | 0.0687785 | 1.80144e-6 | 0.00105104 |
| 8.0 | 0.176166 | 4.71051e-6 | 0.000988587 |
| 8.4 | 0.749826 | 2.78787e-5 | 0.00378784 |
| 8.5 | 4.64701 | 4.93468e-4 | 0.292848 |
| 8.6 | 9.53289 | 0.244592 | 0.963951 |

`t=8 s` 保存场 Umax=2.910889 m/s、kmax=0.103618、ωmax=125692.5 s⁻¹，仍有有限的有界值；到 `t=9 s` Umax≈2.0e27，虽然仍为有限浮点数，已无物理意义。`linear=ok` 与极小质量残差不能证明流场有效。

失稳前压力变化较大的位置在翼型尾缘下游近区，例如全局单元 40074，中心 `(−0.744087, 0.138611, 0.01)`。总压力形式的 3 次校正续算在 t=8.3 s 已在该处出现约 −36.2 Pa。对出口附近速度变化的初步怀疑没有得到单变量试验支持，未据此修改边界条件。

## 2. 单变量定位试算

工作目录：`piso_audit_20260928/`。全部使用原 127,013 单元网格、4 ranks、dt=0.01 s，从原故障运行的 t=8 s U/p/k/ω 读取同一初场。

注意：CSV 续算没有保存旧面通量和 BDF2 双层历史，第一步以 Euler 重建历史。因此这些试算用于同起点的因果比较，不能称为原运行的逐位连续续算。

| 试算 | 变化 | 结果 |
|---|---|---|
| baseline | 原程序，3 次校正 | 重现失稳，人工停止 |
| absolute_pressure | 仅改总压力形式，3 次校正 | 仍失稳，说明该修复单独不够 |
| upwind_probe | 总压力，动量改一阶迎风，3 次校正 | 仍失稳；非二阶迎风单一因素 |
| pressure_corrected | 总压力，压力扩散改无限幅 corrected，3 次校正 | 仍失稳；非限幅器单一因素 |
| no_stress_probe | 临时去掉显式应力余项，3 次校正 | 仍失稳；正式版本已恢复完整应力 |
| green_gauss | 所有梯度改 Green–Gauss，3 次校正 | 仍失稳；正式版保留 leastSquares |
| correctors8 | 总压力，8 次校正，其余同原配置 | 正常退出到 9.5 s；末段 Umax≈2.90，耦合残差≈1e-8 |
| `incremental_correctors8` | 原程序，仅增加至 8 次校正 | 正常退出至 9.5 s；Umax≈2.90，耦合残差≈1.6e-7 |

失败诊断用的 Umax>20 m/s 终止阈值仅存在于试验驱动中；正式求解器不截断速度，也没有引入此经验阈值。

## 3. 实现审计：问题在哪一层

### 3.1 PISO 算法层的主问题

PISO 使用动量矩阵的对角近似校正速度。有限次校正后，即使压力方程使面通量守恒，校正速度仍可能不满足冻结的动量方程。应检查：

`rCoupling = ||b0 − A U − V grad(p)|| / (||b0|| + ||A U|| + ||V grad(p)||)`。

原程序虽输出该指标，但验收时间步时未使用它，故失稳仍报 `converged=true`。本次增加 `couplingTolerance` 和 `maxCorrectors`：至少执行 `nCorrectors` 次，不达标则追加；超过上限返回退出码 2，不写出失败时间步。未设置此可选容限的旧算例保留固定校正次数行为。

实际面通量诊断表明，dt=0.01 s 的启动阶段 CoMax 可达约 240；t≈0.4–1 s 时仍约 50，而体积加权 CoMean 约 0.31。细小近壁网格使局部时间尺度远小于全局弦长对流时间。BDF2 的隐式时间离散不保证有限次 PISO 分裂校正足够准确。

### 3.2 总压力与非正交修正的一致性

原形式将旧压力 p 与压力增量 dp 的限幅非正交通量分别计算后相加。对于 `limitedCorrected`，离散通量一般不满足 `L(p+dp)=L(p)+L(dp)`。

现改为：每次 PISO 校正重新构造 HbyA，构造 `phiH = flux(HbyA) + rAUf * temporalCorrection`，求解总压力方程，然后使用该方程冻结的 `equ::faceFlux` 更新 phi，并用同一总压力更新 U。这样压力矩阵与输出面通量保持一致。配置名 `pressureCorrection` 保持兼容，所绑定未知量现在是总压力 p。

势流初始化仍有独立的齐次边界势函数，不覆盖初始压力。速度与 BDF2 面通量历史权重已核对官方实现，未发现本次故障需要修改其系数。

### 3.3 DSL 和 RANS

对 `equ::ddt` 的 Euler/变步长 BDF2 系数、体积因子、隐式耗散源项符号、矩阵作用、冻结扩散面通量做了针对性检查。本次没有证据证明 DSL 的基础装配符号或 k–ω 方程公式是确定性失稳的根因。

k–ω 继续使用 Wilcox1988m：先解 ω，再以新 ω 作 k 的隐式耗散，生产和扩散系数按当前分步状态计算。模型方程、系数、初始场和边界条件未被更换。本次只新增限幅前 `kMinRaw` / `omegaMinRaw` 诊断；启动 BDF2 的少量 k 负值仍由原有下限处理，因此后续仍需监测，不能据此宣称严格保正。

## 4. 已完成的检查

- `physics_contract_test`：1 / 4 ranks；含倾斜网格上冻结面通量、Courant 数、有限幅压力通量不满足叠加的回归检查。
- `procedural_equation_test`：源项/体积/符号、矩阵代数、变步长 BDF2、扩散及矢量方程通过。
- `piso_coupling_test.py`：1 / 2 / 4 ranks；线性求解成功且守恒但耦合残差不合格时拒绝；增加校正可恢复；串并行速度一致性通过。
- `rans_validation_test.py`：SA / k–ω / k–ε 方程验证、Euler/BDF2 衰减阶次、SIMPLE/PISO 串并行及线性失败/限幅相关检查通过。k–ω BDF2 减半步长误差比约 3.96、3.78。
- `architecture_test.py` 通过；结构检查已允许标准总压力 HbyA 形式，不再强制出现旧增量写法的特定字符串。

这些是离散方程与集成检查；不等于 NACA0012 升阻力、失速或网格/时间精度的实验验证。

## 5. 参考资料

1. Issa, Gosman & Watkins (1986), *The computation of compressible and incompressible recirculating flows by a non-iterative implicit scheme*, JCP. [原文页面](https://www.sciencedirect.com/science/article/pii/0021999186901002)。论文讨论 PISO 分裂误差与时间步的关系；不能由隐式格式推断任意 Co 下少数校正必然足够。
2. OpenFOAM Foundation 10 [`pisoFoam/pEqn.H`](https://cpp.openfoam.org/v10/incompressible_2pisoFoam_2pEqn_8H_source.html)：总压力方程、HbyA、ddtCorr 与压力方程面通量的实际组合。
3. OpenFOAM Foundation 10 [`backwardDdtScheme.C`](https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-10/master/src/finiteVolume/finiteVolume/ddtSchemes/backwardDdtScheme/backwardDdtScheme.C) 与 [`ddtScheme.C`](https://raw.githubusercontent.com/OpenFOAM/OpenFOAM-10/master/src/finiteVolume/finiteVolume/ddtSchemes/ddtScheme/ddtScheme.C)：BDF2 历史权重及旧面通量/速度插值差的耦合修正。
4. OpenFOAM Foundation 10 [`kOmega.C`](https://cpp.openfoam.org/v10/kOmega_8C_source.html)：模型输运次序、生产与耗散项对照。
5. OpenCFD [`CourantNo`](https://doc.openfoam.com/2312/tools/post-processing/function-objects/field/CourantNo/)：本次 CoMax / CoMean 诊断使用其有限体积定义。
6. OpenCFD [压力—速度算法](https://doc.openfoam.com/2306/tools/processing/solvers/pressure-velocity/)：离散动量、面通量和压力方程的关系。
