# PISO 与 k–ω 启动阶段诊断（历史记录）

日期：2026-09-27。网格始终保持原文件，SHA256 为
`ef58e557a8fee0423e193bc5bcd715aec466c562b435c1cc03f07f267f54fe25`。

> **状态更新（2026-09-28）：** 本文记录的是早期启动短算中的初场、湍流下限和方程更新顺序
> 敏感性。它们不能单独解释后来原配置在约 8.5 s 出现的长时失稳。后续同起点对照和从原始
> 初始场开始的 30 s 正式计算见[耦合审计](piso-coupling-audit-20260928.md)；该审计将主要
> 失稳原因定位为 PISO 压力—速度校正不足，并记录了通过完整终点检查的正式结果。

## 结论

早期启动试验发现了几项独立的数值敏感性：均匀速度场在翼面无滑移约束下会产生不守恒初始面通量；ω 被截到极小固定下限时，`mut=rho*k/omega` 可能暴涨；k–ω 更新顺序也应按 ω 后 k 处理。势流初始投影、连续近壁湍流种子及黏度比限制改善了这些启动问题。后续源码审计未发现 DSL 标量输运装配或 BDF2 系数是长时失稳的直接根因；约 8.5 s 的失稳另由 PISO 压力—速度校正不足触发，细节和最终验证见[耦合审计](piso-coupling-audit-20260928.md)。

## 量化证据

- 均匀初始 U、无预投影时，首步前缘局部速度达到约 7.9 m/s；k 达数百。把压力校正从 4 增至 20 次，速度峰值仍约 9.9 m/s，说明单纯增加压力校正不是根因修复。
- 初始速度作势流投影后，初始质量不平衡降到约 `9.2e-14`，首步校正后的速度峰值约 2.0 m/s。
- 即使势流投影后，t=0.5 s 的试算仍有 2 个尾缘单元的 ω 被截到 `1e-10`，相应 `k/omega` 将湍流黏度推到约 `2.4e5 Pa·s`。给 ω 增加与 k/μ 一致的上界后，该次试算最大湍流黏度降至约 `6.7e-2 Pa·s`。
- 使用势流初始投影、连续近壁 k/ω 种子、omega 先于 k 更新、单次未欠松弛湍流校正和黏度比限制后，4 核运行到 t=1.0 s 的 20 个时间步全部通过；最大速度约 2.9 m/s，最大 k 约 `8e-2 m²/s²`。该试算仍是启动检查，不是完整计算验证。

## 耦合依据和实现改动

OpenFOAM Foundation `pisoFoam` 源码在一个物理步中先做动量预测和 PISO 压力校正，再调用湍流模型校正一次；`pimpleFoam` 则在外层流场校正后调用湍流模型。标准 k–ω 源码按 omega 方程、omega 限制、k 方程的顺序更新，k 耗散隐式项使用更新后的 omega，并将 omega 限制到由 k 和最大湍流黏度比决定的下限。

对应改动如下：

- PISO 可在推进时间历史前投影初始通量。
- 新增 `turbulenceCoupling segregated` 模式：每个物理步做一次不欠松弛 RANS 更新；必须使用 `turbulenceRelaxation=1`。瞬态变化残差单独报告，不要求被迭代到稳态容差。
- 保留 `turbulenceCoupling iterated` 作为同一时间层内的非线性收敛模式。
- k–ω 中先解 omega，再用新 omega 组装 k 的隐式耗散；使用显式配置的最大湍流黏度比保护 k/omega。
- 用有限体积算例测试覆盖顺序 BDF2 递推、湍流黏度约束、松弛参数检查和多核一致性。

参考：

- [OpenFOAM Foundation pisoFoam 时间循环](https://cpp.openfoam.org/v4/pisoFoam_8C_source.html)
- [OpenFOAM Foundation PISO 压力校正](https://cpp.openfoam.org/v4/incompressible_2pisoFoam_2pEqn_8H_source.html)
- [OpenFOAM Foundation k–ω 方程与 boundOmega](https://cpp.openfoam.org/v13/kOmega_8C_source.html)
- [OpenFOAM 压力算法控制说明](https://www.openfoam.com/documentation/user-guide/6-solving/6.3-solution-and-algorithm-control)

本报告中的启动问题和修正应与后续长时 PISO 失稳分开理解。完整 `t=30 s` 正式运行现已
以退出码 0 完成；3,000 步、耦合残差和终点场验证结果记录于
[2026-09-28 耦合审计](piso-coupling-audit-20260928.md)。
