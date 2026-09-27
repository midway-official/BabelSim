# PISO 与 k–ω 启动失败分析

日期：2026-09-27。网格始终保持原文件，SHA256 为
`ef58e557a8fee0423e193bc5bcd715aec466c562b435c1cc03f07f267f54fe25`。

## 结论

没有证据指向 DSL 的标量输运组装、BDF2 时间项或 AMG 基础算子错误。更直接的故障来自算例启动场与求解组织：均匀速度场在翼面无滑移约束下产生不守恒初始面通量；ω 截断到很小的固定下限后，`mut=rho*k/omega` 可暴涨；PISO 代码还先解 k、后解 ω，使首轮 k 耗散项使用旧 ω。原始一次 PISO 后再把湍流方程迭代到稳态容差，也把瞬态分步校正和稳态非线性收敛混在了一个接受判据里。

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

DSL 组装测试和 RANS 验证测试通过。完整 `t=30 s` 运行仍须由正式计算结果判定。
