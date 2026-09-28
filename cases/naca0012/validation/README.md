# NACA0012 验证与诊断档案

本目录同时保存当前正式验证证据和早期试算记录。阅读时请先区分正式的
`dt=0.01 s`、BDF2、PISO/k–ω 运行与旧的启动、性能 smoke；后者不能替代正式运行，
也不能直接代表当前算法配置。

## 当前正式计算

从原始初始场完成了 `t=0–30 s` 计算，使用原 127,013 单元网格、4 个 MPI rank、
PISO、Wilcox 1988 k–ω、BDF2、动量二阶 linearUpwind、k/ω 一阶 upwind 和 Hypre AMG。
3,000 个时间步全部通过耦合接受条件，最终快照时间为 30 s。

- [耦合失稳审计与最终验证结果](piso-coupling-audit-20260928.md)：根因、算法修复、
  对照试验、完整日志统计和最终保存场范围。
- [算例配置、网格和边界条件](../README.md)：当前参数及 NACA0012 建模假设。
- 正式运行日志、退出状态和对照材料：[`piso_audit_20260928/`](piso_audit_20260928/)。
- 速度场与流线图：[t=14 s](flow_t14_dt0p01.png)、[t=27 s](flow_latest_t27_dt0p01.png)、
  [终点 t=30 s](flow_t30_dt0p01.png)。
- 完整瞬态场保存在本地 [`../results/piso_komega_bdf2_dt0p01_T30_mpi4_coupled/`](../results/piso_komega_bdf2_dt0p01_T30_mpi4_coupled/)；
  大型逐时刻场文件不纳入版本库，终点元数据与字段完整性结果已摘要记录在审计中。

## 历史诊断

- [`piso-komega-diagnosis.md`](piso-komega-diagnosis.md)：早期启动场和 k–ω 更新顺序检查；
  文末补充了它与最终长时失稳根因的区别。
- [`instability_analysis.md`](instability_analysis.md)：中心/迎风离散和早期短算记录，属于
  最终配置前的历史试验。

## 历史性能与后端测试

- [`performance_analysis.md`](performance_analysis.md)：旧 `dt=0.005/0.05/0.5 s` 短算的
  求解计时和 MPI 对照，不是当前正式配置的性能报告。
- [`petsc_one_step_timing.md`](petsc_one_step_timing.md)：PETSc 后端 `dt=0.5 s` 单步及五步
  AMG 配置筛选。它用于早期线性求解器选择，不是 30 s 正式计算的性能 A/B。
- [`performance_profiles/README.md`](performance_profiles/README.md)：对应原始计数器和分区证据。

旧中心格式试算和已归档的 `dt=0.5 s` 100 步结果仅供追溯。小质量误差、`linear=ok` 或
单独的 `PISO converged=true` 都不足以证明整个流场有效；当前正式结果还经过了
`rCoupling` 接受、最终时刻及全局单元覆盖检查。
