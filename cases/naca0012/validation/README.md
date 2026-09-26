# 试算与根因分析档案

已有的长算运行计数和流场结果来自 PETSc 替换前的 Eigen 求解后端，不是迁移后的回归证据或性能对照。
PETSc 后端已完成一个 2-rank、全尺寸网格、`deltaT=0.5 s` 的启动步计时 smoke；设置、逐方程计时和限制见
[`petsc_one_step_timing.md`](petsc_one_step_timing.md)。完整 NACA0012 长算和物理验证仍未完成。

`instability_analysis.md` 记录中心/迎风受控对照、两次失稳长算，以及最终稳定的 100 步试算。

`summary.json`、`dt_0.01/`、`dt_0.005/`、`previous_central_dt2e-5.json` 记录旧中心格式试算。它们不是当前迎风设置的验收结果。

当前正式算例采用动量、k、omega 均一阶迎风，`dt=0.5`，速度/湍流松弛均为 0.3；压力修正启用跨时间步初值复用，
2 核完成 100 步至 `t=50`。10 个输出时刻均通过有限性、完整单元数和质量误差检查。详见 `current_run.json`、
`completed_run_summary.json`、`dt0.5_100steps.log`、`flow_t30_dt0.5.png` 与 `flow_t50_dt0.5.png`。原 `dt=0.05` 与
`dt=0.005` 记录分别归档为 `current_run_dt0.05.json`、`completed_run_summary_dt0.05.json` 及对应 0.005 文件。
清理掉旧 smoke、失败场和冗余 final 副本，具体路径见 `cleanup_20260924.json`；其余有效时间步对照结果保留。
小质量误差或 `PISO converged=true` 本身不证明场解健康；本次也不构成物理验证。

`performance_analysis.md` 记录 100 步耗时诊断、MPI 图分区修正、2/4/8 核扩展对比、减少非正交修正的短算例 A/B，以及压力方程 ILUT/AMG 对比。原始计数器和各 rank 分区数据在 `performance_profiles/`。
