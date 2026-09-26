# PETSc 后端 NACA0012 AMG 启动步计时

日期：2026-09-24。这里记录全尺寸网格上的短程性能对照和稳定性检查，不代表稳态收敛或湍流模型物理验证。

## 测试条件

- 网格：127,013 个单元，SHA-256：`ef58e557a8fee0423e193bc5bcd715aec466c562b435c1cc03f07f267f54fe25`。
- 软件：PETSc 3.25.5、Open MPI 4.1.2；2 个 MPI rank，`OMP_NUM_THREADS=1`。
- 线性求解设置比较均从同一份算例复制到独立临时目录，单步隐式 Euler，`deltaT=0.5 s`，结果未写回正式算例。
- 基线为迁移后的 `BCGS + Block Jacobi`；比较 `BCGS + Hypre BoomerAMG`、`CG + Hypre BoomerAMG` 和 `CG + PETSc GAMG`。相同压力容差、初值复用和迭代上限。
- 原始单步基线目录为 `/tmp/babelsim-naca-amg-ab/smoke-20260924-191436`；Hypre/GAMG 对照目录见 `/tmp/babelsim-naca-amg-ab/` 下相应时间戳的运行目录。

## 单步 A/B 结果

计时使用每 rank 性能 JSON 的 `solverComputeSeconds`（从求解器回调扣除结果写盘时间）。压力 KSP 时间和 Krylov 迭代数是整个时间步内 16 次压力修正求解的累计值。Hypre 两种组合各重复两次；原 BJacobi 和 GAMG 结果各为一次探索性样本。因此这些数字只用于本算例、本分区下的配置筛选。

| 压力求解器/预条件器 | 重复数 | 求解器计算时间 | 16 次压力求解累计 KSP 时间 | 压力 Krylov 迭代总数 |
|---|---:|---:|---:|---:|
| BCGS + Block Jacobi（基线） | 1 | 9.058 s | 8.117 s | 3,844 |
| BCGS + Hypre BoomerAMG | 2 | 1.981 s（均值） | 0.988 s（均值） | 53（均值） |
| CG + Hypre BoomerAMG | 2 | 1.745 s（均值） | 0.755 s（均值） | 94（均值） |
| CG + PETSc GAMG | 1 | 3.745 s | 2.550 s | 552 |

四种配置的单步场检查均通过，质量误差约 `6.3e-11`。在本次小样本中，`CG + Hypre` 的平均求解器计算时间最低，比 `BCGS + Block Jacobi` 的单次基线低约 81%；这个比例混有配置调优效应，且基线没有重复样本，不能作为通用代码加速或 MPI 可扩展性结论。GAMG 在当前网格与 2-rank 设置下慢于 Hypre。

据此，正式 NACA `solution.bs` 选择压力修正 `KSPCG + PCHYPRE/BoomerAMG`，并启用压力初值复用。上表单步 A/B 的压力相对容差为 `1e-7`；随后 MPI 场敏感性对照显示，收紧线性容差能显著减小 1/2-rank 场差，因此当前正式配置使用压力/动量 `rtol=1e-9`、k/omega `rtol=1e-10`；绝对容差分别为动量 `1e-12`、压力 `1e-14`、湍流 `1e-18`，压力最多 2400 次迭代。CG 要求对称正定系统及兼容的预条件器；本次配置在 NACA 压力方程的检查中通过，但更广泛的矩阵与网格情形仍需单独验证。动量、k 和 omega 维持 BCGS + Block Jacobi。

## 初始 AMG 配置的五步短算

- 运行目录：`/tmp/babelsim-naca-amg-pilot/smoke-20260924-191759`。
- 配置：正式 NACA `solution.bs`，2 rank，`deltaT=0.5 s`，5 步。
- 结果：退出码 0，5/5 步完成，所有线性求解通过；全局 127,013 个单元完整且唯一，字段有限，无湍流下限裁剪，最大质量误差 `6.3322e-11`。
- 端到端运行时间 `19.266 s`。性能 JSON 的 rank 0 `solverComputeSeconds=8.628 s`，五步共 105 次线性求解、559 次 Krylov 迭代，其中压力修正 80 次求解、475 次迭代、累计 KSP 时间 `4.009 s`。rank 0 输出时间约 `6.871 s`，因此端到端墙钟不能直接当作纯求解耗时。

复现命令：

```bash
env OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
  python3 cases/naca0012/run_smoke.py --steps 5 --ranks 2 --dt 0.5 \
    --performance --runs-directory /tmp/babelsim-naca-amg-pilot
```

此处 `PISO converged=true` 和 smoke 通过只说明启动短程的线性求解、守恒和字段完整性检查成功；5 步不足以说明时间统计收敛、网格/时间步独立、升阻力准确或物理模型得到验证。

## 正式容差下的 NACA 1/2-rank 对照

使用上述当前 `solution.bs` 容差，另将全尺寸 NACA 网格分别运行 1 和 2 rank，各 5 步（`deltaT=0.5 s`）。两次均完成 5/5 步、无 inexact 线性求解；每个保存时刻及 final 目录都用 `tools/compare_parallel_results.py --atol 5e-6 --rtol 5e-6` 通过。两种 rank 的网格哈希一致。

| 场 | 全局相对 L2 差最大值（跨五步） | 最大绝对差 |
|---|---:|---:|
| U | `1.20e-10` | `1.31e-8` |
| p | `1.72e-10` | `1.24e-8` |
| k | `1.72e-9` | `1.48e-9` |
| omega | `3.79e-10` | `7.77e-4` |
| mut | `4.11e-9` | `1.02e-12` |

相对 L2 为 `||x₁-x₂||₂ / ||x₁||₂`，由全局 cell ID 对齐后计算；各字段逐点最大绝对差也一并列出，避免单看全局范数掩盖局部差异。较松的原始容差（压力/动量 `1e-7`/`1e-8`、湍流 `1e-8`）下，五步场仍有限且质量守恒，但逐点 `5e-6` mixed-tolerance 检查在 omega 上未通过；收紧到当前设置后所有时刻通过。2-rank 五步 `solverComputeSeconds` 从约 `8.63 s` 增到 `9.27 s`（约 7.4%），换来明显更紧的 rank 一致性。最终压力仍由 CG + Hypre BoomerAMG 求解。

严格容差对照目录：`/tmp/babelsim-naca-amg-tight-i5uzyk8t/rank-1` 与 `rank-2`。五步对照命令示例：

```bash
python3 tools/compare_parallel_results.py \
  /tmp/babelsim-naca-amg-tight-i5uzyk8t/rank-1/results/final \
  /tmp/babelsim-naca-amg-tight-i5uzyk8t/rank-2/results/final \
  --atol 5e-6 --rtol 5e-6
```

## 计数器解释和限制

每方程长期保留 PETSc 矩阵与 KSP/PC；同一时间步中 4 次矩阵模式建立、20 次矩阵数值更新以及 75 次 RHS-only 求解体现了矩阵更新与重复求解分离。该生命周期减少重复装配并复用压力矩阵结构。

当前 JSON 的 `globalReductions` 只数 BabelSim 显式归约，不包含 PETSc KSP/Hypre 内部归约；`sparseMatvecs` 与 `preconditionerApplications` 也尚未覆盖 PETSc 内部调用，不能据此推断总通信次数或 SpMV 数量。`haloBytes` 只统计 BabelSim 的场同步交换，不包括 PETSc/Hypre 内部通信。后续要分析 AMG 的 MPI 通信成本，需要给 PETSc/Hypre 阶段增加直接计时/事件采样，并按 rank 汇总关键路径。

历史 100 步数据来自旧 Eigen 数值后端，与当前 PETSc/Hypre 配置不是等价对照；不能直接比较并宣称性能提升或退化。
