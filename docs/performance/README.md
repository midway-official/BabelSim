# BabelSim 性能测量与 PETSc 后端

性能数据必须同时说明网格、输入配置、MPI rank、硬件/绑定、构建方式和求解状态。配置 A/B 用相同算例和停止条件；
不能把更换预条件器带来的迭代下降称为代码优化，也不能以短程稳定性结果替代物理精度验证。

## 获取计时和性能 JSON

默认构建目录是 `build-petsc/`，可通过 `BUILD=<目录>` 改变。使用 `-performance` 后，应用为每个 MPI rank
写一份 `rank-XXXX.json`，并在求解器日志写汇总。例：

```bash
make -j18
mpirun -np 2 build-petsc/babelsim-solve -case cases/naca0012 \
  -time timing -performance /tmp/babelsim-naca-performance
```

长时间性能比较使用 `tools/benchmark_backend.py`。`complete` 模式保留算例原本的停止条件；驱动对每个 rank
执行 warmup 和重复样本，记录环境、网格/算例哈希、二进制哈希、绑定信息、退出状态、墙钟时间和性能 JSON：

```bash
python3 tools/benchmark_backend.py \
  --case cases/naca0012 --binary build-petsc/babelsim-solve \
  --ranks 1,2,4 --warmup 1 --repeat 3 --mode complete \
  --output /tmp/babelsim-evidence/naca0012-amg
```

正式比较应检查所有样本成功且数值结果满足同一精度条件，再比较关键路径时间。MPI 时间取最慢 rank；不要把各 rank 的阶段时间直接相加。

## PETSc 后端生命周期

每个方程保留长期存在的 PETSc Mat、Vec、KSP 和 PC。稀疏模式建立后，矩阵数值更新、RHS 更新和求解分开执行；
同一时间步内重复压力修正可以只换 RHS 并复用已装配矩阵。`matrixPatternBuilds`、`matrixValueUpdates`、
`rhsOnlySolves`、`linearSolves`、`krylovIterations`、`preconditionerSetups` 和 `trueResidualChecks` 是当前已实现的
后端计数器。`solverComputeSeconds` 从求解器回调时间中扣除已计量结果写出时间。

PETSc 线性设置由 `numerics/solution.bs` 的 `equation.<name>.kspType`、`pcType`、容差、迭代上限及可选 warm start
读取。`pcType hypre` 使用 PETSc 的 Hypre 接口配置 BoomerAMG；`pcType gamg` 使用 PETSc 自带 GAMG，
`amgMaxLevels`、`amgCoarseSize`、`amgSmoothingSteps` 可调 GAMG。CG 用于对称正定系统，必须与保持对称正定性质的
预条件配置搭配。最终应以真实残差和守恒检查确认求解质量。

以下计数器存在边界，解释性能时要特别小心：

- `globalReductions` 只统计 BabelSim 显式归约，不包括 PETSc KSP 或 Hypre 内部 MPI 归约。
- `haloBytes`、`haloExchanges` 统计 BabelSim 的场同步通信，不包括 PETSc/Hypre 矩阵、向量、KSP/AMG 内部通信。
- `sparseMatvecs`、`preconditionerApplications` 还没有覆盖 PETSc 内部调用；其当前值不能用于推断实际 SpMV 数、PC 应用数或总通信量。
- `linearSolveSeconds` 包含求解内部的 SpMV 与通信，所以不要与这些内部阶段时间重复相加。比较 MPI 性能应记录完整墙钟时间和各 rank 关键路径。

## NACA0012 全尺寸网格 AMG 短程对照

2026-09-24 的受控启动步对照使用 127,013 单元网格、2 rank、`OMP_NUM_THREADS=1`、隐式 Euler、`deltaT=0.5 s`，
所有配置只运行 1 步。压力修正的求解/预条件选择及计时见
[`cases/naca0012/validation/petsc_one_step_timing.md`](../../cases/naca0012/validation/petsc_one_step_timing.md)。
`CG + Hypre BoomerAMG` 在两次样本中的平均 `solverComputeSeconds` 为 1.745 s，单次 `BCGS + Block Jacobi`
基线为 9.058 s。这个差异混有配置选择效应，而且基线只有一个样本；它不是通用后端加速或 MPI 扩展性结论。
当前正式 NACA 压力设置据此选择 `CG + Hypre BoomerAMG`。

该正式配置的 5 步短算在 2 rank 下完成：全局 127,013 个单元字段完整且有限，没有失败求解，最大质量误差约
`6.33e-11`。端到端用时 19.27 s，其中 rank 0 的五步 `solverComputeSeconds` 为 8.63 s；本次输出耗时约 6.87 s。
短算只证明当前启动阶段可运行，不说明时间统计收敛、升阻力准确、网格/时间步独立或物理模型已经验证。

完整复现和单步/五步原始数据见
[`NACA0012 PETSc AMG 计时记录`](../../cases/naca0012/validation/petsc_one_step_timing.md)。
