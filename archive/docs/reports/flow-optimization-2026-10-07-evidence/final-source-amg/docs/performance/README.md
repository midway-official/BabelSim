# BabelSim 性能测量与 PETSc 后端

性能数据必须同时说明网格、输入配置、MPI rank、硬件/绑定、构建方式和求解状态。配置 A/B 用相同算例和停止条件；
不能把更换预条件器带来的迭代下降称为代码优化，也不能以短程稳定性结果替代物理精度验证。

2026-10-07—08 的通用后端、网格分发、划分与 CPU AMG 优化，以及工程规模固定工作量的
强扩展、数值一致性和 RSS/PSS 证据，见 [工程优化报告](../reports/flow-optimization-2026-10-07.md)。

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

## 排除初始化的持续运行与二进制 A/B

`iterationComputeSeconds` 是逐 rank 的循环计时窗口：首次调用通用
`monitor::Reporter::iteration` 时启动，之后每个完整迭代更新结束时间，扣除窗口内结果写出。
它排除整个首次迭代及之前的网格读取/构造、场与方程创建、初次矩阵模式建立、初次预条件器准备。
`reportedIterationIntervals` 给出窗口内的完整迭代间隔数；SIMPLE 的 N 次外迭代通常计量 N−1 次，
单外迭代 PISO 的 N 个接受步通常计量后 N−1 步。没有至少两次迭代观测时，窗口为 0，不能用于持续运行结论。
`solverComputeSeconds` 仍保留原定义，作为包含首次迭代的求解器回调计时。

`tools/benchmark_solver_optimization.py` 顺序运行相同输入的不同二进制；每个 variant 拥有独立的
案例结果目录，所有 rank 单线程并绑定物理核。驱动使用互斥锁，避免多个该驱动同时运行。
重复样本按固定随机种子交错，一次只启动一个 mpirun；预热独立保存且不进入中位数。
SIMPLE 固定工作量要求真实达到给定外迭代数，PISO 可要求给定接受步数；失败会停止后续任务。
默认把循环计算不足 60 秒的样本标为短样本，保留原始证据，不称为持续运行。

```bash
python3 tools/benchmark_solver_optimization.py \
  --variant baseline=/tmp/babelsim-baseline/babelsim-solve \
  --variant optimized=/tmp/babelsim-optimized/babelsim-solve \
  --case simple=/tmp/flow-study/simple --throughput simple=384 \
  --case piso=/tmp/flow-study/piso --accepted-steps piso=100 \
  --ranks 1,2,4,8,12 --repeat 3 --warmup 1 \
  --output /tmp/flow-study/evidence --work-root /tmp/flow-study/runs
```

以上命令假定 SIMPLE 输入确实固定执行 384 次外迭代，PISO 输入确实执行 100 个接受步。
两种二进制都必须有相同的循环计时观测代码，且每个案例的离散、容差和停止条件一致。
不同 rank 使用相同输入，才能计算强扩展；最快时间、加速比和并行效率同时报告。
阶段计数和阶段时间仍覆盖完整回调，包括首次迭代，因此不要求与此循环窗口闭合。

强扩展加速比定义为 `S(P) = T(1)/T(P)`，效率定义为 `E(P) = S(P)/P`，
其中分子必须是同一版本、同一网格、同一工作量的单进程时间。
同时优化串行和并行时，要分别报告同进程数的版本加速以及每个版本自身的强扩展效率。
不能借用旧版本较慢的串行时间，或降低串行预条件器效率，使新版本的并行效率显得更高。
工程算例应让最快的并行点也有至少 60 秒有效计算；固定工作量下串行测试通常更长。
当前研究使用 1,024,000 单元畸变三棱柱网格和 512,000 单元非结构 Voronoi 棱柱网格，
同样执行 10 个接受的 PISO 步，计量后 9 步；小网格用于正确性回归。

该驱动按 0.1 秒采样 `mpirun` 及其子进程的累计 RSS，并每秒记录 PSS、进程交换页、
系统换入/换出及日志中已报告的工作量。`peakTreeRssKiB` 是同时累计的 RSS，包含共享页重复计数；
`peakTreePssKiB` 按比例分摊共享映射，未能读取时为 null。两者的峰值均覆盖完整进程生命周期，
包括初始化与输出，不能称为循环计算窗口内的内存峰值。日志可能缓冲，trace 的 `reportedWork`
只能帮助近似识别计算阶段。内存采样不参与求解状态或停止判断。

通用后端的 `preconditionerReuses` 计数使用旧矩阵版本预条件器的逻辑 RHS 求解，
`preconditionerRetries` 计数其中需要重建后重试的次数。`linearSolves` 计数逻辑 RHS，
`krylovIterations` 和 `trueResidualChecks` 计入两个尝试；重建时间计入 `preconditionerSeconds`。
Block Jacobi/ASM 的局部因子准备现在在显式 setup 阶段计时，避免计入首个 RHS 的 solve 时间。

## 通用划分和内存存储

MPI 读取网格时，rank 0 读入完整网格并计算一次全局所有权，然后逐 rank 分发 owned 单元及
stencil 所需 halo 的几何与拓扑。工作 rank 构造自己的局部 Mesh；构造前释放传输缓冲，
root 同时释放完整网格。读取与构造阶段的错误在进入下一阶段前传递给所有 rank。
这个路径保持全局 ID、物理 patch、face 方向和至少三层 ghost 的契约。

可选 METIS 以通信量为目标，同时约束单元数与单元面关联数；连通图请求连续分区，
随后恢复精确 owned 单元配额。选择依据及适用范围见
[Python 网格图模型](partition-model.md)和[实验报告](../reports/partition-model-2026-10-07.md)。
没有 METIS 时继续使用确定性的图生长实现。依赖和静态库链接说明见[构建说明](../../README.md)。

PETSc 装配计划只长期保存 upper/lower face 索引；COO 行列是建立矩阵模式时的临时数据。
owned 单元顺序直接使用 Mesh 中的计划。矩阵值缓存通过缓冲交换更新，避免重复整段复制。
三个向量分量仍顺序求解，共用一对 RHS/解 Vec 和工作数组，保留各分量的收敛记录。
这些存储优化适用于现有标量/向量方程，物理离散与真实残差目标保持原要求。

## PETSc 后端生命周期

每个方程保留长期存在的 PETSc Mat、Vec、KSP 和 PC。稀疏模式建立后，矩阵数值更新、RHS 更新和求解分开执行；
同一时间步内重复压力修正可以只换 RHS 并复用已装配矩阵。`matrixPatternBuilds`、`matrixValueUpdates`、
`rhsOnlySolves`、`linearSolves`、`krylovIterations`、`preconditionerSetups` 和 `trueResidualChecks` 是当前已实现的
后端计数器。`solverComputeSeconds` 从求解器回调时间中扣除已计量结果写出时间。

PETSc 线性设置由 `numerics/solution.bs` 的 `equation.<name>.kspType`、`pcType`、容差、迭代上限及可选 warm start
读取。`pcType hypre` 使用 PETSc 的 Hypre 接口配置 BoomerAMG；`pcType gamg` 使用 PETSc 自带 GAMG，
`amgMaxLevels`、`amgCoarseSize`、`amgSmoothingSteps` 可调 GAMG。CG 用于对称正定系统，必须与保持对称正定性质的
预条件配置搭配。最终应以真实残差和守恒检查确认求解质量。

Hypre 的 CPU 层级默认使用 HMIS、ext+i、每行最多 4 个插值项。
CG、GMRES、FGMRES 使用 1 层 aggressive coarsening；BiCGStab 保留常规粗化和 CF 顺序，
避免在已验证的扩散算子上发生 shadow-residual breakdown。
其他 Krylov 方法串行用 CF 松弛顺序，MPI 用自然顺序，保留对称 SOR/Jacobi 平滑和一次 V-cycle。
高级选项通过 PC 私有的临时 options 数据库应用，不修改全局 PETSc 配置。
该选择覆盖框架中所有 Hypre 方程；具体性能取决于矩阵和分区，需要用固定精度的实际运行验证。

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
