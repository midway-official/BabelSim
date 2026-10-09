# BabelSim 性能测量与 PETSc 后端

性能数据必须同时说明网格、输入配置、MPI rank、硬件/绑定、构建方式和求解状态。配置 A/B 用相同算例和停止条件；
分别注明实现改进、通用后端策略改变和案例参数改变；短程稳定性不能替代物理精度验证。

2026-10-07—08 的通用后端、网格分发、划分与 CPU AMG 优化，以及工程规模固定工作量的
强扩展、数值一致性和 RSS/PSS 证据，见 [工程优化报告](../reports/flow-optimization-2026-10-07.md)。
以 Re3900 圆柱为基准的范数合并、梯度复用、真实残差续算及 PISO 湍流内校正研究，
见 [圆柱优化报告](../reports/cylinder-optimization-2026-10-09.md)。该报告标明本次验证状态。

## 获取计时和性能 JSON

默认构建目录是 `build-petsc/`，可通过 `BUILD=<目录>` 改变。使用 `-performance` 后，应用为每个 MPI rank
写一份 `rank-XXXX.json`，并在求解器日志写汇总。例：

```bash
make -j4
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
多外迭代 PISO 的计时单位是流动外迭代：各已完成物理时间步的最终迭代编号之和减一
应等于 `reportedIterationIntervals`。计时只排除首次外迭代，不能说成排除了整个首物理步。
驱动分别记录 `workUnits`（接受物理步）和 `outerIterations`（流动外迭代）。
`solverComputeSeconds` 仍保留原定义，作为包含首次迭代的求解器回调计时。

`tools/benchmark_solver_optimization.py` 顺序运行相同输入的不同二进制；每个 variant 拥有独立的
案例结果目录，所有 rank 单线程并绑定物理核。驱动使用互斥锁，避免多个该驱动同时运行。
重复样本按固定随机种子交错，一次只启动一个 mpirun；预热独立保存且不进入中位数。
SIMPLE 固定工作量要求真实达到给定外迭代数，PISO 可要求给定接受步数；失败会停止后续任务。
默认把循环计算不足 60 秒的样本标为短样本，保留原始证据，不称为持续运行。
长期研究应把输入、冻结的二进制、工作目录和证据保存在持久磁盘目录中，避免 `/tmp` 清空导致中断。
可指定 `--field-archive-root <持久目录>`，在每个求解进程退出后复制并校验导出的 owned CSV 与字段元数据。
复制发生在计时结束后；逐项清单写在对应测点的 `fields.json` 中。
旧测点的原始字段已经丢失时会明确标记不可用，已有计时记录不能替代字段一致性验证。

```bash
python3 tools/benchmark_solver_optimization.py \
  --variant baseline=/home/midway/flow-study/bin/baseline \
  --variant optimized=/home/midway/flow-study/bin/optimized \
  --case simple=/home/midway/flow-study/simple --throughput simple=384 \
  --case piso=/home/midway/flow-study/piso --accepted-steps piso=100 \
  --ranks 1,2,4,8,12 --repeat 3 --warmup 1 --timeout 7200 \
  --output /home/midway/flow-study/evidence --work-root /home/midway/flow-study/runs \
  --field-archive-root /home/midway/flow-study/fields
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
三棱柱执行 10 个接受的 PISO 步并计量后 9 步，Voronoi 执行 20 步并计量后 19 步；
每种网格的工作量在全部版本、进程数和重复样本上保持一致，小网格用于正确性回归。
若一个正式测点不足最低持续时间，应统一增加该算例在所有版本和进程数上的工作量，
保留原短样本，不能只删除较快样本或只延长某个版本。
单进程和最大进程数都应重复采样，以同版串行中位数作为效率分母并报告波动范围。

该驱动按 0.1 秒采样 `mpirun` 及其子进程的累计 RSS，并每秒记录 PSS、进程交换页、
系统换入/换出及日志中已报告的工作量。`peakTreeRssKiB` 是同时累计的 RSS，包含共享页重复计数；
`peakTreePssKiB` 按比例分摊共享映射，未能读取时为 null。两者的峰值均覆盖完整进程生命周期，
包括初始化与输出，不能称为循环计算窗口内的内存峰值。日志可能缓冲，trace 的 `reportedWork`
只能帮助近似识别计算阶段。内存采样不参与求解状态或停止判断。
`peakSolverRssKiB` 是命令行包含二进制名称的进程所观察到的最大 RSS；
以脚本包装并 `exec` 另一个名称的程序时，这个名称匹配指标可能遗漏实际求解进程，不能用于内存收益判断。
此类诊断仍使用完整进程树 RSS/PSS，并单独记录包装器、库与实际二进制的哈希。

通用后端的 `preconditionerReuses` 计数使用旧矩阵版本预条件器的逻辑 RHS 求解，
`preconditionerRetries` 计数其中需要重建后重试的次数。`linearSolves` 计数逻辑 RHS，
`krylovIterations` 和 `trueResidualChecks` 计入两个尝试；重建时间计入 `preconditionerSeconds`。
已通过真实残差检查的初值不进入 Krylov。递推范数达标但真实范数未达标时可在原迭代预算内
继续至多两次；额外 Krylov 工作和真实残差检查仍计入原逻辑 RHS。最终验收目标保持不变。
Block Jacobi/ASM 的局部因子准备现在在显式 setup 阶段计时，避免计入首个 RHS 的 solve 时间。

## 通用划分和内存存储

MPI 读取网格时，rank 0 读入完整网格并计算一次全局所有权，然后逐 rank 分发 owned 单元及
stencil 所需 halo 的几何与拓扑。工作 rank 构造自己的局部 Mesh；构造前释放传输缓冲，
root 同时释放完整网格。读取与构造阶段的错误在进入下一阶段前传递给所有 rank。
这个路径保持全局 ID、物理 patch、face 方向和至少三层 ghost 的契约。

可选 METIS 以通信量为目标，同时约束单元数与单元面关联数；连通图请求连续分区，
随后恢复精确 owned 单元配额。选择依据及适用范围见
[Python 网格图模型](partition-model.md)和[实验报告](../reports/partition-model-2026-10-07.md)。
没有 METIS 时继续使用确定性的图生长实现。依赖和静态库链接说明见[构建说明](../building.md)。

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

## 当前工程规模证据

以下为 2026-10-07—08 的冻结参考 v5 与优化 v12 对照，不是旧 Eigen/PETSc 对照。
固定网格及工作量，1/12 rank 各三次取中位数，4/8 rank 各一次，合计 48 个正式测点。
计算时间为最慢 rank 的持续循环窗口；全部正式点至少 60 秒。SIMPLE 384 次外迭代计量
383 个间隔；PISO 三棱柱 10 步计量 9 个间隔，Voronoi 20 步计量 19 个间隔。

| 算例 | 单元数 | 12 rank 参考 → 优化计算秒 | 同 rank 版本加速 | 自身强扩展效率参考 → 优化 |
| --- | ---: | ---: | ---: | ---: |
| PISO 畸变三棱柱 | 1,024,000 | 426.788 → 113.601 | 3.757× | 39.51% → 43.02% |
| PISO Voronoi | 512,000 | 185.438 → 85.268 | 2.175× | 41.49% → 45.62% |
| SIMPLE 畸变三棱柱 | 500,094 | 353.050 → 109.925 | 3.212× | 30.37% → 31.73% |

12 rank 进程树生命周期峰值 PSS 中位数分别下降 32.93%、9.38%、35.58%。
PISO 14 组全部 owned 字段对照通过，完整回归通过。工程 SIMPLE 是固定工作量吞吐、
`maxIterations` 且退出码 2，不是稳态收敛；小型已收敛 SIMPLE 的字段一致性由回归覆盖。

并行效率随算例和 rank 变化，当前结果并未达到理想线性扩展；8 rank SIMPLE 效率没有
提升。单机测量也不能外推成多节点扩展。参考 v5 已含此前的若干框架优化，所以表格
不代表全部改动相对最初仓库的收益。
原始数字、重复范围、配置、硬件、审计和完整图见 [版本化报告](../reports/flow-optimization-2026-10-07.md)。

## 运行纪律与结果状态

各驱动持有排他锁并逐项启动 mpirun，不能阻止用户手动启动其他程序。测量期间避免
游戏、后台大负载、另一数值任务或编译，记录绑定及环境。最快点以 1–3 分钟计算窗口为
目标；固定工作量下较慢的串行/参考点自然更长，不能为了让每个点同长改变工作量。
跨版本、rank、重复样本核对相同输入哈希、停止条件、步数/迭代数和全部 rank 报告。

退出码 0 加 `converged` 才接纳完成运行；固定迭代吞吐只接纳退出码 2 加
`maxIterations` 且工作量完整。退出码 1、信号、超时或 `numericalFailure` 均无效，不能
因已存在 JSON 或先前时间步文件而判成功。`--resume` 应使用相同研究配置和哈希，修改
输入、二进制或分析口径用新的证据目录；旧短样本/失败样本仍保留并标出。

求解设置 A/B 可显式使用 `--variant-case variant:case=/absolute/alternative-input`，
驱动保存对应输入清单、哈希和 `settingsOverride` 标记。只对相同输入的二进制比较报告
软件加速；不同迭代安排报告设置收益，并单独检查全部场及守恒。比较不同时间步时，
每个案例应覆盖相同物理时间，但指定各自正确的 `--accepted-steps`，不能以相同步数代替。

当前 Hypre/KSP 内部通信、PCApply 和归约热点需 PETSc profiling 另测，不能用框架
计数器缺项证明它们没有发生。阶段计时有嵌套，不能简单加总成总时间。
