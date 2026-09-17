# BabelSim 后端性能基准

性能实验由显式的 application 选项和 `tools/benchmark_backend.py` 驱动。运行时和
Physics 不打印性能信息，也不判断物理收敛；`-performance <目录>` 只在调用者明确启用时
让应用写出每个 rank 的结构化 JSON。普通求解命令的输出行为保持不变。

## 固定顶盖驱动流案例

`tests/performance/cavity-100k.json` 冻结了 320×320×1、102400 单元、Re=1000、
`rho=1`、`mu=0.001`、无湍流模型以及离散方法。先在仓库外生成案例，避免把大网格和
结果提交到 Git：

```bash
python3 cases/cavity/validation/generate_cavity_case.py \
  /tmp/babelsim-cavity-100k --cells 320 --re 1000 --cluster 0 \
  --convection linearUpwind --gradient greenGauss \
  --velocity-relaxation 0.3 --pressure-relaxation 0.3 \
  --max-iterations 30000
python3 tools/benchmark_backend.py \
  --case /tmp/babelsim-cavity-100k \
  --ranks 1,2,4,6,8,12 --warmup 1 --repeat 3 --timeout 1800 \
  --mode complete --output /tmp/babelsim-evidence/cavity-100k
```

驱动先为每个 rank 执行一次不计入统计的 warmup，然后再执行 `--repeat` 个正式样本。
物理核范围内使用 `mpirun --bind-to core --map-by core`；超过可用物理核但不超过可见
逻辑 CPU 时，使用独立的 `--bind-to hwthread --map-by hwthread` SMT 策略。每次运行的
启动命令、绑定策略和可见 CPU 拓扑都写入 `run.json`/`metadata.json`，并通过
`--report-bindings` 把实际 rank 绑定写入原始日志。超过可见逻辑 CPU 的 rank 会在启动
前拒绝，驱动不会主动添加 `--oversubscribe`。

驱动为每次执行建立独立目录，保存 `stdout.log`、`run.json`、rank 性能 JSON、Git
dirty diff、二进制 SHA-256、编译器/MPI/CPU/线程环境、Makefile 哈希和生效的默认优化
开关信息，以及 `results.csv`。`--resume` 会
跳过已有终态记录（warmup 也有独立记录）；成功返回但缺少性能 JSON 会直接报错。超时会清理整个 MPI 子进程组，
并且只记录为 `timeout`，不会进入收敛排名。

每个 `run.json` 还记录 benchmark 进程树的并发 RSS 峰值（`peakTreeRssKiB`）和观测到的
solver 进程最大 RSS（`peakSolverRssKiB`）。采样器读取 Linux `/proc`，只用于性能证据，
不参与求解器停止、收敛或错误判断；进程退出瞬间的短暂 `/proc` 缺失可能使峰值略保守，
因此内存数据按“观测到的峰值”解释。不同时间运行的各 rank 峰值不能相加为同时峰值。

汇总中的 `phaseTimesByRank`/`localCountersByRank` 先按请求的 `-np` 分组，再按实际
本地 rank 分组；`criticalPathByRank` 对每次运行取所有 rank 的阶段最大值，然后只对这些
关键路径样本求 min/mean/max。这样不会把不同并行度的 rank 0 混在一起，也不会把分别取到
的 rank 最大值相加成虚假的总时间。

rank JSON 还记录 `mpiInitSeconds`（每个进程的 MPI 初始化，一次性启动开支）、
`caseSetupSeconds`（Case、网格读取、分区和配置建立，一次性开支）、`solverSeconds`
（Physics solver 调用）、`solverComputeSeconds`（从 solver 调用中扣除已计量结果写出的
时间）和 `applicationSeconds`（MPI 初始化完成后的应用阶段总时间）。
`outputWrites`/`outputSeconds` 是结果写出的次数和累计 I/O 时间；写出可能发生在多个
时间步，属于 I/O 阶段，不应和每次迭代计算混在一起。`solverSeconds` 包含 solver 内
触发的结果写出，若要比较纯计算，应从同一 JSON 的 `outputSeconds` 及写出策略中说明
这部分成本，不能简单把阶段时间相加。汇总文件的 `phaseTimesByRank` 按 rank 保留
Case/solver/solverCompute/application 阶段的 min/mean/max；MPI 启动和输出计数在 rank
JSON 中单独保留。

每个 rank 的性能 JSON 还包含只读的 `partition` 指标：global/local/owned/ghost cells、
local/owned faces、processor communication faces 和实际邻居 rank 数。汇总文件的
`partitionByRank` 按 MPI size 与 local rank 保留这些值的 min/mean/max，方便检查分区负载、
ghost 比例和通信面增长；它们只用于解释性能，不参与求解或收敛判断。

## 三类测量

当前脚本的 `complete` 模式使用生产求解器的真实停止条件，只有 `converged` 样本可以
用于完整 time-to-solution 排名。`throughput` 只改变证据标签，不能伪造固定迭代或改变
收敛条件；如果需要严格的固定外迭代窗口，应在独立 benchmark application 中实现，不能
修改生产求解器的算法。所有正式样本（包括 `maxIterations` 的 throughput 样本）都在
`wallClockByRank` 中按状态统计；`convergedWallClock` 只保留可用于完整求解排名的样本。

后端计数器包括线性求解、Krylov 迭代、SpMV、halo、全局归约、方程装配、预条件器 setup
和 apply 以及对应时间。`linearSolveSeconds` 包含其内部 SpMV/通信，阶段时间不能简单
相加为关键路径；并行总时间以应用墙钟和各 rank 的最大阶段时间解释。
顶层计数器是各 rank 的最大值，用于表示关键路径上的工作量/字节量；`local` 对象保留
当前 rank 的原始计数。若需要总通信量，应在汇总工具中对 `local.haloBytes` 求和，不能
把顶层 `haloBytes` 当成全局总字节数。

## 当前已实现的后端优化

`HaloExchange` 为每个交换计划缓存组件数缩放后的 MPI count/offset 数组，避免每次 SpMV
重新分配四个临时向量；并提供 `beginFirstLayer`/`finishFirstLayer`，让分布式 Krylov
先发起 `MPI_Ialltoallv`，计算本地 interior SpMV，再等待并计算远程耦合。Krylov 临时
布局不再每轮清零：owned 值全部覆盖，跨分区耦合所需的第一层 ghost 由交换完整覆盖。
分布式 SpMV 默认使用后端私有的行式 CSR 视图；`CSR_SPMV=0` 可切换回 Eigen 列式
稀疏乘，`ASYNC_HALO=0` 可切换回阻塞 halo，两个开关都只影响后端。CSR 视图不替代
Eigen 矩阵，预条件器和 AMG 仍使用原有实现。这些机制位于 `src/parallel` 和
`src/algebra`，公共 Physics DSL 不可见。
串行 `PreparedLinearSolver` 也复用后端 CSR SpMV；首次 `compute` 建立 pattern，后续
`factorize` 只复制连续系数并检查 pattern，避免每个外迭代重新分配稀疏结构。该视图同样
是内部实现，Eigen 仍保留给 IC、ILUT、AMG 和其他因子化操作。串行 A/B 可用
`SERIAL_CSR_SPMV=0` 恢复 Eigen SpMV。
CSR 热循环使用连续数组指针和 GCC/Clang 的最多 8 次循环展开提示；这是后端编译优化，
不改变 Physics DSL、稀疏模式或每行累加顺序。该优化在 102400 单元 pilot 中的 A/B 数值
见 [`backend-audit-2026-09.md`](/home/midway/BabelSim/docs/performance/backend-audit-2026-09.md)。

优化验收必须同时满足：

1. `make -j4 test`、所有 MPI/工作流回归和 `git diff --check` 通过；
2. 与同一网格、同一配置的基线全局 cell ID 对齐，场相对 L2 ≤ 1e-5、L∞ ≤ 1e-4，
   并通过现有更严格数值合同和 Ghia 检查；
3. 至少三次独立重复，报告中位数、范围和变异系数；
4. 只在同一配置下比较代码 A/B，分别报告配置调优、单进程代码收益和 MPI 加速；
5. 没有测量证据时不声称达到 15% 端到端加速或某个 MPI 效率目标。

## 计时和复现审计

编译参数由 Makefile 的 `OPTFLAGS` 冻结并写入实验元数据。`-O3`、本机指令集、LTO 和
自动向量化属于后端实现选择；Physics 不应包含 SIMD intrinsic、CSR、MPI 或原始存储
指针。任何 SoA/对齐存储改动都必须保持 `Field` 的值语义、global cell ID、边界方向、
owned/ghost 同步合同，并在标量、矢量、张量和串行/MPI 路径分别做 A/B。
