# 计算后端静态审计（2026-09-17）

本文只记录当前源码可直接验证的事实、待验证的性能假设和已经测到的结果。历史报告
不作为本次基线。基线提交为 `23997abd589f0dbca7034b2c52ddf2f0f098264f`，本次
工作区差异由基准驱动写入 `metadata.json`。

## 实际调用图

```text
Physics main.cpp
  -> public Case / Field / math / equ / diagnostics
  -> detail::FvmExecution
  -> ComputeBackend
  -> SparseAssembly + Eigen solver
  -> (distributed) DistributedLinearSolver + HaloExchange + MPI
```

Physics 看不到 `MPI_Comm`、CSR、Eigen 或 `MeshStorage` 原始数组。`RunTime`/`Case` 只
提供 DSL 执行桥和只读计数快照，应用通过 `-performance` 选择导出 JSON。

## 源码事实

| 热点 | 当前事实 | 层次 |
| --- | --- | --- |
| `Field<T>` | [`include/babelsim/field.h`](/home/midway/BabelSim/include/babelsim/field.h:57) 的值、边界迹和通量均为连续 `std::vector`；`Vec3/Tensor3` 为 AoS 值类型。 | public DSL |
| 网格几何/拓扑 | [`include/babelsim/mesh.h`](/home/midway/BabelSim/include/babelsim/mesh.h:94) 将几何、面、邻接、owned/ghost 和 global ID 分开存储。 | public geometry / private storage |
| 稀疏结构 | [`src/backend/eigen_assembly.cpp`](/home/midway/BabelSim/src/backend/eigen_assembly.cpp:35) 建立一次结构并更新系数；分布式后端另分 interior/boundary SpMV 结构。 | backend |
| halo 计划 | [`src/parallel/parallel_context.cpp`](/home/midway/BabelSim/src/parallel/parallel_context.cpp:546) 在构造时解析 global ID 请求；值交换只应复用该计划。 | parallel |
| Krylov SpMV | [`src/algebra/distributed_solver.cpp`](/home/midway/BabelSim/src/algebra/distributed_solver.cpp:511) 先 halo，再 interior SpMV，再处理边界和远程耦合。 | algebra/backend |
| 预条件器 | 同文件的 `computePreconditioner` 区分 IC、ILUT 和分布式 AMG；AMG 粗层仍有全局聚合和归约。 | backend |
| 观测 | [`include/babelsim/solver_control.h`](/home/midway/BabelSim/include/babelsim/solver_control.h:31) 保存求解、SpMV、halo、归约、装配、预条件器和结果写出计数及时间；application JSON 另列 MPI/Case/solver 阶段。 | generic runtime data |

## 已实现的改变

1. `HaloExchange` 为每个计划缓存常用组件数（标量、向量、张量）的 MPI count/offset
   数组。原先 `exchange` 每次构造四个临时 `std::vector<int>`；现在只在该组件布局首次
   使用时建立，后续 SpMV 不再分配这些数组。
2. `beginFirstLayer`/`finishFirstLayer` 使用 `MPI_Ialltoallv` 和持久的 send/receive
   buffer。`KrylovHalo` 在通信进行时计算 interior 矩阵乘，finish 后才计算远程项。
   `BABELSIM_ASYNC_KRYLOV_HALO=0` 可编译出同一算法的阻塞 A/B 路径。
3. `KrylovHalo` 不再每轮清零整个局部 cell 缓冲区。owned 值全部覆盖，跨分区耦合所需
   的第一层 ghost 由计划完整接收；其余更深 ghost 从不进入 Krylov 远程项。
4. 应用新增显式 `-performance <directory>`。它写 rank JSON，包含状态、墙钟、线性/SpMV/
   halo/归约/装配/预条件器计数和阶段时间；没有该选项时运行时不打印性能数据。
5. 分布式 Krylov 的 interior/boundary SpMV 增加私有行式 CSR 视图。Eigen 稀疏矩阵仍由
   预条件器和 AMG 使用；CSR 只缓存行偏移、列号和 Eigen 系数位置，方程更新时覆盖数值。
   `CSR_SPMV=0` 可恢复原 Eigen SpMV，便于同一配置 A/B。
6. 串行 PreparedLinearSolver 同样缓存后端 CSR 视图；稀疏矩阵首次建立 pattern，后续
   factorize 只复制值并检查结构，避免每个外迭代的稀疏对象赋值和临时分配；
   `SERIAL_CSR_SPMV=0` 保留 Eigen A/B 回退。
7. Application 级性能 JSON 增加 MPI 初始化、Case 构造、solver 调用和应用总时长，并在
   Case 内累计结果写出次数/时间；benchmark 汇总将计算、启动和 I/O 作为独立阶段，
   不与后端包含式计时相加。
8. benchmark 驱动在每次 MPI 运行期间采样 `/proc` 进程树，记录并发 RSS 峰值和 solver
   进程最大 RSS，并把实际使用的关键环境变量写入 `run.json`；这只提供观测数据，不
   参与求解停止或收敛判断。
9. 私有 CSR SpMV 的热循环缓存了原始数组指针，并在 GCC/Clang 下请求最多 8 次循环
   展开。该提示只作用于后端双精度行乘，保持每行累加顺序和矩阵模式；编译器不支持该
   提示时会退化为同一标量循环。
10. `Mesh::partitionInfo()` 为通用观测接口返回 global/local/owned/ghost cell、面和
    processor-neighbor 统计；性能 JSON 和 benchmark 的 `partitionByRank` 按 MPI size
    与 local rank 保存这些数据。该接口不参与方程组、场同步或收敛决策。
11. IC/ILUT 的因子化仍由 Eigen 完成；`src/algebra/inplace_preconditioner.h` 只在
    backend 内继承 Eigen 因子类型并复用输出向量及置换 scratch，保持原有缩放、置换和
    三角代入顺序。`INPLACE_PRECONDITIONER=0` 选择 Eigen 原始 `solve()`，便于逐配置回退。
12. `src/backend/algebraic_multigrid.cpp` 的每个 AMG level 预分配 `product`，将
    `right_hand_side - A*x` 改成复用向量的原地更新；分布式 AMG 的 fine-level residual
    同样复用已有 buffer。V-cycle 两条入口路径都会完整覆盖输出，因此 `apply` 不再
    预先清零输出向量。层级、聚合、平滑次数、Galerkin 矩阵和粗层求解器均未改变。

## 性能假设和实测结果

假设是 count/offset 分配和无效清零属于每次 Krylov 迭代的固定内存/分配成本，非阻塞
交换可以与 interior SpMV 重叠。使用生成的 320×320×1（102400 cells）Re=1000 算例，
将 `maxIterations` 临时设为 1 只做 pilot，保持同一网格、离散和线性配置：

| 路径 | ranks | 重复 | 墙钟中位数 | solver 中位数 | SpMV 中位数 | 结论 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| async + 缓存布局 | 2 | 3 | 4.861 s | 0.921 s | 0.221 s | pilot，未收敛 |
| blocking + 缓存布局 | 2 | 3 | 4.859 s | 0.936 s | 0.223 s | pilot，未收敛 |
| async + Eigen SpMV | 2 | 5 | 4.953 s | 0.962 s | 0.235 s | pilot，未收敛 |
| async + 行式 CSR SpMV | 2 | 5 | 4.939 s | 0.926 s | 0.200 s | pilot，未收敛 |

该样本的端到端差异约为噪声范围，不能声称达到 15% 加速；它只证明新路径在 100k
网格上运行且计数器完整。相同 5 次样本的 solver 阶段约下降 3.7%，SpMV 阶段约下降
15%，但墙钟被网格加载、MPI 启动和输出成本主导，不能把它表述为端到端同幅度收益。
固定系统回放和完整生产收敛实验仍需按
[`tests/performance/cavity-100k.json`](/home/midway/BabelSim/tests/performance/cavity-100k.json)
执行。

串行 CSR 索引改为 32 位连续数组后，在同一个 160²、50 外迭代、BiCGSTAB+ILUT/CG+IC
配置上重新进行了 3 次 warmup + 3 次正式样本的 A/B。`SERIAL_CSR_SPMV=1` 与
`SERIAL_CSR_SPMV=0` 只改变 Krylov 的 SpMV 实现，Krylov 迭代次数保持一致：

| ranks | 路径 | 墙钟中位数 | solver 中位数 | SpMV 中位数 |
| ---: | --- | ---: | ---: | ---: |
| 1 | 行式 CSR | 11.55 s | 9.00 s | 1.27 s |
| 1 | Eigen 回退 | 11.85 s | 9.32 s | 1.60 s |
| 2 | 行式 CSR | 8.60 s | 5.85 s | 1.42 s |
| 2 | Eigen 回退 | 8.61 s | 5.95 s | 1.44 s |

该结果支持串行 CSR 作为默认实现，但不能推断所有稀疏模式都获得同样收益；Eigen
回退仍用于回归和新网格形状的 A/B。当前热点计数显示预条件器 apply 通常约占 solver
时间的 55--60%，SpMV 约占 20--25%，因此下一阶段应优先针对 IC/ILUT/AMG 的数据布局和
通信归约做剖析，而不是仅优化装配。

随后对同一台机器上的 102400 单元 pilot（固定 `maxIterations=1`、同一 Re=1000 cavity
和线性配置）做了“指针缓存 + 循环展开”局部 A/B。两组都经过一次 warmup 和 3 次正式
样本；A 组只缓存指针，B 组再启用 8 次展开。结果如下，均为后端 rank 0 计时中位数：

| ranks | 路径 | 墙钟中位数 | solver 中位数 | SpMV 中位数 | Krylov 迭代 |
| ---: | --- | ---: | ---: | ---: | ---: |
| 1 | 指针缓存 | 5.903 s | 2.523 s | 0.421 s | 949 |
| 1 | 指针缓存 + 展开 | 4.894 s | 1.814 s | 0.251 s | 949 |
| 2 | 指针缓存 | 6.118 s | 1.398 s | 0.326 s | 978 |
| 2 | 指针缓存 + 展开 | 5.007 s | 0.987 s | 0.232 s | 978 |

这是固定工作量的 throughput pilot，所有样本状态为 `maxIterations`，因此不能解释为
收敛时间或物理精度提升。相同配置下 solver 阶段约减少 28--29%，SpMV 阶段约减少
26--40%；墙钟仍包含 MPI/Case/输出等一次性开支。`make -j4 test` 和 `make test-mpi`
在该改动后通过，现有数值合同未发现迭代数或场结果变化。

AMG 也做了独立的 160²、50 外迭代 A/B。标量压力方程将平滑步数从 2 调到 1 时，solver
中位数从约 6.03 s 降到约 5.79 s（约 4.5%），但 Krylov 迭代从约 15.7k 增到约 22.1k；
刷新间隔 1/2/4 的差异低于重复样本噪声。默认配置暂不修改，以保持现有数值和配置语义；
这项结果仅作为 AMG 成本--收敛速度的调参证据。

最终代码还在同一台 12 物理核/24 逻辑 CPU 主机上完成了 102400 单元、固定 1 外迭代的
1/2/4/6/8/12 ranks pilot，每组 5 次正式样本。以 `criticalPathByRank.solverSeconds`
的均值计算（只表示计算阶段，不包含 Case/网格一次性开支）：

| ranks | solver 均值 | 相对 1 rank 加速 | 效率 |
| ---: | ---: | ---: | ---: |
| 1 | 3.190 s | 1.00 | 100.0% |
| 2 | 1.784 s | 1.79 | 89.4% |
| 4 | 1.153 s | 2.77 | 69.2% |
| 6 | 1.014 s | 3.15 | 52.4% |
| 8 | 0.843 s | 3.78 | 47.3% |
| 12 | 0.895 s | 3.57 | 29.7% |

这是固定工作量的 `throughput` pilot，所有样本状态为 `maxIterations`，不能当作完整
物理解收敛时间。1/2/4/6/12 组的墙钟变异系数仍受主机频率和 MPI 启动噪声影响（约
7.8%--18.9%），因此只作为当前硬件上的阶段性扩展性证据；没有据此宣称达到 4 rank
70% 或 8 rank 50% 的工程目标。solver 之外的 Case/网格/分区开支约 2.6--2.9 s，
说明端到端优化必须继续降低一次性并行初始化和分区成本，或单独报告 amortized cost。
该 pilot 观测到的单 rank solver RSS 约为 201 MiB（4--12 ranks）到 322 MiB（1 rank）；
进程树并发 RSS 会随 rank 数增加，因为当前 Case/网格数据在各 rank 保留一份，不能把
不同 rank、不同时间的峰值相加成同时内存占用。

为检查大问题规模下的存储和通信拐点，随后运行了 499849 与 1000000 单元的固定一次
SIMPLE 外迭代 pilot。50 万单元的 solver 关键路径为 33.54/16.62/9.13 s（1/2/4 ranks），
总墙钟为 40.30/28.04/18.52 s；100 万单元为 85.13/43.14/24.28/16.93 s
（1/2/4/8 ranks），总墙钟为 97.14/66.18/43.73/42.70 s。100 万单元 8 ranks 的
进程树峰值 RSS 约 12.90 GiB，已经接近当前主机的可用内存，故未启动 12/24 ranks。
这些样本都返回 `maxIterations`，只用于规模、内存和阶段成本分析，不能用于完整收敛排名。

在同一机器上又对 IC/ILUT apply 做了后端 A/B。500k 单元、`maxIterations=1`、相同
2003 次 Krylov 迭代的 1-rank 三次正式样本中，Eigen 原始路径的 solver 均值为 44.58 s、
预条件器 apply 均值为 36.98 s；原地路径分别为 15.40 s 和 8.01 s。原始路径墙钟中位数
为 47.58 s（35.93--69.86 s），原地路径为 21.80 s（21.74--22.08 s），因此基线噪声
较大，不能只用这个 pilot 推出固定百分比的端到端收益。2-rank 单次样本的 solver 为
23.67/10.63 s，apply 为 18.04/4.96 s（Eigen/原地）。

为了检查是否改变数值轨迹，完整 4096 单元 cavity 的两条路径都在 SIMPLE 第 2756 次
外迭代收敛，线性求解次数、Krylov 迭代和预条件器调用数相同；按 global cell ID 比较
`U`、`p` 的最大绝对差均为 0。solver 阶段从 33.24 s 降到 28.77 s，预条件器 apply
从 16.26 s 降到 11.63 s。该证据支持原地路径作为默认后端实现，但仍保留编译期开关和
Eigen 回退。

AMG 工作区 A/B 使用约 50k 单元的固定一次外迭代 pilot。1-rank solver 为 0.611/0.589 s，
预条件器 apply 为 0.171/0.159 s；2-rank solver 为 0.288/0.295 s，apply 为 0.138/0.145 s
（原 AMG/工作区复用），差异与运行噪声同量级。因此当前只记录“减少临时向量”的实现事实，
不宣称 AMG 已获得稳定速度提升；随后去掉无效输出清零的 65536 单元复测中，1-rank
solver 均值为 0.599 s、2-rank 为 0.283 s，仍不足以单独归因出稳定收益。粗层构造和
全局粗矩阵归约仍是后续独立热点。

## 尚未验证的风险

- OpenMPI 的非阻塞集体是否在当前环境真正推进，需要用等待时间、SpMV 时间和 profiler
  分离验证，不能仅凭 `MPI_Ialltoallv` 名称判断重叠。
- AMG 的全局粗层复制、refresh 频率和分区质量尚未完成系统 A/B；串行 AMG 与分布式
  AMG 不能合并成同一条排名曲线。
- `Vec3/Tensor3` 仍为 AoS。若后续采用 SoA/AoSoA 或对齐分配，必须保留 public Field
  值语义、global ID 对齐、边界方向和 poisoned-halo 测试，并分别验证标量/矢量/张量。
- 当前 pilot 是 `maxIterations`，不属于完整收敛结果；超时和数值失败不计入最快配置排名。
- 行式 CSR 目前仅覆盖双精度标量 Krylov SpMV（串行与分布式）；预条件器仍使用 Eigen，
  尚未证明在所有网格形状、非结构网格和高阶稀疏模式上都优于 Eigen。IC/ILUT 的原地
  wrapper 依赖 Eigen 因子类型的受保护成员，升级 Eigen 时必须重新检查 ABI 和操作顺序；
  默认路径继续保留 `INPLACE_PRECONDITIONER=0`、`CSR_SPMV=0`/`SERIAL_CSR_SPMV=0` 回退
  和数值 A/B 检查。
