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

AMG 也做了独立的 160²、50 外迭代 A/B。标量压力方程将平滑步数从 2 调到 1 时，solver
中位数从约 6.03 s 降到约 5.79 s（约 4.5%），但 Krylov 迭代从约 15.7k 增到约 22.1k；
刷新间隔 1/2/4 的差异低于重复样本噪声。默认配置暂不修改，以保持现有数值和配置语义；
这项结果仅作为 AMG 成本--收敛速度的调参证据。

## 尚未验证的风险

- OpenMPI 的非阻塞集体是否在当前环境真正推进，需要用等待时间、SpMV 时间和 profiler
  分离验证，不能仅凭 `MPI_Ialltoallv` 名称判断重叠。
- AMG 的全局粗层复制、refresh 频率和分区质量尚未完成系统 A/B；串行 AMG 与分布式
  AMG 不能合并成同一条排名曲线。
- `Vec3/Tensor3` 仍为 AoS。若后续采用 SoA/AoSoA 或对齐分配，必须保留 public Field
  值语义、global ID 对齐、边界方向和 poisoned-halo 测试，并分别验证标量/矢量/张量。
- 当前 pilot 是 `maxIterations`，不属于完整收敛结果；超时和数值失败不计入最快配置排名。
- 行式 CSR 目前仅覆盖双精度标量 Krylov SpMV（串行与分布式）；预条件器仍使用 Eigen，
  尚未证明在所有网格形状、非结构网格和高阶稀疏模式上都优于 Eigen。默认路径必须继续
  保留 `CSR_SPMV=0`/`SERIAL_CSR_SPMV=0` 回退和数值 A/B 检查。
