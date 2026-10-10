# 性能测量与后端观察

[开发指南](README.md) · [运行手册](../user-guide/running-and-parallel.md) · [图分区](partitioning.md)

本章说明如何测量当前软件，而不提供跨机器性能承诺。固定工作量、数值成功条件、输入、程序、进程数和计时边界应一起记录。已有版本的数字保存在可选的 archive/docs 中，不纳入产品操作说明。

## 获取每个进程的计数器

```bash
python3 docs/examples/create_case.py --template piso --output work/performanceCase
mpirun -np 2 build/babelsim-solve -case work/performanceCase -time observation \
  -performance work/performanceCase/performance/observation
```

输出目录中每个 rank 有一份 JSON，应用日志也包含汇总。上述小算例用于理解记录格式，不能用 5 步粗网格测点判断工程规模并行加速。

| 计数/时间 | 解释 |
| --- | --- |
| matrixPatternBuilds | 稀疏模式建立次数 |
| matrixValueUpdates | 当前矩阵数值更新次数 |
| rhsOnlySolves | 矩阵不变、只换 RHS 的逻辑求解 |
| linearSolves、krylovIterations | 逻辑 RHS 数与实际 Krylov 工作 |
| preconditionerSetups、preconditionerReuses | PC 准备与旧版本 PC 复用 |
| preconditionerRetries | 重建旧 PC 后重试次数 |
| trueResidualChecks | 当前系统真实残差检查 |
| outputWrites、outputSeconds | 已计量结果写出 |
| haloBytes、haloExchanges | BabelSim 场同步，不包括 PETSc/Hypre 全部通信 |
| globalReductions | BabelSim 显式归约，不包括 KSP/Hypre 内部全部归约 |

并行阶段时间取最慢进程的关键路径，不直接相加所有 rank 的时间。计数器覆盖其实现范围，不能将 haloBytes 当作整个程序网络流量。

## 完整运行与持续循环窗口

墙钟时间包括初始化、首次矩阵/PC 和 I/O。solverComputeSeconds 从回调时间扣除已计量写出，仍包含首次迭代。iterationComputeSeconds 是通用 iteration 报告观测的持续循环窗口：排除首次完整迭代及此前初始化，并扣除窗口内已计量写出；reportedIterationIntervals 记录有效间隔数。

至少有两次适当的迭代观测才存在持续窗口。heat 等仅用 record 的报告路径不应推断具有相同 iteration 窗口。多外迭代 PISO 的工作单位还包括流动外迭代，不能将首个外迭代排除解释为整个首物理步排除。

阶段计数可能覆盖完整回调，不能强求它们与排除首步的循环窗口完全闭合。报告必须说明选择的时间字段与工作单位。

## 顺序测量相同输入

普通后端测量工具：

```bash
python3 tools/benchmark_backend.py \
  --case work/performanceCase --binary build/babelsim-solve \
  --ranks 1,2,4 --warmup 1 --repeat 3 --mode complete \
  --output work/performanceEvidence
```

工具保留输入/程序信息与重复样本。小型教程仍仅用于检查流程；正式研究应统一增加所有版本和进程数的工作量，保证最快点也有足够的有效计算时长。不要只延长慢版本或删除较快的短样本。

对不同二进制进行 A/B 比较使用 tools/benchmark_solver_optimization.py，详细参数可通过该 Python 工具的 --help 查询。每个 variant 使用独立工作目录，输入与停止条件保持一致，先做预热，再重复采样。采用持久证据目录并保存实际 owned 场，不能只留下时间 JSON。

## 强扩展与设置比较

同版本、同网格、同工作量下，`S(P)=T(1)/T(P)`，`E(P)=S(P)/P`。各版本分别使用自己的串行基准；同进程数版本加速另行报告。比较不同时间步时保持相同物理时间并记录不同步数；不能用相同步数代表相同物理工作。

更换 KSP/PC、容差、外层校正或湍流安排后，先比较解与接受条件，再描述设置收益。软件实现优化与数值设置变化分开报告。固定迭代 SIMPLE 吞吐可以记录迭代耗尽，但必须注明退出码 2 与未达稳态条件；它不代表收敛时间。

## 内存与进程树

记录完整生命周期峰值 RSS/PSS，包括初始化、运行库与输出。进程树 RSS 同时累加共享页可能重复计数，PSS 按比例分摊；无法读取的指标记录为空，不能当作 0。窗口计算时间与全生命周期内存峰值是不同口径。

包装脚本应记录实际二进制与运行库；仅按程序名字匹配进程可能遗漏被 exec 的求解进程。观察完整进程树并说明采样间隔、交换页和后台负载。

## PETSc 生命周期的解释

长期 Equation 使用长期 PETSc Mat/Vec/KSP/PC，分离模式、矩阵更新、RHS 更新和重复 solve。PC 可以复用旧矩阵版本，A 始终更新到当前装配；失败重建/重试仍检查当前 A 的真实残差。参数定义见[数值参考](../reference/numerical-controls.md)。

linearSolves 计逻辑 RHS；两个尝试的 Krylov 工作与真实残差检查均计入实际计数。已经满足真实残差的初值可以不进入 Krylov。不要用“少一次 solve 调用”直接推断完成了更少数学工作。

## 报告验收

保留源码和二进制哈希、输入清单、编译优化、MPI/线程与核绑定、硬件、工作量、状态、重复次数、分布和关注量比较。超时、信号、配置错误或 numericalFailure 的点不能因存在 JSON 就接纳为成功测点。物理准确、并行一致和性能是分别验证的结论。
