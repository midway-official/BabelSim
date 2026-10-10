# 08 运行与并行计算

[用户手册](README.md) · 上一章：[数值设置与收敛](numerics.md) · 下一章：[结果处理与继续计算](results-and-restart.md)

## 启动命令

```bash
build/babelsim-solve -case work/heat -time baseline
mpirun -np 4 build/babelsim-solve -case work/heat -time mpi4
```

求解器由 `case.bs` 中的 `solver` 选择。`-time` 是结果运行名；省略时使用 `output.bs` 的 `timeName`。正式求解程序没有 `--help`、`--version` 或 `--list-solvers` 选项；完整参数表见[命令行参考](../reference/command-line.md)。

同一次 MPI 运行使用相同的可执行程序与共享输入。分区、幽灵层和场交换由框架完成。`ghostLayers` 默认 3，当前支持的最小值也是 3。处理器数不应超过算例能够有效分配的单元规模；小算例用于检查功能，不用于评价并行加速。

## 日志与退出码

直接保存标准输出和错误输出，并保留进程退出码：

```bash
mkdir -p work/heat/logs
build/babelsim-solve -case work/heat -time logged \
  > work/heat/logs/logged.log 2>&1
run_status=$?
printf 'solver exit: %s\n' "$run_status"
```

| 退出码 | 含义 | 后续检查 |
| --- | --- | --- |
| `0` | 求解器返回完成状态 | 核对物理结束时刻、输出和关注量 |
| `1` | 输入、启动或应用错误 | 查看最后错误信息与路径、参数 |
| `2` | 未收敛或数值失败 | 查看失败时间步/迭代及各收敛指标 |

MPI 启动器可能将进程失败包装成额外错误信息或自己的状态；应结合日志中的 BabelSim 状态判断。程序的文件存在与输出目录名不能替代退出码。失败步通常不写成已接受结果，已有最后快照可能仍是更早的成功步。

## 流动日志的阅读

`mass` 表示离散质量指标，`dU` 表示速度变化，`rU`、`dP` 和 `rCoupling` 等描述方程或校正状态，`linU`、`linP` 描述线性求解状态。PISO 同时输出 Courant 数和校正次数。不同求解器的日志项不同，不能把线性残差当成全部物理耦合误差。

`converged=true` 和 `accepted=true` 只针对代码中配置的数值条件。检查是否达到 `endTime`、是否存在每个应保存的时刻、是否发生截断/模型限幅，并比较工程关注量。瞬态计算一般不会要求每一步的 `U` 与上一步完全不变。

## 性能记录

```bash
mpirun -np 4 build/babelsim-solve -case work/heat -time perf4 \
  -performance work/heat/logs/perf4
```

该目录保存每个进程的性能 JSON。统计包含构建、装配、求解或通信计数，具体字段以当前输出为准。用相同网格、物理参数、时间步和收敛阈值比较进程数，同时记录每组是否完成相同工作量。高级方法见[性能指南](../developer-guide/performance.md)。

## 长时间计算

在能维持会话的终端或作业系统运行，并使用独立运行名、输入快照和日志。集群调度系统的 MPI 启动方式应遵循本机部署规则。不要同时让两个任务写入同一个结果运行目录。

结果写出量近似与保存时刻数、输出场数和单元数成正比。通过 `writeInterval` 和 `writeFields` 控制存储，保证关注时刻有数据。仅保存稀疏时刻会限制频谱、周期和瞬态结构的分析；为了减少磁盘占用而增大间隔时，要记录这一限制。

## 并行一致性

```bash
python3 tools/compare_parallel_results.py \
  work/heat/results/baseline/0.05 work/heat/results/mpi4/0.05 \
  --atol 1e-8 --rtol 1e-8
```

工具按全局单元编号合并拥有者数据后比较。重复全局编号、缺失字段或几何不一致不是正常浮点差异。切换 METIS 或 MPI 进程数会改变分区与运算顺序，应采用合适的容差并检查误差分布；不能假定二进制结果完全相同。
