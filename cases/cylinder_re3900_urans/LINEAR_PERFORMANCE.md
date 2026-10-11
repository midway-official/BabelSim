# Re=3900 单时间步线性求解性能优化

2026-10-11 完成验收与交付。测量目录保留为 `validation/linear_speed_20261010`。
以当前 49,350 单元圆柱网格、dt=0.01、8 个绑定物理核的 MPI 进程为工程基准，
60 步计算窗口中位数从 **80.613434 s 降至 73.332042 s**，时间减少 **9.03%**，
加速比 **1.0993**。窗口按完成的 60 个物理时间步折算，由 **1.343557 s/步降至
1.222201 s/步**。所有样本均完成规定步数并通过原有外层、耦合、质量和湍流收敛检查。

硬件为 AMD RYZEN AI MAX+ 392，当前环境为 WSL，MPI 使用 8 个核；每进程
OMP/OpenBLAS/MKL 线程数为 1。框架使用 GCC 11.4、PETSc 3.25.5、本机构建优化，
包括 O3、march=native、LTO、fast-math，并保留 finite 检查、关闭 FMA contraction。
本结果是这台机器、该网格和配置上的测量，不代表多节点扩展或其他矩阵上的统一收益。

三组正式样本在一次预热后顺序交错运行。计时期间没有并行编译、测试或场后处理。

| 样本 | 原程序及配置 / s | 优化程序及配置 / s | 同组时间减少 |
| --- | ---: | ---: | ---: |
| 1 | 84.447164 | 74.850639 | 11.36% |
| 2 | 80.613434 | 73.332042 | 9.03% |
| 3 | 79.683566 | 71.775529 | 9.92% |
| 中位数 | 80.613434 | 73.332042 | 9.03% |

计时字段为 `max_rank(iterationComputeSeconds)`：从首次外迭代完成后开始，扣除
窗口中已计量的结果输出。它排除网格读取、生成、划分、初始场读取、首次矩阵模式
和首次外迭代等初始化工作；每个后续求解所需的矩阵数值更新、PC 重建、通信及
收敛检查均计入。只在末尾输出场，末尾输出发生在最后一次 iteration 计时之后。
**这里是持续计算窗口的每物理步摊销成本，并非排除整个首步后的逐步计时序列。**
窗口包含 60 步中除首个外迭代之外的工作，日志报告与监测开销也包含在内。

持续基准从正式运行 t=71.5 的完整保存场建立新的初场，关闭一次性潜势初始化。
所有对照使用完全相同的初场；首步为 Euler，之后为 BDF2。它恢复的是场值，
没有恢复 BDF2 历史，因此不是原正式运行的历史续算。模拟时间标签为 0 至 0.6。

本次默认配置只调整以下四项线性控制：

| solution.bs 配置键 | 原值 | 当前值 |
| --- | ---: | ---: |
| equation.momentum.relativeTolerance | 1e-12 | 1e-10 |
| equation.pressureCorrection.relativeTolerance | 1e-11 | 1e-8 |
| equation.pressureCorrection.absoluteTolerance | 1e-14 | 1e-13 |
| equation.pressureCorrection.preconditionerRebuildInterval | 8（默认） | 64 |

动量继续使用 BiCGStab/Block Jacobi，压力继续使用 CG/Hypre；二者保留热启动。
dt、网格、离散、校正次数、外层容差及 k/omega 输运设置沿用原值。
更长的 PC 重建间隔仅允许复用预条件器，矩阵 A 始终更新；迭代增长触发的提前
重建和失败后重建重试仍生效。各方程仍按当前 A 的原始真实残差验收：

```text
||b - A*x|| <= max(atol, rtol * max(||r0||, ||b||))
```

通用代码修改仅在 `src/backend/petsc_linear_system.cpp`：分布式 CG 使用 PETSc
single-reduction 变体合并两个内积归约；热启动通过 `VecNormBegin/End` 合并
`||b||` 与 `||b-A*x0||` 的归约。实现根据线性求解器类型和 communicator 工作，
不检查物理方程名、圆柱几何或算例规模。最终真实残差检查、原目标尺度、迭代预算
和残差重启逻辑继续使用。新增独立制造解测试覆盖两种初始残差尺度及 1e-12、1、
1e12 的 RHS 尺度。

**收益主要来自线性配置，不能把 9.03% 全部归因于通用代码。**仅更换代码而保持
原配置的 30 步初筛为 44.414 s，原基准为 44.421 s；这组单样本没有建立独立的
整步实现加速证据。其最大场差异为 8.10e-12，数值验收通过。

60 步完整回调的工作计数与阶段计时如下。阶段计数包含完整回调，计时窗口排除
首个外迭代，因此不能将阶段时间强行相加为上述窗口时间。

| 指标 | 原版 | 优化版 |
| --- | ---: | ---: |
| 接受物理步数 | 60 | 60 |
| 全部外迭代数 | 476 | 475 |
| 总 Krylov 迭代 | 39,002 | 28,622 |
| 压力 Krylov 迭代 | 19,468 | 10,841 |
| 压力 PC 建立次数 | 159 | 34 |
| 总 PC 失败重试次数 | 0 | 0 |
| 线性求解阶段最慢 rank 时间中位数 / s | 21.846 | 15.515 |
| PC 阶段最慢 rank 时间中位数 / s | 1.631 | 0.390 |

数值验收预算在筛选前固定：全域/尾流、普通/面积加权场相对 L2 差异不超过 1e-5，
力系数绝对差不超过 1e-4。按 global ID 合并 owned 字段，检查完整覆盖、有限值
和湍流场正值，比较 U、p、k、omega、mut 及独立计算的涡量。三组 60 步最终场
均通过，最大相对差为 **5.14949e-8**，力系数最大绝对差为 **3.83280e-8**。
优化版接受步的最大质量指标为 3.61487e-12、耦合指标为 1.38403e-11、湍流指标
为 1.16410e-7，分别满足 1e-8、1e-7、1e-6 的原有条件。

另从算例原始初场、保留潜势初始化完成了 30 步冷启动对照。两版全部收敛，最大
场差异 **6.18965e-9**，力系数差 **2.80278e-10**。计算窗口为原版 50.683 s、
优化版 45.758 s。这只有一组短程样本，主要用于启动正确性检查，不作为正式
重复性能结论。以上场与力的比较证明所测试窗口内的数值一致性，不建立长期湍流
统计、网格/时间步无关性或实验物理准确性。

筛选中拒绝了 k/omega 热启动、组合宽容差方案及 12-rank 候选，因为其 k 场
相对差约 1.20e-5 至 1.30e-5，超过相同预算。编译曾与部分配置初筛短暂重叠，
这些筛选计时没有纳入正式性能结论。误差门槛没有在观察结果后放宽。

验证已完成：`make petsc-check`（串行、2/4 rank）、`make test`、`make test-mpi`、
`make test-piso-coupling`，以及三组成熟场和一组冷启动的 owned 场验收。

完整证据保存在本算例的 `validation/linear_speed_20261010`，包括：

- `final-ab-audit.json`、`final-ab-validation-01/02/03.json`：正式统计和逐组场验收。
- `final-ab/records.json`、`final-ab/study.json`：原始报告、输入/程序哈希、机器、绑定与环境。
- `cold-ab-audit.json`、`cold-ab-validation-01.json`：原初场验证。
- `bin/baseline`、`bin/optimized`、`source/`、`binary-source.json`：冻结程序与源码。
- `original-case/numerics/solution.bs`：调整前的默认线性配置。
- `test.log`、`test-mpi.log`、`test-piso-coupling.log`、`petsc-check.log`：回归记录。

原始 owned 场、日志及输入副本保存在 `results/linear_speed_20261010`。已经冻结的
正式运行 `formal_dt0p01_optimized_v3` 保持原程序、原输入及暂停状态。当前默认
配置用于后续创建的新输入副本。新构建程序为 `build/babelsim-solve`。

复测使用冻结输入与程序，选择尚未使用的输出及 work-root 目录：

```bash
bench_root=cases/cylinder_re3900_urans/validation/linear_speed_20261010
python3 tools/benchmark_solver_optimization.py \
  --variant baseline="$bench_root/bin/baseline" \
  --variant optimized="$bench_root/bin/optimized" \
  --case mature="$bench_root/input-final-base" \
  --variant-case optimized:mature="$bench_root/input-final-optimized" \
  --accepted-steps mature=60 --ranks 8 --repeat 3 --warmup 1 \
  --minimum-compute 60 \
  --output "$bench_root/rerun" \
  --work-root cases/cylinder_re3900_urans/results/linear_speed_20261010/rerun
```

测量后顺序运行 `tools/validate_cylinder_optimization.py` 比较对应两版的
`results/benchmark-r8-n01`（以及 n02、n03）。使用不同的证据输出路径，保留
已有测量。性能测量期间不要并行执行这些后处理任务。
