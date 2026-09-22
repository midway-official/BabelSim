# Polyhedral FVM 并行一致性验证

验证提交：`24f4dcff43f3ff143726a31722ad1b44d5e09de5`  
分支：`codex/polyhedral-fvm`  
验证日期：2026-09-23

这次验证把“算子一致性”和“线性求解器一致性”分开。算子层要求同一全局单元在
1/2/4 rank 上得到相同的面值、梯度、通量、散度和离散贡献；求解器层要求各 rank
达到同一真实残差合同后，按 global cell id 对齐的最终场满足固定误差界限。没有把
不同 MPI 进程数的 Krylov 迭代轨迹或位模式相同作为验收条件，因为默认 IC/ILUT
是每个 rank 的局部因子，分区改变后其矩阵和填充顺序本来就会改变。

## 已执行的门禁

`make test-mpi` 通过。关键输出如下：

| 层次 | 证据 |
|---|---|
| 公开算子 | `parallel_math_test` 在 1/2/4 rank 均报告 `23 operations x 3 diffusion x 2 gradient methods, poisoned halos, maxError=0` |
| 数值合同 | `numerical_contract_test` 在 2/4 rank 均通过 30 个 solver/preconditioner/scale 组合 |
| halo/拓扑 | `parallel_domain_test` 通过 3 层 ghost、异步第一层交换、面 owner 发布者和分区面算子 |
| polyhedral 图分区 | `parallel_unstructured_test` 通过乱序 hex 图分区和 halo |
| MPI SIMPLE | `parallel_simple_test` 的 1/2/4 rank 外迭代均为 153，中心速度相同 |
| 热传导 | 1/2 rank 结果按 global id 对齐，最大绝对差 `4.1302393e-12` |

`make test-mpi-poiseuille` 也通过。native v3 网格的 SIMPLE 收敛外迭代为 1/2/4
rank=`1156/1151/1152`；按 global id 比较最终场得到：

| 比较 | `max |ΔU|` | `max |Δp|` |
|---|---:|---:|
| 1 → 2 rank | `6.5397226e-7` | `1.7120154e-6` |
| 1 → 4 rank | `7.6763065e-7` | `1.7498935e-6` |

两组都满足测试固定的 `5e-6 + 5e-6 max(|a|,|b|)` 逐分量界限，且 post 读取 v3
结果并生成 VTK/Tecplot 成功。

## 扭曲 polyhedral 短时矩阵

使用当前提交的 `tests/simple_parallel_consistency_test.py`，选取 384 单元三维 cube、
480 单元扭曲三维 polyhedral case，以及相同几何但 4 层 ghost 的 case；每个 case
分别运行 steady/Euler/BDF2，使用 1/2/4 rank。总计 27 次启动，瞬态每次比较 5 个
物理时间步，全部收敛并通过字段、global-id 覆盖、分片几何和物理时间检查。

最大差异（同时覆盖每个 case 的所有比较时间）：

| case | 算法 | `max |ΔU|` | `max |Δp|` |
|---|---|---:|---:|
| cube | steady | `1.3651e-11` | `3.2535e-12` |
| cube | Euler | `8.7204e-13` | `3.1641e-12` |
| cube | BDF2 | `6.9516e-13` | `2.4482e-12` |
| warped / warped4 | steady | `1.5665e-11` | `3.5098e-12` |
| warped / warped4 | Euler | `1.9105e-12` | `5.6657e-12` |
| warped / warped4 | BDF2 | `1.1108e-12` | `3.5952e-12` |

`warped` 与 `warped4` 的结果相同到本次比较精度，说明额外 ghost 层没有被错误地当作
物理数据参与结果或输出。

## 算子实现的差异来源

1. 公开算子在执行前使 cell/face 输入的 halo 有效；测试先故意污染 ghost，再调用
   `grad`、`interpolate`、`flux`、`div`、`laplacian`、`reconstruct` 及其组合。1/2/4
   rank 与串行参考的最大误差为 0，说明 owner/neighbour 面方向、第一层 halo 和三层
   非正交中间量同步没有引入可见的分区误差。
2. Green--Gauss 和最小二乘梯度都按局部 cell/face 拓扑计算；扭曲网格上的 corrected
   扩散、LinearUpwind 和组合通量也通过同一 global-id 结果比较。面访问仍由
   `owner` 和至多一个 `neighbour` 决定，因而没有六面体固定邻居的隐式分支。
3. `SparseAssembly` 只生成 owned 行。局部矩阵保留 owned-owned 项，跨分区面的系数
   由 Krylov halo 在 `A*x` 中按 ghost cell 追加；相同 row 的重复子面在固定面顺序中
   先求和。这保持了 CSR/LDU 热路径，不需要把变长面拓扑带入每次 SpMV。

## 求解器实现的差异来源

默认分布式路径使用 CSR SpMV：先发起第一层 Krylov halo，在本地 interior rows 上
计算，等待后计算 boundary rows 和 remote coupling。每次 SpMV 计数一个 halo exchange；
因此异步路径改变的是通信/计算重叠，不是离散算子的加法顺序。

线性求解器的全局内积和残差范数使用 `MPI_Allreduce(MPI_SUM)`。浮点加法不满足结合律，
不同 rank 数可能产生最后几位差异，但该因素不是本次迭代数变化的主因。为隔离因素，
对同一 32 单元 heat case 做了默认 ILUT、无预条件器和 AMG 三组 A/B：

| 预条件器 | rank | Krylov | SpMV | halo | global reductions | 中位 wall time (s) |
|---|---:|---:|---:|---:|---:|---:|
| ILUT | 1 | 5 | 15 | 0 | 0 | 0.474 |
| ILUT | 2 | 16 | 39 | 49 | 77 | 0.463 |
| ILUT | 4 | 24 | 55 | 65 | 110 | 0.483 |
| None | 1 | 82 | 172 | 0 | 0 | 0.438 |
| None | 2 | 82 | 172 | 182 | 339 | 0.474 |
| None | 4 | 82 | 172 | 182 | 339 | 0.497 |
| AMG | 1 | 5 | 15 | 0 | 0 | 0.450 |
| AMG | 2 | 5 | 35 | 45 | 46 | 0.476 |
| AMG | 4 | 5 | 35 | 45 | 46 | 0.481 |

无预条件器时 1/2/4 rank 的迭代数完全相同，最终 `T` 的最大差为
`1.67e-16`；因此 MPI 归约顺序虽会影响舍入位，但没有导致收敛轨迹分叉。默认 ILUT
的 `5/16/24` 变化来自每 rank 只对本地块做 ILUT：分区后远程耦合不在本地因子中，
局部图、行列顺序、drop/fill 结果都随 rank 数改变。AMG 的全局粗层使本例迭代数恢复为
`5/5/5`，代价是每次预条件应用多做粗层 SpMV 和全局归约。

扭曲 480 单元 steady case 的默认 ILUT 性能计数也显示同一趋势：

| rank | Krylov | SpMV | halo bytes | global reductions | 中位 wall time (s) |
|---:|---:|---:|---:|---:|---:|
| 1 | 21382 | 25164 | 0 | 3160 | 0.745 |
| 2 | 23393 | 27968 | 134,960,168 | 78,728 | 0.749 |
| 4 | 28204 | 33404 | 171,505,616 | 93,766 | 0.800 |

这是小算例，wall time 受进程启动和计时噪声影响，不能外推强扩展；但它足以说明
“rank 一致的最终场”和“rank 一致的迭代数”是两个不同目标。当前 1/2/4 运行仍在
固定外层收敛合同下完成，输出场满足扭曲 case 上述误差界限。

## 当前策略和后续使用建议

- 默认继续使用局部 ILUT/IC、CSR SpMV、第一层异步 halo 和现有 Allreduce。它保留了
  polyhedral 拓扑改造的性能目标，不为逐位复现引入每次迭代的全局排序或串行归约。
- 并行一致性验收比较 global-id 对齐后的最终场、真实残差、有限性、收敛状态和关键
  计数；不要求不同 rank 数的局部预条件器产生相同 Krylov 迭代轨迹。
- 需要更强 rank 一致性的研究/回归算例可选择 AMG。上面的 A/B 说明它能显著稳定迭代
  数和字段误差，但应在目标规模上重新测量粗层内存、归约比例和总时间后再设为默认。
- 只有在出现跨 rank 的最终场超出固定误差界限时，才优先检查面 owner/neighbour、
  halo 有效性、局部矩阵的 remote coupling 和真实残差；不要仅因迭代数不同就改变
  离散算子或关闭通信重叠。

复现命令：

```bash
make test-mpi
make test-mpi-poiseuille
python3 tests/simple_parallel_consistency_test.py \
  --cases cube warped warped4 --modes steady euler bdf2 \
  --output /tmp/babelsim-simple-mpi-20260923-representative
python3 tools/benchmark_backend.py --case cases/heat --ranks 1,2,4 \
  --repeat 3 --warmup 0 --output /tmp/polyhedral-consistency-perf-20260923 \
  --binary build/babelsim-solve
```

