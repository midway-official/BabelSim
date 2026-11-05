# Re=3900 时间步求解性能：2026-10-11 进一步优化

以包含上一轮优化的 HEAD `622f6d9016c3ac52f025f59bb11a8397c0435918` 为新基线，
49,350 单元、dt=0.01、8 个绑定物理核的 60 步持续窗口中位数从
**69.588227 s 降至 65.079547 s**，时间减少 **6.4791%**，加速比 **1.06928**。
按完成的物理时间步折算，从 **1.159804 s/步降至 1.084659 s/步**。
这是本轮相同机器与输入的实测结果；上一轮报告的 9.03% 不计入本轮收益。

最终交付只调整 `numerics/solution.bs` 中一项线性控制：

```text
equation.pressureCorrection.relativeTolerance 1e-6
```

原值为 `1e-8`。压力继续使用 CG/Hypre、absoluteTolerance=1e-13、PC 重建间隔 64
和热启动。动量及 k/omega 的线性控制、网格、dt、离散方法、校正次数及外层接受条件
沿用新基线。物理求解器与通用后端源码均保持该基线版本。

## 计时与对照

硬件为 AMD RYZEN AI MAX+ 392，环境为 WSL，8 MPI ranks 使用
`--bind-to core --map-by core`；OMP/OpenBLAS/MKL 每进程均为 1 线程。
冻结程序与运行库采用本机构建，GCC 11.4、PETSc 3.25.5、O3、march=native、LTO、
fast-math，保留 finite 检查并关闭 FMA contraction。完整机器、编译、绑定、输入与
程序哈希保存在 `validation/linear_speed_20261011/final-ab/study.json`。

三组程序/配置各进行一次预热与三次正式运行，按固定随机顺序交错执行。
计时期间没有编译、测试、场后处理或其他求解作业；各路使用同一进程树内存监测器。
所有正式样本均完成 60 个物理时间步、通过求解器接受条件，且有效窗口超过 60 秒。

| 正式样本 | 当前基线 / s | 仅配置调整 / s | 配置与增量代码 / s | 仅配置的配对时间减少 |
| --- | ---: | ---: | ---: | ---: |
| 1 | 67.260991 | 64.000881 | 63.817893 | 4.8470% |
| 2 | 69.588227 | 67.651603 | 66.222461 | 2.7830% |
| 3 | 69.595783 | 65.079547 | 66.805430 | 6.4892% |
| 中位数 | 69.588227 | 65.079547 | 66.222461 | 中位数比值：6.4791% |

主指标为 `max_rank(iterationComputeSeconds)`。该窗口在首个外迭代结束后开始，
扣除已计量的窗口内结果写出，排除网格生成/读取/划分、初场读取、首次矩阵模式与
首次外迭代。后续每次矩阵数值更新、PC 建立/复用、线性求解、通信和收敛检查均计入。
每步值为窗口除以 60，**不是排除整个首步后的逐步延迟序列**。
日志与监测开销仍在窗口内。阶段计数覆盖完整求解回调，不能与窗口时间强行加总。

基准从暂停正式运行 t=71.5 的完整保存场建立共同的新初场，关闭一次性潜势初始化。
只恢复场值，首步采用 Euler，随后 BDF2；没有恢复 BDF2 历史，模拟时间标签为 0–0.6。
因此该基准是成熟流场上的新时间推进窗口。

## 工作量与数值验收

| 完整回调指标 | 当前基线 | 最终配置 |
| --- | ---: | ---: |
| 接受物理步数 | 60 | 60 |
| 全部外迭代 | 475 | 476 |
| 逻辑线性求解数 | 5,670 | 5,678 |
| 总 Krylov 迭代 | 28,622 | 23,259 |
| 压力 Krylov 迭代 | 10,841 | 5,460 |
| 压力 PC 建立次数 | 34 | 36 |
| 最慢 rank 线性求解时间中位数 / s | 14.598 | 10.732 |
| 最慢 rank PC 阶段时间中位数 / s | 0.361 | 0.376 |

线性工作下降是主要收益来源；PC 次数没有下降。外层多执行了一次迭代，
因此不能把收益描述为减少物理时间步或放松外层接受条件。
每个系统仍按当前矩阵 A 的原始真实残差验收：

```text
||b - A*x|| <= max(atol, rtol * max(||r0||, ||b||))
```

筛选前固定的误差预算为：全域/尾流、普通/面积加权场相对 L2 差异不超过 `1e-5`，
壁面力系数估计绝对差不超过 `1e-4`。按 global ID 合并 owned 字段，独立验证完整覆盖、
有限值、湍流非负及 omega 正值，并比较 U、p、k、omega、mut 与重建涡量。

三组最终配置均通过。相对本轮基线最大场差异 **2.34834e-7**，壁面力系数最大绝对差
**1.20370e-7**。另核对初场、网格、离散方法与控制区间的哈希，并与上一轮更严格线性
容差的原始参考场比较，三组最大场差异 **2.46085e-7**、力系数差 **8.20419e-8**，
避免多轮配置调整累计越过同一个误差预算。

对完成且 converged 的物理步，最终配置最大质量指标为 `1.20001e-10`、耦合指标
`3.56172e-11`、湍流变化 `2.39447e-7`、湍流残差 `1.14272e-7`，满足原有
`1e-8`、`1e-7`、`1e-6`、`1e-6` 条件。中间外迭代日志没有用作物理步计数。

原始初场、保留潜势初始化的 30 步检查也已完成。两版全部接受且收敛，最大场差异
**2.53209e-7**、力系数差 **1.52658e-7**，通过相同预算。短程窗口为原配置 46.100 s、
新配置 43.000 s；它只有一组样本，主要用于启动正确性检查，不作为正式重复性能结论。
结果写入 `cold-validation.json`。
上述比较仅支持所测窗口内的数值一致性，不证明长期湍流统计或网格/时间步无关性。

## 通用代码实验与拒绝的设置

代码实验限定于 PETSc 线性系统：利用已计算的 `r=b-A*x` 求增量 `A*delta=r`，
再加回当前 x，消除热启动重复矩阵乘法。所有实际 KSP 调用从零开始，真实残差重启和
PC 失败重试保持原始 RHS、尺度与验收目标。实验通过四类 Krylov、多 RHS 尺度、
大初猜消去、真实残差重启、warm/cold 切换、PC 重试及框架/MPI/PISO 回归。

相同设置的 5 步 PETSc profile 中，MatMult 从 6,674 次降至 6,424 次，VecNorm 从
4,757 次降至 4,607 次。但需要两个额外 owned 工作向量、复制与加回；30 步单独代码
筛选没有显示整步收益。正式重复测量中，代码与配置的中位数为 66.222 s，高于仅配置的
65.080 s；相对于配置的配对变化方向也不一致。因此本轮不采用该代码实验。
源码快照、可逆 patch、制造解测试与冻结程序完整保存在证据目录，生产源码已恢复。

其他筛选结论：

- 12 ranks 的 30 步为 30.735 s，但尾流 k 相对差 `1.30210e-5` 超限；k/omega 容差
  收紧到 `1e-10` 后仍为 `1.03369e-5`，均拒绝。
- k/omega 改 GMRES 没有时间优势，最大场差 `2.17832e-5` 超限，拒绝。
- pressure rtol=`1e-5` 通过误差预算，但外迭代从 251 增至 301，30 步时间升至
  37.474 s；说明减少 Krylov 迭代不足以保证整步加速。
- pressure rtol=`1e-7` 和 pressure CG/GAMG 通过场验收，但单样本没有显示时间优势；
  4/6 ranks 较慢，没有进行额外场验收。

## 交付与复测

最终构建的 `build/babelsim-solve` 与冻结基线程序 SHA256 相同：
`7d2214f3b8dbbf8987ef7726bc8466124ed7f17d05496b2008bba0c737613ef6`。
配置变更应用于后续创建的新输入副本。已经冻结并暂停的
`formal_dt0p01_optimized_v3` 输入、程序及进程状态保持原样。

生产源码恢复后重新完成构建、`make petsc-check`、`make test`、`make test-mpi`、
`make test-piso-coupling`。日志以 `final-` 为前缀保存在本轮证据目录。

原始数据为 `validation/linear_speed_20261011/final-ab/{records,study,summary}.json`，
`final-ab-audit.json` 检查三路输入、程序哈希、接受步数、计时、重复样本与数值验收；
`audit_final.py` 可以重新生成该审计。各组场比较保存在
`final-settings-validation-01/02/03.json` 与 `strict-reference-validation-01/02/03.json`。
原始 owned 场位于 `results/linear_speed_20261011`。

选择新的输出与工作目录复测最终两路对照：

```bash
bench_root=cases/cylinder_re3900_urans/validation/linear_speed_20261011
python3 tools/benchmark_solver_optimization.py \
  --variant current="$bench_root/bin/current" \
  --variant settings="$bench_root/bin/current" \
  --case mature="$bench_root/input-final-current" \
  --variant-case settings:mature="$bench_root/input-final-optimized" \
  --accepted-steps mature=60 --ranks 8 --repeat 3 --warmup 1 \
  --minimum-compute 60 \
  --output "$bench_root/rerun" \
  --work-root cases/cylinder_re3900_urans/results/linear_speed_20261011/rerun
```

计时结束后再运行 `tools/validate_cylinder_optimization.py` 比较各组相同标签的最终场。
已保存的三路审计脚本针对本次三路完整证据；新的两路复测应独立保留其统计与验收。
