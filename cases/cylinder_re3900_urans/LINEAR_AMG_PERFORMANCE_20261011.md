# Re=3900 时间步求解性能：AMG 预条件器比较

2026-10-11 已完成 30 步筛选、60 步一次预热与三组交错正式对照，以及对应场验收。
**最终选择压力 CG/Hypre BoomerAMG、动量和 k/omega BiCGStab/Block Jacobi 的混合
PC 配置。**该配置已位于 `numerics/solution.bs`，本轮无需再调整生产求解配置。
正式 60 步窗口中位数为 **77.344848 s**，全部 AMG 为
**102.683831 s**，全部 AMG 时间增加 **32.7610%**。
按 60 个接受物理步摊销为 **1.289081 s/步** 和
**1.711397 s/步**。结论限于本轮机器、矩阵、初场和所测候选。

本轮测试 AMG 用于动量和湍流输运系统是否能降低整步计算时间，不修改物理求解器
或后端代码。所有候选使用同一冻结程序，改变范围为 PC 类型、Krylov 类型、PC
重建间隔，以及单独列出的 k/omega 热启动设置。原有线性容差、外层接受条件、
离散方法和物理参数均保持本轮基线值。

## 工程基准与计时

网格为 49,350 个单元，dt=0.01，8 个绑定物理核的 MPI 进程；硬件为
AMD RYZEN AI MAX+ 392，环境为 WSL。MPI 使用 `--bind-to core --map-by core`，
OMP/OpenBLAS/MKL 每进程各为 1 线程。冻结程序为
`validation/amg_pc_20261011/bin/current`，各批 `study.json` 记录的 SHA256 均为：

```text
7d2214f3b8dbbf8987ef7726bc8466124ed7f17d05496b2008bba0c737613ef6
```

成熟场输入继承前一轮 `validation/linear_speed_20261011/input-pressure6`。
其来源是暂停正式运行 t=71.5 的完整保存场，关闭一次性潜势初始化；各候选使用
共同的初始场值。首步为 Euler，之后按配置使用 BDF2，没有恢复 BDF2 历史，
因此是成熟场上的新推进窗口，不是原正式运行的 BDF2 历史续算。
30 步窗口的时间标签为 0–0.3，正式 60 步窗口为 0–0.6。

计时主指标为 `max_rank(iterationComputeSeconds)`：在首个外迭代结束后开始，
扣除窗口内已计量的结果输出，排除网格生成/读取/划分、初场读取、首次矩阵模式
及首次外迭代等初始化。后续矩阵数值更新、PC 建立/复用、线性求解、通信和
收敛检查均计入，日志和监测开销仍在窗口内。
正式窗口除以 60 所得到的每步值是**完成物理步数的摊销成本，不是逐时间步的
延迟序列，也没有排除整个首个物理步**。阶段计数及阶段时间覆盖完整求解回调，
不能与排除首个外迭代的窗口时间强行加总。

最终采用的线性配置如下，所有候选的 absolute/relative tolerance 均沿用此表。

| 方程 | KSP / PC | absoluteTolerance | relativeTolerance | PC 间隔 | 热启动 |
| --- | --- | ---: | ---: | ---: | --- |
| 动量 | bcgs / bjacobi | 1e-15 | 1e-10 | 8（默认） | true |
| 压力修正 | cg / hypre | 1e-13 | 1e-6 | 64 | true |
| k 输运 | bcgs / bjacobi | 1e-14 | 1e-8 | 8（默认） | false |
| omega 输运 | bcgs / bjacobi | 1e-14 | 1e-8 | 8（默认） | false |

压力已使用 AMG；本轮“全部 AMG”指进一步把动量、k、omega 的 PC 改为 AMG。
PC 间隔同时限制矩阵年龄与非零迭代 RHS 工作量，**不是物理时间步数**；A 始终
按当前有效系数更新。Hypre 的现有粗化、插值及按 Krylov 类型选择的策略，
GAMG 的插值平滑含义与默认 Chebyshev 条件，见
[后端说明](validation/amg_pc_20261011/backend-notes.md)。

## 30 步候选筛选

三批筛选均为单次运行、无预热。第一批和第二批分别包含新测的本批基线，
基线窗口为 39.310524 s 与 41.508502 s，说明批次间存在运行状态波动。
因此各批单样本用于筛选方向，不把跨批绝对时间之差当作优化收益。第三批没有
再次运行基线，热/冷初猜候选只作为本批成本观察；其场验收使用第二批相同初场
和控制的基线输出。正式结果由后续交错重复对照建立。

下表最大场差为 U、p、k、omega、mut 和重建涡量在全域/尾流、普通/面积加权
相对 L2 指标中的最大值。时间为完整 30 步的计算窗口。

| 批次 / 候选 | 动量、k、omega 的改动 | 窗口 / s | 外迭代 | 最大场差 | 场验收 |
| --- | --- | ---: | ---: | ---: | --- |
| screen / baseline | 本轮混合 PC 基线 | 39.310524 | 252 | 参考 | — |
| screen / momentum_amg | 仅动量 PC=hypre | 42.708702 | 255 | 1.66476e-7 | 通过 |
| screen / transport_amg | 仅 k/omega PC=hypre | 50.641851 | 252 | 9.93750e-6 | 通过 |
| screen / all_amg | 三者 PC=hypre，间隔 8 | 49.573512 | 255 | 9.93292e-6 | 通过 |
| screen / all_amg_lag64 | 三者 PC=hypre，间隔 64 | 50.163226 | 255 | 1.00228e-5 | 不通过 |
| screen2 / baseline | 本轮混合 PC 基线 | 41.508502 | 252 | 参考 | — |
| screen2 / bjacobi_lag32 | 三者维持 bjacobi，间隔 32 | 42.148016 | 252 | 1.98347e-7 | 通过 |
| screen2 / gmres_hypre_lag32 | 三者 gmres/hypre，间隔 32 | 46.336739 | 255 | 1.38588e-5 | 不通过 |
| screen2 / bcgs_gamg_lag32 | 三者 bcgs/gamg，间隔 32，插值平滑 1 | 未完成 | 0 | 无最终场 | 首步失败 |
| screen3 / bjacobi_warm_lag32 | 三者 bjacobi/间隔 32；k/omega 热启动 | 38.687161 | 252 | 1.20252e-5 | 不通过 |
| screen3 / hypre_warm_lag32 | 三者 hypre/间隔 32；k/omega 热启动 | 43.131024 | 255 | 1.20186e-5 | 不通过 |
| screen3 / hypre_cold_lag32 | 三者 hypre/间隔 32；k/omega 保持冷初猜 | 54.761124 | 255 | 9.94446e-6 | 通过 |

除 GAMG 候选外，各候选均完成 30 个接受物理步，求解器状态为 `converged`。
场验收不通过的候选仍被拒绝；收敛不能替代与基线的数值一致性检查。
表中 all_amg、transport_amg 和 hypre_cold_lag32 的最大差接近 `1e-5` 门槛，
30 步通过不能预先证明 60 步仍通过。

GAMG 候选在首个物理步内返回 `numericalFailure`，没有接受物理步。已保存的
omega 方程报告为 392 次 Krylov 迭代，最终真实残差约 15.12957，框架状态为
`numericalFailure`。这份记录没有提供具体的负 KSP reason，不据此确定某个
PETSc breakdown 原因；该候选不能计入有效 30 步性能样本。

数据出处为 `validation/amg_pc_20261011/screen{,2,3}/records.json`、各批
`study.json` 和目录根下对应的 `screen*-validation.json`。候选配置分别保存在
`screen-plan.json`、`screen2-plan.json`、`screen3-plan.json` 及各冻结输入中。

## 迭代减少与求解成本

第一批同批对照显示，AMG 的线性求解阶段成本抵消了迭代数减少。
以下为完整求解回调的各方程计数和最慢 rank 阶段时间，不能直接加总为窗口时间。

| 方程 | 基线迭代 | all_amg 迭代 | 基线线性求解 / s | all_amg 线性求解 / s |
| --- | ---: | ---: | ---: | ---: |
| 动量 | 1,195 | 621 | 0.626 | 2.237 |
| 压力修正 | 2,885 | 2,896 | 2.781 | 2.656 |
| k 输运 | 3,864 | 1,594 | 1.709 | 4.991 |
| omega 输运 | 4,560 | 2,292 | 1.931 | 6.694 |
| 全部线性系统 | 12,504 | 7,403 | 7.033 | 16.577 |

PC 阶段最慢 rank 时间从 0.257 s 增至 2.103 s。把全部 AMG 的重建间隔增至
64 后，PC 建立次数由 278 降至 53、PC 时间降至 0.499 s，但线性求解阶段仍为
17.733 s，窗口为 50.163 s，未在该单样本中获得时间优势。
这些样本没有 PC 失败重试，不能把增加的时间解释为重试造成。

第三批热启动对照也有相同现象：Block Jacobi 与 Hypre 分别为 5,309 和 4,552
次 Krylov 迭代，但线性求解阶段分别为 4.586 s 和 9.517 s；两者都因场差超限
未被接受。因此本轮不以“AMG 迭代更少”或“热启动更快”作为交付依据。

## 精度预算与正式 60 步对照

筛选前固定的预算为：全域/尾流、普通/面积加权场相对 L2 差异不超过 `1e-5`，
壁面力系数估计绝对差不超过 `1e-4`。现有校验按 global ID 合并 owned 字段，
检查覆盖、有限值、湍流场非负及 omega 正值，并独立重建涡量。30 步全部已完成
场比较的候选均满足力系数预算，最大绝对差为 1.80626e-7；被拒绝者主要因 k 场
相对差超限，不能在看到结果后放宽预算。

每个线性系统继续使用当前 A、原始 b 的未预条件真实残差合同：

```text
||b-A*x|| <= max(atol, rtol * max(||r0||, ||b||))
```

配置热启动时，r0 来自当前字段初猜，RHS 范数仍保留在尺度中；r0 小于 RHS
范数时与冷初猜尺度相同，r0 更大时按原合同取较大的初始残差尺度。
外层、耦合、质量及湍流接受条件没有因更换 PC 而改变。

正式对照采用本轮混合 PC 基线与 `all_amg`（间隔 8、k/omega 冷初猜），同一冻结
程序和共同初场，60 步、一次预热、三组按固定顺序交错运行。
该对照用于确认短程筛选方向及更长窗口的精度，正式测量不与编译、测试、场后处理
或其他求解作业并行。

| 正式样本 | 混合 PC 基线 / s | 全部 AMG / s | 全部 AMG 时间增加 |
| --- | ---: | ---: | ---: |
| 1 | 76.571390 | 105.723481 | 38.0718% |
| 2 | 77.344848 | 101.811169 | 31.6328% |
| 3 | 77.991372 | 102.683831 | 31.6605% |
| 中位数 | 77.344848 | 102.683831 | 中位数比值：32.7610% |
| 窗口除以 60 / s 每物理步 | 1.289081 | 1.711397 | 同上 |

六个正式样本均完成 60 个接受且收敛的物理时间步，所有 rank 的状态为 converged，
窗口均超过 60 秒；最短正式窗口为 76.571390 s。三组计数各自一致，基线每组
476 次外迭代、5,678 次线性求解，全部 AMG 为 479 次外迭代、5,702 次线性求解。
两路均无 PC 失败重试。减少 Krylov 的收益不足以抵消更高的线性阶段成本：

| 完整回调指标 | 混合 PC 基线 | 全部 AMG |
| --- | ---: | ---: |
| Krylov 总迭代 | 23,259 | 13,781 |
| PC 建立次数 | 515 | 511 |
| 最慢 rank PC 时间中位数 / s | 0.488653 | 4.389854 |
| 最慢 rank线性求解时间中位数 / s | 14.125130 | 34.058381 |

全部 AMG 的三组最终场相对差均为 **1.016762104e-05**，最大位置为尾流 k 的
普通相对 L2，略超过 `1e-5` 预算；壁面力估计最大绝对差为 **1.730691994e-07**，
满足 `1e-4` 预算。owned 场覆盖完整、有限值、k/mut 非负、omega 正值均已验证。
因此全部 AMG 虽满足求解器收敛/物理步接受条件，仍未通过本轮场一致性验收；
不因超限幅度较小而修改预算。

最终选定基线的接受物理步最大质量指标为 `1.20001e-10`、耦合指标为
`3.56172e-11`、速度变化为 `9.91784e-8`、湍流变化/残差为
`2.39447e-7`/`1.14272e-7`，符合原接受条件。当前生产文件与冻结基线的全部
`equation.*` 控制相同；生产文件保留原有一次性潜势初始化，成熟场基准关闭该初始化。
这些数值推进与场一致性证据不证明长期湍流统计、网格/时间步无关性或实验物理准确性。

`audit_amg.py` 已执行成功：输入/程序哈希、计划控制差异、接受步数、外迭代区间、
计时、全部数值比较均完整自洽，状态为 `complete_with_rejections`，不存在遗漏或
矛盾。该状态表示拒绝证据也已检查，不能解释为所有候选通过精度验收。
审计比较已记录的 SHA；只重新计算小型配置字典的哈希，没有重新计算二进制、
大型网格和初场文件的哈希。输入副本和程序的实际哈希由基准驱动器冻结及核对。
最终选择和当前配置校验记录位于 `selection.json`。

## 证据与复测

本轮证据根目录为 `validation/amg_pc_20261011`；原始 owned 场和日志位于
`results/amg_pc_20261011`。主要文件如下：

- `bin/current`：各候选共用的冻结程序。
- `screen{,2,3}/records.json`、`study.json`：30 步原始阶段计时、工作量、输入与程序哈希。
- `screen*-validation.json`：各已完成候选的 owned 场和壁面力估计校验。
- `backend-notes.md`：通用后端实现和适用边界。
- `final-all_amg-validation-01/02/03.json`：三组正式场/力比较。
- `selection.json`、`amg-audit.json`：最终选择及完整证据审计。
- `input-final-baseline`、`input-final-all_amg`：正式 60 步冻结输入。
- `final-ab/records.json`、`study.json`：正式对照证据，另有 `summary.json` 的完整正式统计。
- `audit_amg.py`：本轮审计脚本，已运行并生成 `amg-audit.json`。

复测时选用尚未使用的输出及 work-root 目录，保留已有测量：

```bash
bench_root=cases/cylinder_re3900_urans/validation/amg_pc_20261011
python3 tools/benchmark_solver_optimization.py \
  --variant baseline="$bench_root/bin/current" \
  --variant all_amg="$bench_root/bin/current" \
  --case mature="$bench_root/input-final-baseline" \
  --variant-case all_amg:mature="$bench_root/input-final-all_amg" \
  --accepted-steps mature=60 --ranks 8 --repeat 3 --warmup 1 \
  --minimum-compute 60 \
  --output "$bench_root/rerun-amg" \
  --work-root cases/cylinder_re3900_urans/results/amg_pc_20261011/rerun-amg
```

计时结束后再使用 `tools/validate_cylinder_optimization.py` 比较对应两路相同
`benchmark-r8-n01/n02/n03` 标签的最终场，独立保存新的统计和验收证据。
