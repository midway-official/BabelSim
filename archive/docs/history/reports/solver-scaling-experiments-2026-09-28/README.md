# 求解阶段强扩展实验与复现

> 历史开发记录：保留原报告的接口、命令、数据与结论，只用于追溯对应时期的工作。
> 不作为当前实现或性能的使用说明；现行文档见 [文档中心](../../../../../docs/README.md)，归档范围见 [历史索引](../../README.md)。

本目录保存 2026-09-28 在当前 CPU/PETSc 后端上进行的实验。性能指标唯一采用每次运行所有 rank 中最大的 `solverComputeSeconds`；预热不计入中位数。不使用进程外部墙钟、网格读取、Case 初始化或输出时间计算收益。`run.json` 中仍保留外部墙钟等诊断信息。

主分析见 [强扩展优化报告](../strong-scaling-analysis-2026-09-28.md)。主程序没有因本实验而修改；实验输入使用独立目录，正式 NACA 配置不变。

## 证据结构

- `screen/`：单次 4-rank 配置筛选。三个 PISO 试验误将 `maxIterations` 设置为 5（它实际是单步 Picard 外迭代上限），因此未达到外迭代收敛要求。这些失败不用于有效加速比或推荐配置；输入副本与失败日志保留。
- `flow-screen/`：修正为 `maxIterations=1`，五步、dt=0.01，最少 8 个 PISO corrector、最多 50 个、耦合容限 1e-5；比较压力 Hypre/GAMG 与动量 ASM。全部条件以该目录输入副本为准。
- `scaling101/`：101306 单元热方程的 Hypre/BJacobi/GAMG1，以及同网格五步 PISO 的 Hypre/GAMG1，1/2/4/8 rank。每点一次预热、三次正式测量，候选与 rank 顺序按固定随机种子交错。
- `scaling500/`：500094 单元热方程，Hypre 与 BJacobi，1/4/8 rank，每点一次预热、三次正式测量。
- `naca/`：当前 NACA 正式输入的独立五步副本，仅改变压力 PC；4 rank，一次预热、三次正式测量。它是配置 A/B，不具有单 rank 基线，不能据此计算 NACA 的强扩展效率。
- `petsc-profile/`：单独启用 PETSc 事件采样的诊断运行，不混入无 profiling 的性能曲线。PETSc 的 PCApply 含调用 Hypre 的成本，但 PETSc 消息/归约统计不能保证完整覆盖 Hypre 内部 MPI。

正式 study 下有主机、二进制 SHA-256、输入文件哈希、MPI 绑定与全部命令；每个 `ranks-P/repeat-NN` 下有 `stdout.log`、每 rank 性能 JSON 与 `run.json`。`summary.csv` 用候选自身的 1-rank 时间计算加速比和效率。没有单 rank 基线时对应列留空。

生成输入的小文件随证据保存；大网格由脚本再生并以哈希核验。CSV 场结果位于 `study.json` 中记录的 `/tmp/babelsim-solver-scaling-20260928/` 独立算例目录，未将大体积重复输出放入版本库。每次运行的输出哈希保存在 `run.json`，正式场对比结果另存 `validation.json`。

## 可执行复现

在仓库根目录执行。`CASE_ROOT` 和 `EVIDENCE_ROOT` 应指向新的目录；准备脚本拒绝覆盖已有实验。

```bash
make -j4 all
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PETSC_OPTIONS
STUDY_ROOT=$(mktemp -d /tmp/babelsim-scaling-reproduce-XXXXXX)
CASE_ROOT="$STUDY_ROOT/cases"
EVIDENCE_ROOT="$STUDY_ROOT/evidence"
python3 tools/prepare_solver_scaling_cases.py --output "$CASE_ROOT" --include-naca

python3 tools/run_solver_scaling_study.py \
  --case "heat101-hypre=$CASE_ROOT/heat101-hypre" \
  --case "heat101-bjacobi=$CASE_ROOT/heat101-bjacobi" \
  --case "heat101-gamg1=$CASE_ROOT/heat101-gamg1" \
  --case "piso101-hypre=$CASE_ROOT/piso101-hypre" \
  --case "piso101-gamg1=$CASE_ROOT/piso101-gamg1" \
  --ranks 1,2,4,8 --repeat 3 --warmup 1 \
  --output "$EVIDENCE_ROOT/scaling101"

python3 tools/run_solver_scaling_study.py \
  --case "heat500-hypre=$CASE_ROOT/heat500-hypre" \
  --case "heat500-bjacobi=$CASE_ROOT/heat500-bjacobi" \
  --ranks 1,4,8 --repeat 3 --warmup 1 \
  --output "$EVIDENCE_ROOT/scaling500"

python3 tools/validate_solver_scaling_study.py "$EVIDENCE_ROOT/scaling101"
python3 tools/validate_solver_scaling_study.py "$EVIDENCE_ROOT/scaling500"

python3 tools/run_solver_scaling_study.py \
  --case "naca-hypre=$CASE_ROOT/naca-hypre" \
  --case "naca-gamg1=$CASE_ROOT/naca-gamg1" \
  --ranks 4 --repeat 3 --warmup 1 \
  --output "$EVIDENCE_ROOT/naca"

python3 tools/validate_solver_scaling_study.py "$EVIDENCE_ROOT/naca" --reference-rank 4
```

网格和输入哈希应与原证据匹配（路径元数据不要求相同）。不要并发运行不同性能实验；场对比也放在计时实验结束后执行。短程对比不代替完整物理时间、网格与时间步精度验收。

## 无需修改代码的 PETSc profiling

本次实测确认，虽然主程序使用 `PetscInitializeNoArguments()`，仍能通过 `PETSC_OPTIONS` 环境变量启用初始化阶段读取的 `-log_view`。这不意味着 `-ksp_type` 或 `-pc_type` 能覆盖程序内的显式 KSP/PC 配置。

```bash
PETSC_OPTIONS="-log_view :$EVIDENCE_ROOT/petsc-events.log" \
mpirun --bind-to core --map-by core -np 4 \
  build-petsc/babelsim-solve \
  -case "$CASE_ROOT/heat101-hypre" -time profile \
  -performance "$EVIDENCE_ROOT/profile-ranks"
```

`PCApply`、`KSPSolve`、各层 MatMult 等事件可能嵌套，不能直接相加。外部 Hypre 的内部阶段需额外 profiling。正式 timing 仍关闭 profiler，以同一计时边界比较。
