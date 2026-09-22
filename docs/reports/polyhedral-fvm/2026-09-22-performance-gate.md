# 2026-09-22 线性求解性能门禁

本次门禁把同一 `heat` 算例的旧版 `c6c9dfa` 和当前 `codex/polyhedral-fvm` 实现交错执行 7 组；每组依次测量 MPI 1、2、4 ranks。两边分别使用各自 checkout 中的同一算例版本（旧版 v2、当前版迁移后的 v3），release 编译选项、MPI 绑定、线性求解配置、输出开关和收敛停止准则保持一致。每次运行没有额外 warmup，因为每个进程由 benchmark driver 独立启动；求解器本身必须报告 converged 才计入统计。

复跑命令（先准备旧版 worktree 并分别构建两个 binary）：

```bash
git worktree add --detach /tmp/babelsim-baseline-20260922 c6c9dfa
make -j2 all                                      # 在旧版 worktree
make -j2 all                                      # 在当前工作区
python3 tools/benchmark_backend.py --case cases/heat \
  --ranks 1,2,4 --repeat 7 --warmup 0 \
  --output /tmp/polyhedral-current --binary build/babelsim-solve
```

交错运行的完整记录已压缩为 [2026-09-22-performance-gate.csv](2026-09-22-performance-gate.csv)；原始每次 JSON/日志位于本机 `/tmp/polyhedral-perf-gate-20260922/`。中位 wall time 和求解器计数如下：

| ranks | c6c9dfa 中位数 (s) | 当前 v3 中位数 (s) | 当前/旧版 | Krylov / SpMV / halo（两边相同） |
|---:|---:|---:|---:|---:|
| 1 | 0.462633 | 0.451746 | 0.976 | 5 / 15 / 0 |
| 2 | 0.456323 | 0.460435 | 1.009 | 16 / 39 / 49 |
| 4 | 0.464634 | 0.479247 | 1.031 | 24 / 55 / 65 |

该小型算例没有显示稳定的求解迭代或通信退化；1/2 ranks 在测量噪声内，4 ranks 中位数约高 3.1%，达到计划中的 2% 回归告警线附近。它不能证明目标生产规模的性能零损失，尤其不能替代缓存外大矩阵、固定分区的代数 replay。源码已把面顶点/拓扑访问限制在 setup 和输出层；若目标硬件上的大规模重复测量仍稳定超过告警线，应继续定位 setup、输出或 MPI rank imbalance 的来源。

独立的冻结 `A/b/x0` 重放、重复子面矩阵汇总、dense oracle 和 update 覆盖测试在 `tests/assembly_solver_test.cpp` 中执行；该测试使用固定 2x2 非对称矩阵和 BiCGSTAB 真实残差，不改变生产求解器容差。
