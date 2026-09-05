# Re=1000 顶盖驱动流性能案例

这是一个可复现的性能案例，不替换 `cases/cavity` 的回归设置。网格为 `64×64×1`，
采用双曲正切壁面加密（`cluster=1.5`）、`linearUpwind` 对流、最小二乘梯度和中心型
扩散，流体参数为 `rho=1`、`mu=0.001`，因此以顶盖速度和腔体长度归一化的 Reynolds
数为 1000。两套线性方案都使用相同的 SIMPLE 松弛、物理外迭代上限和收敛判据。

当前 `numerics/solution.bs` 保留第二套 GMRES+AMG 设置。复现实验时只替换该文件最后两
行，然后为每个进程数使用不同的 `-time` 标签：

```text
# 方案 A：已有基线
vectorSolver bicgstab ilut 1e-12 1e-7 1000
scalarSolver cg incompleteCholesky 1e-12 1e-7 1000

# 方案 B：重启 GMRES + AMG
vectorSolver gmres amg 1e-12 1e-7 1000 gmresRestart=30 amgMaxLevels=12 amgCoarseSize=48 amgSmoothingSteps=2
scalarSolver gmres amg 1e-12 1e-7 1000 gmresRestart=30 amgMaxLevels=12 amgCoarseSize=48 amgSmoothingSteps=2
```

运行示例：

```bash
/usr/bin/time -f 'wall_seconds=%e\npeak_kib=%M\nexit_status=%x' -o /tmp/re1000-np4.time mpirun -np 4 build/babelsim-solve -case cases/cavity/benchmark-re1000 -time perf-np4
```

同一方案应分别使用 `-np 1`、`-np 2`、`-np 4`，并使用不同的结果标签。`results/` 是
运行产物，不属于输入案例；不要将它作为源代码提交。
