# Re=1000、128×128 腔体吞吐基准

这个 Case 是 MPI/HPC 吞吐基准，不是 Ghia 精度验证算例。

- 网格：`128×128×1`，双曲正切壁面加密（`cluster=1.5`）；
- 物理：不可压层流，`rho=1`、`mu=0.001`，即 `Re=1000`；
- 离散：线性中心面插值、最小二乘梯度、二阶线性迎风对流、正交扩散；
- 线性后端：GMRES(30) 加局部聚合 AMG 预条件器；
- 终止条件：质量残差 `1e-8`、速度相对变化 `8e-4`。

最后两项刻意面向稳定、可重复的吞吐比较，使 1、2、4 个 MPI rank 都在十分钟内完成近似相同的外迭代工作量。它们不适合替代收紧 `continuityTolerance`、`velocityTolerance` 与线性容差后的精度/网格无关性验证。

在项目根目录编译后，分别运行：

```bash
build/babelsim-solve -case cases/cavity/benchmark-re1000-n128 -time throughput-np1
mpirun -np 2 build/babelsim-solve -case cases/cavity/benchmark-re1000-n128 -time throughput-np2
mpirun -np 4 build/babelsim-solve -case cases/cavity/benchmark-re1000-n128 -time throughput-np4
```

输出位于 `results/<time>/rank-*`。并行输出可独立后处理：

```bash
build/babelsim-post -case cases/cavity/benchmark-re1000-n128 -time throughput-np4 -format vtk tecplot
```

本机实测数据、热点和强缩放解释见 `docs/reports/re1000-n128-mpi-performance.md`。
