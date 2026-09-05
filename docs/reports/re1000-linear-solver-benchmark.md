# Re=1000 顶盖驱动流：线性后端性能对照

## 实验条件

实验使用 [benchmark-re1000](../../cases/cavity/benchmark-re1000) 案例：`64×64×1`、
双曲正切壁面加密 `cluster=1.5`、`Re=1000`、`linearUpwind` 对流、最小二乘梯度和中心
型扩散。两套方案的线性容差、SIMPLE 松弛和外迭代终止条件相同。测试主机为 AMD RYZEN
AI MAX+ 392，12 个物理核心/24 个逻辑 CPU；实验只使用 1、2、4 个 MPI rank，没有把
超线程当作物理并行度。

## 方案

| 名称 | 速度方程 | 压力方程 |
| --- | --- | --- |
| 基线 | `BiCGSTAB + ILUT` | `CG + IncompleteCholesky` |
| AMG | `GMRES + AMG` | `GMRES + AMG` |

AMG 方案使用 `gmresRestart=30`、`amgMaxLevels=12`、`amgCoarseSize=48`、
`amgSmoothingSteps=2`。每个 run 都返回成功状态，SIMPLE 外迭代均为 3651 次；尾部日志
中的 `dU≈9.996e-7`、`converged=true` 与案例的 `velocityTolerance=1e-6` 一致。

## 实测结果

墙钟时间由 `/usr/bin/time` 测量，包含启动、网格读取、矩阵准备、线性求解和结果写出。
`peak_kib` 是该进程的峰值常驻内存；MPI 结果列出单个 rank 的峰值。

| 方案 | MPI rank | wall (s) | peak (KiB) | 状态 | 外迭代 |
| --- | ---: | ---: | ---: | --- | ---: |
| 基线 | 1 | 68.51 | 28704 | 收敛 | 3651 |
| 基线 | 2 | 40.73 | 25536 | 收敛 | 3651 |
| 基线 | 4 | 24.78 | 25152 | 收敛 | 3651 |
| GMRES+AMG | 1 | 94.35 | 39300 | 收敛 | 3651 |
| GMRES+AMG | 2 | 114.83 | 28792 | 收敛 | 3651 |
| GMRES+AMG | 4 | 217.20 | 25344 | 收敛 | 3651 |

基线从 1 到 4 rank 的实测加速约为 `2.76×`。在本小规模、强耦合 SIMPLE 案例上，当前
轻量聚合 AMG 的绝对时间没有优于 ILUT：AMG 每个外迭代都要为变化的动量/压力系数更新
Galerkin 粗层并进行局部 V-cycle，4 rank 时局部问题较小，层级准备和全局同步开销占比
反而上升。因此不能把“实现了 AMG”误写成“该案例已经获得 AMG 加速”。

这组结果仍验证了三个重要性质：

- 新后端没有改变 Physics、FVM 或 `eqn/math/solve` API；
- 1/2/4 rank 使用相同的全局收敛判据并成功完成，未出现额外外迭代；
- GMRES 的基向量、Hessenberg/Givens 工作区和 AMG 层级在后端准备阶段分配/复用，MPI
  matvec 与全局归约仍保持在计算后端。

## 结论与后续优化

当前实现适合作为可替换、可验证的 AMG/GMRES 后端基线；对大规模三维扩展，应该在更大
局部问题上重新调节粗层阈值、平滑次数和重启维数，并进一步实现全局粗网格/并行聚合，
再与 ILUT 做强扩展和弱扩展对照。小网格结果只说明正确性和当前开销组成，不代表 AMG
在所有 CFD 问题上都更快。
