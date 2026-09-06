# Re=1000、256×256 线性后端吞吐矩阵

这是固定问题规模的性能 Case：`256×256×1`、双曲正切壁面加密（`cluster=1.5`）、
`Re=1000`、二阶线性迎风对流、最小二乘梯度、线性面插值和正交扩散。

`numerics/solution.*.bs` 保存五种可复现的线性方案：

- `baseline`：速度 `BiCGSTAB+ILUT`，压力 `CG+IncompleteCholesky`；
- `bicgstab-ilut`：速度和压力均为 `BiCGSTAB+ILUT`；
- `gmres-ilut`：速度和压力均为 `GMRES(30)+ILUT`；
- `bicgstab-amg`：速度和压力均为 `BiCGSTAB+AMG`；
- `gmres-amg`：速度和压力均为 `GMRES(30)+AMG`。

CG+ILUT 不是 BabelSim 支持的组合，也不适合压力/动量方程的通用假设；基线中的 CG
只用于近似对称正定的压力校正方程。

`velocityTolerance=3e-3` 以及线性相对残差 `3e-3` 是固定工作量的吞吐判据，不是精度验证阈值。完整的 1/2/4-rank
运行时间、峰值内存和热点解释见 `docs/reports/re1000-n256-linear-solver-mpi.md`。
