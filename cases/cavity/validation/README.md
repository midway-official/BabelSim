# 方腔专项验证工具

本目录只保存二维顶盖驱动方腔的案例生成、Ghia 定量校验和绘图程序。通用 Solver、网格、
并行和后处理实现不依赖本目录。

## 文件职责

- `generate_cavity_case.py`：从 `cases/cavity` 生成指定 Re、网格数、壁面加密强度和迎风格式的案例；
- `validate_cavity.py`：读取串行 `U.csv` 或自动合并 `rank-*/U.csv`，计算 Ghia 中心线误差；
- `plot_cavity_validation.py`：绘制多个 Re 的 BabelSim 曲线与 Ghia 散点；
- `plot_cavity_convergence.py`：绘制 Ghia 误差和相邻网格整线差异；
- `plot_cavity_fields.py`：绘制速度、流线和涡量；
- `plot_cavity_history.py`：从 SIMPLE 日志绘制收敛历史；
- `plot_cavity_schemes.py`：比较同一 Re 下不同格式、网格或参数的中心线。

完整数据、误差定义、运行命令、图表和结论见
[`docs/reports/cavity-ghia-validation.md`](../../../docs/reports/cavity-ghia-validation.md)及
[`PDF 报告`](../../../docs/reports/cavity-ghia-validation.pdf)。

## 最小复现

内置 64² Re=100 案例使用 4 个 MPI rank 运行并自动检查 Ghia RMS：

```bash
make validate-cavity
```

生成新的 128²、Re=1000 壁面加密案例并使用 8 个 rank：

```bash
python3 cases/cavity/validation/generate_cavity_case.py \
    /tmp/cavity-Re1000-N128 --cells 128 --re 1000 --cluster 1.5

TMPDIR=/tmp mpirun -np 8 build/babelsim-solve \
    -case /tmp/cavity-Re1000-N128 -time validation

python3 cases/cavity/validation/validate_cavity.py \
    /tmp/cavity-Re1000-N128/results/validation --re 1000
```

脚本要求 Python 3 和 NumPy；绘图脚本另需 Matplotlib。生成器拒绝覆盖已有输出目录，避免误删结果。
