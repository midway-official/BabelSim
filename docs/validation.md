# BabelSim 数值验证与回归测试

本文只记录已收敛或明确标注为回归的结果。二维 Ghia 腔体和解析 Poiseuille 是
定量精度验证；小尺度三维和非正交腔体用于验证几何、守恒、对称性与重构流程，
尚不宣称与外部高分辨率基准完全一致。

## 1. 自动测试

```bash
make test
make test-mpi
make test-mpi-poiseuille
make validate-cavity
make validate-poiseuille
```

`make test` 覆盖：

- `nz=1` 二维退化、三维六面体几何、非正交与偏斜量；
- scalar/vector/tensor Field、FixedValue、FixedGradient、ZeroGradient、
  InletOutlet、Symmetry/Mirror；
- Green-Gauss/Least-Squares 梯度、插值、重构、通量、散度、中心/迎风对流、
  正交/非正交扩散、Laplacian、Euler/BDF2；
- LDU 装配、可变扩散系数、稀疏结构复用与串行线性求解；
- `case.bs`、物性/数值字典、花括号 `.field`、通用 scalar/vector/tensor 结果写出与
  result reader；
- 独立的 `MomentumInterpolation` 和 `PressureCorrection` 基本代数检查；
- 三维仿射非正交制造解上的两个专用算子精确性检查；
- 小通道、二维腔体、三维腔体、二维扭曲和三维扭曲非正交腔体的 SIMPLE 回归。

代表性回归结果：

| 案例 | 网格 | SIMPLE 迭代 | 相对质量不平衡 | 其他门限 |
|---|---:|---:|---:|---|
| 小通道 | 9 | 49 | `1.59e-16` | `max(abs(w)) < 1e-13` |
| 二维腔体 | 12 × 12 × 1 | 137 | `1.61e-16` | 中心 `u=-0.149236` |
| 三维腔体 | 6 × 6 × 6 | 57 | `1.59e-15` | 中心 `u=-0.114134`，`max(abs(w))=0.0192033` |
| 非正交腔体 | 10 × 10 × 1 | 164 | `9.77e-9` | `max(abs(k_nonorth))/abs(Sf)=0.412398` |
| 三维非正交腔体 | 5 × 5 × 5 | 79 | `1.86e-8` | `centreU=-0.0940683`，`max(abs(w))=0.0108384` |

## 2. Re=100 顶盖驱动腔体

运行：

```bash
make validate-cavity
```

案例为稳态不可压 SIMPLE，网格 `64 × 64 × 1`，一阶迎风对流、正交扩散、速度收敛
阈值 `1e-6`。通过 `babelsim-solve -case cases/cavity -time validation` 启动，在
2355 次迭代后收敛：

```text
相对质量不平衡 = 1.0811e-14
相对速度变化   = 9.9889e-07
```

与 Ghia 中心线样本比较：

| 指标 | 数值 |
|---|---:|
| 水平速度中心插值 | `-0.19697720` |
| 垂直速度中心插值 | `0.05108802` |
| 水平速度 L∞ 误差 | `0.011111441` |
| 水平速度 L2 误差 | `0.005884307` |
| 垂直速度 L∞ 误差 | `0.007404026` |
| 垂直速度 L2 误差 | `0.003986802` |

该结果对应当前 64²、一阶迎风、指定松弛与线性容差；不应将它表述为高阶离散或
网格无关解。

## 3. Poiseuille 通道流

运行：

```bash
make validate-poiseuille
```

案例网格为 `62 × 70 × 1`，共 4340 个 cell。通过 case 驱动启动器在 865 次迭代后
收敛，出口速度与单位平均速度的抛物线比较为：

| 指标 | 数值 |
|---|---:|
| 相对质量不平衡 | `2.6603e-15` |
| 最大速度 | `1.49850478` |
| L∞ 误差 | `0.001189098` |
| L2 误差 | `0.000645703` |

结果来自 `cases/poiseuille/results/validation/rank-0000/U.csv`，验证脚本直接读取
通用 field 输出的 `value0` 分量；不再依赖旧的 U/p 合并示例格式。

## 4. MPI 框架级验证

运行：

```bash
make test-mpi
make test-mpi-poiseuille
```

`test-mpi` 验证的不是 SIMPLE 外层包装，而是通用并行基础：

- `Mesh` 的 owned/ghost/global-ID 映射和两层 halo；
- cell scalar/vector/tensor 与 interface face field 的 HaloExchange；
- 分区接口处的 Green-Gauss、Least-Squares、修正插值、修正通量、散度、中心对流和非正交扩散；
- 分区接口处的 `MomentumInterpolation` 与 `PressureCorrection`；
- 仅生成 owned 行的 SparseAssembly；
- halo matvec、全局点积/范数的分布式 Krylov 求解；
- 1/2/4 rank 二维腔体及 2 rank 三维腔体。

12 × 12 × 1 腔体结果：

| rank 数 | SIMPLE 迭代 | 相对质量不平衡 | 中心 `u` |
|---:|---:|---:|---:|
| 1 | 137 | `1.14e-14` | `-0.149236` |
| 2 | 137 | `1.49e-14` | `-0.149236` |
| 4 | 137 | `1.94e-14` | `-0.149236` |

2 rank 三维腔体在 57 次迭代后收敛，相对质量不平衡 `3.09e-15`，中心
`u=-0.114134`，`max(abs(w))=0.0192033`，并通过中心面对称性检查。

`test-mpi-poiseuille` 使用同一个 `cases/poiseuille` 分别运行 1、2、4 rank，按
global ID 比较通用 rank 结果：

| 比较 | 迭代 | 最大 `U` 差异 | 最大 `p` 差异 | 比较容差 |
|---|---:|---:|---:|---|
| 1 rank 与 2 rank | 865 / 865 | `6.15e-7` | `1.54e-6` | `atol=rtol=5e-6` |
| 1 rank 与 4 rank | 865 / 1131 | `4.45e-4` | `5.56e-5` | `atol=rtol=1e-3` |

4 rank 使用各进程局部预条件块，Krylov 与 SIMPLE 的迭代路径可与串行不同；因此
比较使用明确的绝对/相对容差，不要求逐 bit 一致。所有比较前会检查 global ID 无
重复、无遗漏。随后 `babelsim-post` 读取 rank 文件并生成原始六面体 VTK/Tecplot
文件，验证结果格式独立于求解器内存。

## 5. 内存安全与性能检查

以下命令用于串行 AddressSanitizer/UndefinedBehaviorSanitizer 回归：

```bash
make BUILD=build-sanitize \
  CXXFLAGS='-std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -Wall -Wextra -Wpedantic -Wshadow -DOMPI_SKIP_MPICXX=1 -DMPICH_SKIP_MPICXX=1' test
```

当前已执行 case IO、通用结果写出和算子测试的 ASan/UBSan 检查。性能设计重点是
预计算几何与稀疏模式、复用预条件器/工作向量、避免迭代中的大对象分配，以及只在
需要 ghost 数据的阶段通信。后续应在代表性三维网格上补充并行规模、内存带宽和
预条件器性能剖析，而不是把小网格回归时间误作 HPC 可扩展性结论。

## 6. 当前未完成项

- 更高分辨率三维腔体的外部基准比较与网格收敛研究；
- 非正交流动算例的独立精度基准；
- 瞬态全流程的时间阶验证；
- xyz/图分区、vertex halo、可扩展全局预条件器与 MPI-IO；
- 任意非结构化网格和 GPU 后端。

这些是扩展方向，不影响当前已验证的三维结构化有限体积与框架级 MPI 基础。
