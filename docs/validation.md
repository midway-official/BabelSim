# BabelSim 数值验证与回归测试

本文只记录已收敛或明确标注为回归的结果。二维 Ghia 腔体和解析 Poiseuille 是
定量精度验证；小尺度三维和非正交腔体用于验证几何、守恒、对称性与重构流程，
尚不宣称与外部高分辨率基准完全一致。

## 1. 自动测试

### 1.1 框架维护边界精简验收（2026-08-31）

基线为 `58442b5`。本轮删除重复物理入口、统一离散方程命名、清理 SIMPLE 反向依赖、
合并重复配置/解析以及分离应用分派后，执行：

```bash
make -j2 all test test-workflow test-mpi test-mpi-poiseuille
```

上述目标全部通过。最后一次头文件依赖收紧后再次运行相同命令，未修改数值容差、
压力/速度线性配置、非正交修正次数或并行停止规则。

| 检查 | 结果 |
| --- | --- |
| 生产 src/include 文件 | `64 → 58`，净减 6 个文件 |
| 生产代码行 | `8724 → 8424`，净减 300 行；排除空行和整行注释，不包含测试/文档 |
| 架构检查 | 全部 58 个源文件/头文件的项目包含图无环；无缺失头、底层反向包含算法状态或应用分派 |
| 公开头自包含 | 所有 include/babelsim 头只用 `-Iinclude`（另加 Eigen 路径）独立编译，严格警告加 `-Werror` 通过 |
| Case/方法选择 | 四种算子按 Field 覆盖值及默认值测试通过；重复配置在 2 rank 上明确失败 |
| 单节点 SIMPLE | 通道 49、二维腔体 137、三维腔体 57、二维/三维非正交腔体 164/79 次 |
| MPI SIMPLE | 1/2/4 rank 腔体均 137 次，Poiseuille 均 865 次；停止次数未改变 |
| MPI 场比较 | Poiseuille 1 对 4 rank：U `6.1705172e-7`、p `1.5343525e-6`；Heat 1 对 2 rank：T `3.7990135e-11` |
| 真实用户流程 | Heat/Transport/双标量耦合 1/2/4 rank、逐时间结果、失败路径、实际 ParaView PVD/VTU 读取通过 |
| 定向内存检查 | 新构建的 ASan/UBSan Case 生命周期、时间历史、SIMPLE 测试通过 |

场差异按相同 global ID 对齐后取所有 cell/分量的最大绝对差；不是迭代残差。
本轮生产文件没有块注释，代码行数由非空且非整行 `//` 的行统计。
没有增加新的框架类，也没有删除可复用工作数组以换取表面上的代码减少。

ASan/UBSan 使用 `-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer`，
运行时设 `ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1`。
关闭泄漏检测，因此只报告这些测试未发现地址/未定义行为错误，不声称已经证明无泄漏，
也不把单节点检查扩大为完整 MPI sanitizer 覆盖。时间历史测试中的一次
`vector equation failed` 是刻意设置的线性不收敛用例，测试同时确认不得误报成功。

静态检查和回归不能证明“所有潜在抽象泄漏不存在”；
[architecture.md](architecture.md) 逐项记录剩余维护接口及允许依赖，尤其是 Field 的
维护者存储 API、算法初始化/日志对 RunTime 的使用，以及 FVM 执行协调器的职责边界。
本轮没有重新测量 HPC 加速比，也没有把相同小网格结果当作大规模并行性能证明。

### 1.2 前一阶段 Solver 工作流验收（2026-08-31）

Solver 工作流阶段另外确认：

- SIMPLE 公开头由 187 行降到 40 行；Heat、SIMPLE、transport 主入口分别为
  20、18、21 行（含 include、空行和命名空间）。这只是可读性指标，不替代数值验收。
- 新双场耦合例子只新增一个普通函数；在 1/2/4 rank 上对每个保存时间层检查
  隐式 Euler 的解析递推结果，不要求新 reader、执行类或矩阵实现。
- 检查公共头文件闭包未引入 RunTime/MPI/Eigen；检查 Case scalar/vector/tensor 输入、
  稳定引用、中间场不输出、错误配置、失败步不保存及旧分区冲突。
- ParaView PVDReader 实际读取了 `[2, 4, 6, 8, 10]` 时间序列、32 个 cell 和温度值，
  不仅检查文件后缀或 XML 语法。
- 1/2/4 rank 小腔体均为 137 次，Poiseuille 均为 865 次；Poiseuille 1 对 4 rank
  的逐 global-ID 最大绝对差为 U 约 `6.17e-7`、p 约 `1.53e-6`。
  文件比较容差收紧为 `atol=rtol=5e-6`。Heat 1 对 2 rank 最大差约 `3.80e-11`。
- ASan/UBSan 下的时间历史、SIMPLE 单次迭代和 Case 生命周期测试通过；
  此项关闭了泄漏检测，不能据此宣称 MPI 库或整个程序已证明无泄漏。

测试还修正了原有 BDF2 首步缺少 Euler 启动、时间步内多次 solve 推进历史、以及分布式
Krylov 提前停止但真实残差略高于容差的边界问题。没有通过放宽线性容差或额外 SIMPLE
外迭代来让测试通过。

```bash
make test
make test-workflow
make test-mpi
make test-mpi-poiseuille
make validate-cavity
make validate-poiseuille
```

`make test` 与 `make test-workflow` 合起来覆盖：

- `nz=1` 二维退化、三维六面体几何、非正交与偏斜量；
- scalar/vector/tensor Field、FixedValue、FixedGradient、ZeroGradient、
  InletOutlet、Symmetry/Mirror，以及压力修正齐次边界映射；
- Green-Gauss/Least-Squares 梯度、插值、重构、通量、散度、中心/迎风对流、
  正交/非正交扩散、Laplacian、Euler/BDF2；
- Runtime 隐藏的离散 LDU 装配、Field 系数时间/扩散项、稀疏结构复用与串行线性求解；
- `case.bs`、物性/数值字典、花括号 `.field`、通用 scalar/vector/tensor 结果写出与
  result reader；
- Heat、标量输运和 SIMPLE 对同一 `fvm/fvc/solve` 路径的复用检查；
- Case 场引用稳定性、确定析构顺序、时间历史按物理步推进，以及矢量方程全分量收敛；
- 新双场耦合 Solver 的单文件开发路径、1/2/4 rank 一致性、多时间步写出、
  失败路径、XML VTU/PVD 和真实 ParaView 读取；
- `fvc::subtract` 的 cell 梯度修正和 face 扩散通量修正，确保 SIMPLE 不以手写 Field 索引
  绕过 Runtime 的非正交、halo 与连续存储路径；
- 小通道、二维腔体、三维腔体、二维扭曲和三维扭曲非正交腔体的 SIMPLE 回归。

代表性回归结果：

| 案例 | 网格 | SIMPLE 迭代 | 相对质量不平衡 | 其他门限 |
|---|---:|---:|---:|---|
| 小通道 | 9 | 49 | `1.77e-16` | `rho=2` 缩放对流；`max(abs(w)) < 1e-13` |
| 二维腔体 | 12 × 12 × 1 | 137 | `1.65e-16` | 中心 `u=-0.149236` |
| 三维腔体 | 6 × 6 × 6 | 57 | `1.59e-15` | 中心 `u=-0.114134`，`max(abs(w))=0.0192033` |
| 非正交腔体 | 10 × 10 × 1 | 164 | `9.74e-9` | `max(abs(k_nonorth))/abs(Sf)=0.412398` |
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
| 相对质量不平衡 | `2.6099e-16` |
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
- 分区接口处的修正插值、通量、散度、中心对流和非正交扩散；
- 仅生成 owned 行的 SparseAssembly；
- halo matvec、全局点积/范数的分布式 Krylov 求解；
- 1/2/4 rank 二维腔体、2 rank 原生网格通道流、2 rank 三维腔体，以及 1/2 rank 热传导。
- 2 rank 对流--扩散标量输运；该 Solver 与串行版本使用完全相同的 `fvm::ddt +
  fvm::div == fvm::laplacian + source` 源码。

文件型启动器还通过 `readDistributedMesh()` 验证根 rank 解析、尺寸/patch 广播、
局部顶点点对点传输和接收端重建；`parallel_channel_test` 使用该入口，不再先在
每个进程创建全局 Mesh。uniform `.field` 初值直接绑定各 rank 的局部 Field，因而
求解阶段不存在全局 Field 副本。

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
| 1 rank 与 4 rank | 865 / 865 | `6.17e-7` | `1.53e-6` | `atol=rtol=5e-6` |

Poiseuille 案例的压力方程使用 BiCGSTAB+ILUT，使分区后的内层线性方程也满足其全局残差
容差；因此 `SimpleIterationResult::converged` 只会在所有速度分量、所有压力非正交循环
均线性收敛，且连续性/速度外迭代判据同时满足时变为真。SIMPLE 外迭代的连续性、速度变化量
和健康状态由 `MPI_Allreduce` 的全局值共同决定，因此 1/2/4 rank 在该算例的停止点一致，
结果差异仅来自分区矩阵乘法、归约顺序和局部预条件器的浮点误差。所有比较前
会检查 global ID 无重复、无遗漏。随后 `babelsim-post` 读取 rank 文件并生成原始
六面体 VTK/Tecplot 文件，验证结果格式独立于求解器内存。

`test-mpi` 同时用 `cases/heat` 比较 1 与 2 rank 的瞬态温度场。五个 Euler 步后最大
`T` 差异为 `3.80e-11`（`atol=rtol=1e-9`），表明新的热 Solver 通过同一 Runtime 自动
使用 halo、分布式矩阵和全局线性收敛，而不是在 Solver 外额外包一层 MPI。

## 5. 内存安全与性能检查

以下命令用于串行 AddressSanitizer/UndefinedBehaviorSanitizer 回归：

```bash
make BUILD=build-sanitize \
  CXXFLAGS='-std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -Wall -Wextra -Wpedantic -Wshadow -DOMPI_SKIP_MPICXX=1 -DMPICH_SKIP_MPICXX=1' test
```

此前完整 `make test` 的 ASan/UBSan 记录覆盖 Field/算子、case IO、通用结果写出、Heat、标量输运、
SIMPLE、二维/三维与非正交回归；它是前一版本的验证记录，本轮重新执行的范围见 1.1 节。性能设计重点是
预计算几何与稀疏模式、复用预条件器/工作向量、避免迭代中的大对象分配，以及只在
需要 ghost 数据的阶段通信。后续应在代表性三维网格上补充并行规模、内存带宽和
预条件器性能剖析，而不是把小网格回归时间误作 HPC 可扩展性结论。

本次优化前后在同一开发机上对 `64 × 64 × 1` 腔体做了一次粗测（`-O3 -march=native`，
单 rank，包含完整 2355 次外迭代）：串行墙钟时间由 `22.34 s` 降至 `20.41 s`，约
`8.6%`。该数字只用于确认优化方向，受 CPU 频率和系统负载影响，不代表跨平台基准；
MPI 扩展性仍应使用更大三维网格和多次重复测量评估。

## 6. 当前未完成项

- 更高分辨率三维腔体的外部基准比较与网格收敛研究；
- 非正交流动算例的独立精度基准；
- 瞬态全流程的时间阶验证；
- xyz/图分区、vertex halo、可扩展全局预条件器与 MPI-IO；
- 任意非结构化网格和 GPU 后端。

这些是扩展方向，不影响当前已验证的三维结构化有限体积与框架级 MPI 基础。
