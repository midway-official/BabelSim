# BabelSim 数值验证与回归测试

本文只记录已收敛或明确标注为回归的结果。二维 Ghia 腔体和解析 Poiseuille 是
定量精度验证；小尺度三维和非正交腔体用于验证几何、守恒、对称性与重构流程，
尚不宣称与外部高分辨率基准完全一致。

## 1. 自动测试

### 具体 Solver 接口私有化与瞬态 SIMPLE 回归（2026-09-05）

本轮删除 `include/babelsim/simple.h` 和 `include/babelsim/simple_control.h`，并确认不存在
`transient_simple.h`。Heat、Transport、稳态 SIMPLE 与瞬态 SIMPLE 均不再发布具体 Solver
接口；前两者保持单文件方程入口，后两者在各自 Physics 目录中保留私有算法文件。
公共 `Case/Field/eqn/math/solve`、FVM、Runtime、MPI、代数和存储实现均未改动。

执行：

```bash
make -j4 test
make test-workflow
make test-external
make test-mpi
make test-mpi-poiseuille
```

上述目标全部通过。架构检查覆盖 69 个生产源文件/头文件：Physics 对 `internal/`、Runtime、
MPI、矩阵、原始 Field 存储的直接依赖为零，并拒绝重新加入具体 Solver 公共头。仓库外测试
区分两类证据：普通新 Solver 只复制公开 SDK；完整稳态 SIMPLE 则作为私有模块维护对象复制，
不把它伪装成公共二次开发接口。

另用临时 16×16 腔体 Case 验证 `transientSimple`，不向仓库增加专用 Case：Euler 三个物理
时间步在 1/2/4 rank 均收敛且各 rank 校正次数一致；以 1 rank 为参考，最终场最大绝对差为：

| 对照 | U | p |
| --- | ---: | ---: |
| 1 对 2 rank | 2.5650228e-10 | 6.6396855e-10 |
| 1 对 4 rank | 4.1818399e-10 | 2.8137813e-9 |

同一临时 Case 的 BDF2 三个物理时间步也正常收敛。该测试验证当前小网格的时间循环、历史、
全局停止行为与并行一致性，不替代高雷诺数瞬态基准、时间步无关性或物理精度验证。

### eqn/math 命名迁移回归（2026-08-31）

以 `dfc73e5` 为改名前基线，统一公开命名空间、头文件、实现文件、测试及文档。
本次仅调整名称和注释，不改变数值公式、离散格式、容差、内存布局或通信顺序。

```bash
make -j2 all test test-mpi test-workflow test-external
```

上述目标全部通过。新增独立包含 `math.h`、`eqn.h` 的普通 g++ 编译检查，
以及旧命名空间不可用的两项负向检查；原有十项维护接口隔离检查继续通过。

| 检查 | 本次结果 |
| --- | --- |
| 源码边界 | 65 个生产源文件/头文件的包含图无环；公开旧名称无残留 |
| 显式同步 | 22 项 math 运算 × 3 种扩散方法，污染 halo 后与串行对照；1/2/4 rank 最大差异为 0 / 7.42e-14 / 6.94e-14 |
| SIMPLE 回归 | 1/2/4 rank 腔体均 137 次；2 rank Poiseuille 为 865 次；三维及非正交测试通过 |
| 改名前后场对照 | 同进程数、同配置的腔体（1/2/4 rank）、三维腔体（2 rank）、Poiseuille（2 rank）的 U/p，以及 Heat（2 rank）的 T，与预先保存的结果最大绝对差均为 0 |
| 外部 Solver | 方程驱动、耦合和矢量示例只使用公开头文件；完整 SIMPLE 作为私有模块维护回归，均通过 1/2/4 rank 测试 |
| 完整工作流 | Heat/transport/耦合时间序列、失败路径、ParaView PVD/VTU 实际读取通过 |

场对照按 global ID 对齐，逐分量计算最大绝对差，以 `atol=rtol=1e-12` 验收；
这里的零差异只指本次相同进程数的改名前后结果，不表示不同进程数必须逐 bit 一致。
本次没有重做性能或 sanitizer 测试。以下各节为先前阶段的历史证据，不替代本次回归。

### 1.0 公开契约与仓库外二次开发验收（2026-08-31）

本轮基线为 `be18f32`，以下是本轮新增证据；后续小节保留前几轮的历史记录，
其中旧文件数、旧维护 API 状态不代表当前版本。

执行了强制重新编译，避免复用旧对象文件：

```bash
make -B -j2 all test test-workflow test-external test-mpi test-mpi-poiseuille
```

随后加入点值 Field 映射、输入/响应别名检查与应用失败返回码测试，再运行完整普通回归、
工作流、外部接入和 MPI 回归。最后应用返回码收紧后重跑 all/test/test-workflow/test-external，
均通过；没有放宽已有数值容差或增加 SIMPLE 外迭代。

| 验收 | 本轮结果 |
| --- | --- |
| 真正仓库外构建 | 临时目录仅复制 include/、预编译库和外部 Solver；没有 -Isrc 或框架源文件 |
| Solver 编译隔离 | 普通 g++、C++17、-Wall -Wextra -Werror 编译；不需要 MPI/Eigen 头 |
| 新 Solver | 方程驱动输运、算法驱动双标量耦合、矢量源/响应示例均在 1/2/4 进程检查解析值 |
| Field/Case | 非均匀位置源、U→动能点值映射、张量派生输出、面矢量/张量创建、声明阶段与稳定引用 |
| 方程公开能力 | 矢量 Field 体源、标量欠松弛、Poisson 参考规范、V/aP 响应及非法别名拒绝 |
| 负向编译 | 原始 Field 指针/索引、Mesh 数组/分区/替换、离散矩阵入口、已移除的单面局部核、内部头等十项被拒绝 |
| 启动失败路径 | 负返回码与其他进程的 0 不会合并为成功；空表和重复名称拒绝；失败不写成功结果 |
| 结果读取隔离 | 普通 g++ 链接并实际读取并行结果，未链接 MPI |
| 显式同步契约 | 故意污染 ghost/重复面值，17 项 math 运算对照串行；1/2/4 进程最大差异为 0 / 7.42e-14 / 6.94e-14 |
| 源码边界 | 67 个 src/include 文件的包含图无环；FVM 执行实现不包含 RunTime；Application 无 MPI/格式实现 |
| 数值回归 | 通道/2D腔体/3D腔体/2D非正交/3D非正交迭代仍为 49/137/57/164/79 |
| MPI SIMPLE | 腔体 1/2/4 rank 均 137 次；Poiseuille 均 865 次 |
| MPI 场差异 | Poiseuille 1 对 4 rank：U 6.1705172e-7、p 1.5343525e-6；Heat 1 对 2 rank：T 3.7990135e-11 |
| 后处理 | 原时间序列/VTU/Tecplot 工作流及真实 ParaView PVDReader 通过 |
| 定向 ASan/UBSan | 新构建的公开方程、Case 生命周期及 2 rank 非正交 math 测试通过；后者差异 7.05e-14 |

Sanitizer 使用 `-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer`，
运行设置 `ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1 TMPDIR=/tmp`。
关闭泄漏检查，不能据此宣称全程序或 MPI 库无泄漏；没有覆盖所有算法和进程规模。

SIMPLE 过程中允许某一轮线性容差未满足而继续下一外迭代；现在公开 solve 会记录这类诊断。
最终 converged 必须同时满足数值健康、该轮全部线性求解及外迭代判据。
不能把日志中的中途失败解释成“所有线性调用都通过”，也不能把最终符合容差的结果误作失败。
场差异按 global ID 对齐后取 cell/分量的最大绝对差，不是残差。

本轮没有做大规模 HPC 扩展性或性能加速比测试。原数组布局、几何缓存、面遍历及线性后端
保持不变；新维护访问函数为内联访问，数学描述不复制大型矩阵。仍存在离散数组分配，
不宣称零分配、任意多区域、通用约束系统或张量 PDE 求解已经实现。

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

- 当时 SIMPLE 公开头由 187 行降到 40 行；当前版本已进一步删除全部具体 Solver 公开头，
  SIMPLE 只保留 Physics 模块私有算法接口。Heat、SIMPLE、transport 主入口分别为
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
- Green-Gauss/Least-Squares 梯度、插值、重构、通量、散度、中心/一阶/二阶迎风对流、
  正交/非正交扩散、Laplacian、Euler/BDF2；
- Runtime 隐藏的离散 LDU 装配、Field 系数时间/扩散项、稀疏结构复用与串行线性求解；
- `case.bs`、物性/数值字典、花括号 `.field`、通用 scalar/vector/tensor 结果写出与
  result reader；
- Heat、标量输运和 SIMPLE 对同一 `eqn/math/solve` 路径的复用检查；
- Case 场引用稳定性、确定析构顺序、时间历史按物理步推进，以及矢量方程全分量收敛；
- 新双场耦合 Solver 的单文件开发路径、1/2/4 rank 一致性、多时间步写出、
  失败路径、XML VTU/PVD 和真实 ParaView 读取；
- `math::subtract` 的 cell 梯度修正和 face 扩散通量修正，确保 SIMPLE 不以手写 Field 索引
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

## 2. 顶盖驱动腔体

运行：

```bash
make validate-cavity
```

内置案例为稳态不可压 SIMPLE，网格 `64 × 64 × 1`，二阶 `linearUpwind` 对流、
中心型正交扩散、速度收敛阈值 `1e-6`。通过
`mpirun -np 4 babelsim-solve -case cases/cavity -time validation-mpi4` 启动，在
2522 次迭代后收敛：

```text
相对质量不平衡 = 1.2725e-14
相对速度变化   = 9.9929e-07
```

与 Ghia 中心线样本比较：

| 指标 | 数值 |
|---|---:|
| 水平速度中心插值 | `-0.20818258` |
| 垂直速度中心插值 | `0.05740706` |
| 水平速度 L∞ 误差 | `0.004068828` |
| 水平速度 L2 误差 | `0.001775332` |
| 垂直速度 L∞ 误差 | `0.008431023` |
| 垂直速度 L2 误差 | `0.004352236` |

该目标是快速回归，不是网格无关性证明。Re=100/400/1000、32²/64²/128²、壁面加密、
格式对照、MPI 一致性和可视化结果见
[二维顶盖驱动方腔 Ghia 验证报告](reports/cavity-ghia-validation.md)及其
[PDF 版本](reports/cavity-ghia-validation.pdf)。

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
- 2 rank 对流--扩散标量输运；该 Solver 与串行版本使用完全相同的 `eqn::ddt +
  eqn::div == eqn::laplacian + source` 源码。

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
