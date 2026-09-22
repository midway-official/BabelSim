# BabelSim

BabelSim 是一个面向 CFD、PDE 与多物理场计算的紧凑 C++17 有限体积框架。它以
TaihoCFD 的已验证数值算法为参考，但不依赖 TaihoCFD 的输入格式、运行时数据结构
或求解器实现。

框架的核心边界是：

```text
Mesh      空间在哪里，以及单元、面、顶点如何连接
Field     空间中存放什么数据
geometry  将网格几何量以普通 Field 暴露给 Physics
Operator  数据之间进行什么数学运算
Method    该运算采用什么离散方式
Equation  表达要求解的数学方程
equ       向绑定未知场的 Equation 立即加入离散项；显式组装后求解
math      描述并计算已有场上的数学量
Case      提供命名场、物性、配置与显式结果写出
time      显式时间推进和历史场；求解器使用普通 C++ 循环
monitor   通用观测值报告、输出进程和打印周期；不判断收敛
RunTime   内部时间推进、运行域与计算后端生命周期；普通 Solver 不构造它
FVM       数值前端：表达式解释、时间历史、离散方程与数值工作区
Backend   计算后端：整场同步、全局归约、稀疏装配和线性求解
Physics   如何组合方程与算子解决具体物理问题
```

`Equation` 的数学表达与内部离散 LDU 系统明确分离。Heat 与 SIMPLE 的 Physics 源码不接触
MPI、halo、CSR/LDU、Eigen 或 Field 底层存储；FVM 只经粗粒度 ComputeBackend 接口调用
这些能力，默认 Eigen/MPI 后端负责具体实现。

Solver Programming Model 正式分为两种组织方式：Heat、Diffusion、Poisson 和标量输运
采用 **Equation-driven**，核心源码就是一个或少量 PDE；SIMPLE 和耦合算法采用
**Algorithm-driven**，用多个 Equation 与 Correction 直接表达算法流程。两者共享同一套
Field、`equ/math`、离散、线性代数和 MPI Runtime，不建立两套 Framework。
稳态与瞬态 SIMPLE 各自在独立 main.cpp 展开所有算法步骤和 Rhie–Chow，不共享 SIMPLE 实现。

当前实现使用显式 face-based 的三维非结构 polyhedral 网格；每个 face 有一个 owner 和至多一个
neighbour，face 顶点数与 cell 面数都可变。薄域问题仍是三维层，不会维护独立的二维算子或二维求解器。
网格在构建时预计算体积、逆体积、中心、面积向量、单位法向、正交系数、非正交修正向量、偏斜量和
插值权重，以少量内存换取迭代热点中的计算速度。旧 hex 算例通过一次性 v2→v3 工具迁移，运行时
reader 只接受显式 face 的 v3 格式。

已实现：

- 通用三维 polyhedral cell/face/vertex 拓扑、边界 patch、可变长度 CSR 连接和只读 range；
- 连续存储的 scalar/vector/tensor Field 与通用边界条件；
- Gradient、Interpolation、Flux、Divergence、Convection、Diffusion、Laplacian、
  TimeDerivative 等有限体积算子；
- 三维非正交/偏斜修正：Least-Squares、修正 Green--Gauss、修正面插值、非正交
  扩散与压力法向梯度，以及面通量和中心对流的一致重构；
- LDU 方程、稀疏装配、串行与分布式 CG/BiCGSTAB，以及只作为预条件器的 AMG；
- MPI Krylov 使用按 global cell ID 的稀疏 halo matvec 和融合归约；MPI AMG 具有跨 rank 的图聚合粗网格；
- 框架级 MPI：局部 owned/ghost cell、halo exchange、分布式 matvec 与全局归约；
- 分布式网格读取：rank 0 解析原生 `.mesh`，并由并行层按单元邻接图构造每个 rank
  的 owned+ghost 局部 Mesh；
- 不可压 SIMPLE（稳态 `simple`、瞬态 `transientSimple`）与瞬态 PISO（`piso`）；
  动量插值（Rhie–Chow）与压力修正作为各自求解器的私有数值步骤，不下沉到通用层；
- RANS 湍流：Wilcox 1988 k-ω、标准 k-ε、Spalart–Allmaras，经 `rans::load` /
  `effectiveViscosity` / `deviatoricStressRemainder` / `solveTransport` 与动量方程耦合；
- 过程式 `equ/math` 编程模型：方程只有显式组装一层（`equ::createEquation` + 逐项装配），
  离散项写入绑定未知量的方程、场运算立即求值，`equ::solve` 只求解已组装方程，没有延迟表达式；
- `heat` 常物性瞬态热传导入口；通用方程 API 同时支持常数或 Field 系数；
- `transport` 瞬态对流-扩散 Solver，使用标量 Field、`equ::ddt/div/laplacian` 与边界；
- 对流支持一阶迎风、梯度重构的二阶 `linearUpwind` 和中心格式，扩散保持中心型有限体积离散；
- 原生 case/mesh/field 文件、通用并行结果写出与独立 VTK/Tecplot 后处理。

完整的 DSL、运行时 API、全部算子、配置键与文件格式见
[DSL 与运行时用户手册](docs/dsl-runtime-manual.md)；[架构与维护边界](docs/architecture.md)
说明公开接口的分层约束：Equation 只描述绑定未知量和离散项，内部离散存储与 Physics 分开，
SIMPLE 状态只归算法，线性控制只归运行配置。`make test-architecture` 自动检查项目头依赖和
分层约束，文档总入口见 [docs/](docs/README.md)。

## 构建与运行

依赖：C++17 编译器、Eigen 3、MPI-3 实现和 GNU Make。默认配置面向 GCC 工具链：
`mpic++` 调用 GCC，`gcc-ar` 归档 LTO 对象，正式程序与显式构建的测试使用相同 ABI。

默认 `make` 只构建 `build/libbabelsim.a`、`build/babelsim-solve` 和 `build/babelsim-post`，
不编译或运行任何测试。测试与验证必须另行显式调用 `make test*` / `make validate*` 的具体目标。

默认优化为 `-O3 -march=native -mtune=native -flto=auto -ffat-lto-objects
-ffast-math -fno-finite-math-only -ffp-contract=off -DNDEBUG`：启用本机 CPU 优化、
跨文件优化、浮点重结合与倒数优化，关闭融合乘加和调试断言。仍使用 double，
保留 NaN/Inf 检查及显式参数/收敛校验；不保证严格 IEEE 运算顺序或逐位一致。
快速数学还可能改变极小数、舍入和溢出行为，因此既有数值验证结论不能直接替代本配置的验证。
这些选项的含义参见 [GCC 优化选项](https://gcc.gnu.org/onlinedocs/gcc-11.4.0/gcc/Optimize-Options.html)。

`-march=native` 产物只适合本机及支持相同指令集的节点；异构集群应改用共同 CPU 基线并重编译。
fat LTO 让静态库保留普通机器码，外部程序可用 `-fno-lto` 禁用链接时优化。
若需要不放宽浮点规则的构建，可执行 `make clean`，然后使用
`make -j4 OPTFLAGS='-O3 -march=native -mtune=native -flto=auto -ffat-lto-objects'`。
修改命令行编译选项后应先清理；`make clean` 仅删除选定的 BUILD 目录，不删除 Case 结果。

计算后端采用构建期替换，避免运行时注册和热循环虚分派。框架维护者可令
`COMPUTE_BACKEND_SOURCES='src/backend/other.cpp ...'`；该源文件组实现内部
`makeComputeBackend()` 工厂及所需代数能力即可，并会整体排除默认 Eigen 装配/求解源码。
普通 Solver 作者不需要看到或选择这个接口，且替换后端不应修改 Physics、`equ/math` 或
FVM 离散源码。当前没有承诺动态插件或稳定后端 ABI。

```bash
make -j4  # 只构建，不测试；按可用内存调整并行编译数

# 串行腔体
build/babelsim-solve -case cases/cavity

# MPI 通道流
mpirun -np 4 build/babelsim-solve -case cases/poiseuille

# 瞬态热传导：串行与 MPI 使用同一 Solver 源码
mpirun -np 2 build/babelsim-solve -case cases/heat

# 独立读取网格和并行结果，输出 ParaView/Tecplot 文件
build/babelsim-post -case cases/poiseuille -format vtk tecplot

# 自动扫描 results/<time>，生成每个时刻的 VTK 和 ParaView 时间序列 post/series.pvd
build/babelsim-post -case cases/heat -time all -format vtk
```

`-time <名称>` 为独立运行命名：最终状态保留在该目录，瞬态序列位于其 `<物理时间>/` 子目录。例如：

```bash
build/babelsim-solve -case cases/poiseuille -time serial
mpirun -np 2 build/babelsim-solve -case cases/poiseuille -time mpi2
python3 tools/compare_parallel_results.py \
  cases/poiseuille/results/serial cases/poiseuille/results/mpi2 \
  --atol 5e-6 --rtol 5e-6
```

## 案例目录

每个案例自包含，不照搬 OpenFOAM 的复杂层级：

```text
cases/poiseuille/
├── case.bs                    # 选择求解器与各文件的相对路径
├── mesh/poiseuille.mesh       # 显式 Hex 几何、拓扑与 patch 名称和角色
├── fields/initial/U.field     # 初值与 U 的边界条件
├── fields/initial/p.field     # 初值与 p 的边界条件
├── physics/simple.bs          # 密度、黏度等物性
├── numerics/methods.bs        # 时间、方程/项/算子的具名离散格式
├── numerics/solution.bs       # SIMPLE/线性求解控制
├── control.bs                 # 时间区间与步长
├── output.bs                  # 结果目录、时刻名与可选字段筛选
└── results/<time>/rank-0000/  # 运行生成，不纳入 Git
```

每个 MPI rank 仅写出 owned cell 的 `U.csv`、`p.csv` 与 `metadata.bs`；ghost cell
不会写出。`babelsim-post` 按 global ID 检查完整性并合并为真实 polyhedron 的 VTK
XML `.vtu`；`-time all -format vtk` 还会产生 ParaView 可直接打开的
`post/series.pvd`。Tecplot `FEBRICK` 输出保留给八顶点六面体兼容结果；任意面数的网格使用 VTK
输出查看。

库代码需要从已存在的全局网格分区时仍可使用 `decompose()`；启动器和文件型并行程序
应使用 `readDistributedMesh(path, parallel)`。该接口在 rank 0 读取网格，Parallel 层
根据单元邻接关系建立 owned+ghost 几何、processor patch 和 global ID 映射；Physics 与
方程层不接触该分区细节。

## 测试与验证

```bash
make test                 # 几何、算子、case/field IO、Heat/标量输运/SIMPLE、通用输出
make test-architecture    # 头依赖与分层门禁（Physics 不得越界依赖）
make test-workflow        # 新 Solver 单函数开发、双场耦合、时间序列、真实 ParaView 读取
make test-external        # 仓库外 Solver 构建、1/2/4 进程、负向 API 编译、无 MPI 结果读取器
make test-mpi             # MPI 网格、halo、算子、线性求解、SIMPLE 与标量输运
make test-mpi-heat        # 1/2 rank 热传导场比较
make test-mpi-poiseuille  # 1/2/4 rank 案例启动器、结果比较、后处理
make test-rans            # RANS 方程与模型常数
make validate-cavity      # Re=100、二阶迎风的 Ghia 腔体快速回归
make validate-poiseuille  # 收敛的 Poiseuille 解析解比较
```

非正交离散的标准入口是 `interpolation corrected`、`gradient leastSquares`、
`diffusion corrected`；网格角度较大时可改用 `diffusion limitedCorrected`，并通过
`nonOrthogonalCorrections` 设置压力修正右端项的显式迭代次数。串行测试还包含
三维扭曲腔体，MPI 测试包含分区界面上的修正通用算子和 SIMPLE 私有耦合路径。

新增独立 Solver：自己的一个 C++ 源文件，用一行 `SolverRegistration` 注册名称/函数，
通用 main 调用 `runApplication(argc, argv)`，再准备 Case。
只链接公开头和预编译库，不修改 BabelSim 核心或内置启动器。若希望加入内置命令，
则新增 `src/physics/<name>/main.cpp`，在该文件注册自己，Makefile 自动收集，不再修改启动器名单。
不需要专用 Case reader、RunTime、并行输出代码、注册宏或 Solver 基类。
Heat、transport 的完整入口各自是一个短函数；SIMPLE 主循环明确列出五个算法步骤。
输入场、命名中间场及其生命周期由 Case 管理，`solver.h` 不包含 Runtime/MPI/代数实现头。
矢量场源、方程欠松弛、标量参考规范和动量对角响应均有公开数学入口。
Field 原始指针/索引及 Mesh 缓存/分区修改已限制到内部维护接口；按位置定义场可用 evaluate。
Case 的 validate 只校验，`time::start` / `start` / `setTime` 才关闭声明阶段；文件加载场
默认输出，程序创建的场使用 `create*Field`，已创建场使用 `existing*Field`，派生 cell 场
可通过 output(field) 选择输出。
公开 math 统一为整场同步契约，结果读取头和实现均不再要求 MPI。

默认瞬态结果按 `output.bs` 中 `writeInterval` 保存（省略时每步写出），最终时刻总会保存。
`-time mpi4/all` 后处理命名运行的完整序列；ParaView 打开对应的 `post/mpi4/series.pvd`。
内置求解器注册名为 `heat / transport / simple / transientSimple / piso`（RANS 由动量方程
求解器按 `physics` 字典的 `turbulenceModel` 启用）。内置方程通过稳定名称绑定
`equation.<name>.*` 离散与线性配置；每个实际使用的方程都必须写出完整线性配置。
详见 [按方程与算子配置](docs/numerical-configuration.md)。

文档入口见 [docs/README.md](docs/README.md)：

- [DSL 与运行时用户手册](docs/dsl-runtime-manual.md)：写求解器的唯一手册。最小可运行示例、
  Case 与全部配置键、场/网格文件格式、Field/geometry/math/equ 全部算子、时间与历史、
  线性求解契约、诊断与监视、结果与后处理、并行边界、内置求解器与 RANS、开发检查清单、
  当前格式迁移规则。
- [内置求解器手册](docs/solvers.md)：每个内置求解器与 RANS 模块的方程、配置键、运行方式、
  验证证据与验证边界。
- [架构与维护边界](docs/architecture.md)：分层、依赖禁令、所有权、维护流程与验收命令。
- [验证与维护检查](docs/validation.md)：串行/并行验证入口与新 Solver 的最低验收线。
- [性能工具说明](docs/performance/README.md)：`-performance` JSON、构建开关与 benchmark 驱动。

历史证据归档在 [docs/reports/](docs/reports/)：多 Reynolds 数、网格无关性、格式与 MPI 对照见
[Ghia 方腔验证报告](docs/reports/cavity-ghia-validation.md)及其
[PDF 版本](docs/reports/cavity-ghia-validation.pdf)；线性后端性能对照见
[分布式后端优化报告](docs/reports/backend-performance-optimization.md)。归档报告中的接口与
性能数字反映当时提交状态，当前后端只提供 CG/BiCGSTAB 与 IC/ILUT/AMG 预条件组合。
