# BabelSim

BabelSim 是一个面向 CFD、PDE 与多物理场计算的紧凑 C++17 有限体积框架。它以
TaihoCFD 的已验证数值算法为参考，但不依赖 TaihoCFD 的输入格式、运行时数据结构
或求解器实现。

框架的核心边界是：

```text
Mesh      空间在哪里，以及单元、面、顶点如何连接
Field     空间中存放什么数据
Operator  数据之间进行什么数学运算
Method    该运算采用什么离散方式
Equation  表达要求解的数学方程
eqn       构造待求解方程的项，包括隐式项与已知源项
math      描述并计算已有场上的数学量
Case      提供命名场、物性、时间循环与自动结果序列
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
Field、`eqn/math`、离散、线性代数和 MPI Runtime，不建立两套 Framework。

当前实现使用统一的三维结构化六面体网格；二维问题是 `nz=1` 的退化三维网格，
不会维护独立的二维算子或二维求解器。网格在构建时预计算体积、逆体积、中心、面积向量、单位法向、
正交系数、非正交修正向量、偏斜量和插值权重，以少量内存换取迭代热点中的计算速度。

已实现：

- 三维结构化正交/非正交网格、边界 patch、cell/face/vertex 拓扑；
- 连续存储的 scalar/vector/tensor Field 与通用边界条件；
- Gradient、Interpolation、Flux、Divergence、Convection、Diffusion、Laplacian、
  TimeDerivative 等有限体积算子；
- 三维非正交/偏斜修正：Least-Squares、修正 Green--Gauss、修正面插值、非正交
  扩散与压力法向梯度，以及面通量和中心对流的一致重构；
- LDU 方程、稀疏装配、串行与分布式 Krylov 线性求解；
- 框架级 MPI：局部 owned/ghost cell、halo exchange、分布式 matvec 与全局归约；
- 分布式网格读取：rank 0 解析原生 `.mesh`，只发送各 rank 局部几何；不会在每个
  rank 重复保留全局 Mesh/Field；
- 不可压 SIMPLE；动量插值与压力修正作为其私有、具有独立数值语义的数值步骤；
- 轻量 `eqn/math` 编程模型：方程表达式只在 `solve()` 时离散，场运算按需执行，不复制大矩阵；
- `heat` 常物性瞬态热传导入口；通用方程 API 同时支持常数或 Field 系数；
- `transport` 瞬态对流-扩散 Solver，复用标量 Field、`eqn::ddt/div/laplacian` 与边界；
- 对流支持一阶迎风、梯度重构的二阶 `linearUpwind` 和中心格式，扩散保持中心型有限体积离散；
- 原生 case/mesh/field 文件、通用并行结果写出与独立 VTK/Tecplot 后处理。

维护者以 [架构与文件边界](docs/architecture.md) 为准：数学 EquationDefinition 与
DiscreteEquation 存储分开，SIMPLE 状态只归算法，线性控制只归运行配置。
Heat/Transport 不再维护重复的库式入口。`make test-architecture` 自动检查项目头依赖和分层约束。

## 构建与运行

依赖：C++17 编译器、Eigen 3、MPI-3 实现和 GNU Make。默认配置面向 GCC 工具链：
`mpic++` 调用 GCC，`gcc-ar` 归档 LTO 对象，正式程序与显式构建的测试使用相同 ABI。

默认 `make` 只构建 `build/libbabelsim.a`、`build/babelsim-solve` 和 `build/babelsim-post`，
不编译或运行任何测试。测试与验证必须另行显式调用 `make test*` / `make validate*` 的具体目标。

默认优化为 `-O3 -march=native -mtune=native -flto=auto -ffat-lto-objects
-ffast-math -fno-finite-math-only -ffp-contract=fast -DNDEBUG`：启用本机 CPU 优化、
跨文件优化、浮点重结合与倒数优化、融合乘加，并关闭调试断言。仍使用 double，
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
普通 Solver 作者不需要看到或选择这个接口，且替换后端不应修改 Physics、`eqn/math` 或
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
├── mesh/poiseuille.mesh       # 几何、拓扑尺寸、patch 名称与角色
├── fields/initial/U.field     # 初值与 U 的边界条件
├── fields/initial/p.field     # 初值与 p 的边界条件
├── physics/simple.bs          # 密度、黏度等物性
├── numerics/methods.bs        # 算子默认格式及可选的 Field 覆盖
├── numerics/solution.bs       # SIMPLE/线性求解控制
├── control.bs                 # 时间区间与步长
├── output.bs                  # 结果目录与时刻名
└── results/<time>/rank-0000/  # 运行生成，不纳入 Git
```

每个 MPI rank 仅写出 owned cell 的 `U.csv`、`p.csv` 与 `metadata.bs`；ghost cell
不会写出。`babelsim-post` 按 global ID 检查完整性并合并为原始六面体网格的 VTK
XML `.vtu` 或 Tecplot `FEBRICK` 文件；`-time all -format vtk` 还会产生
ParaView 可直接打开的 `post/series.pvd`。

库代码需要从已存在的全局网格分区时仍可使用 `decompose()`；启动器和文件型并行程序
应使用 `readDistributedMesh(path, parallel)`。该接口在 rank 0 读取网格、广播尺寸和
patch 元数据、点对点发送局部顶点；接收方只持有本地 owned+ghost 几何和索引映射。

## 测试与验证

```bash
make test                 # 几何、算子、case/field IO、Heat/标量输运/SIMPLE、通用输出
make test-workflow        # 新 Solver 单函数开发、双场耦合、时间序列、真实 ParaView 读取
make test-external        # 仓库外 Solver 构建、1/2/4 进程、负向 API 编译、无 MPI 结果读取器
make test-mpi             # MPI 网格、halo、算子、线性求解、SIMPLE 与标量输运
make test-mpi-heat        # 1/2 rank 热传导场比较
make test-mpi-poiseuille  # 1/2/4 rank 案例启动器、结果比较、后处理
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
Case 的 validate 只校验，start/loop 才关闭声明；派生 cell 场可通过 output(field) 选择输出。
公开 math 统一为整场同步契约，结果读取头和实现均不再要求 MPI。

默认瞬态结果按 `output.bs` 中 `writeInterval` 保存（省略时每步写出），最终时刻总会保存。
`-time mpi4/all` 后处理命名运行的完整序列；ParaView 打开对应的 `post/mpi4/series.pvd`。
案例名称现为 `heat/simple/transport`，线性配置统一为 `scalarSolver/vectorSolver`。
所有 Case 的 `solution.bs` 必须同时填写这两项，缺项报错，不使用隐式默认选择。

详细设计与数据结构见 [架构说明](docs/architecture.md)，新增物理模型/求解器的流程
见 [Solver 开发指南](docs/solver-development.md)，数学表达规则见
[eqn/math 说明](docs/eqn-math.md)，两个内置 Solver 的对照说明见
[热传导](docs/heat-solver.md)、[SIMPLE](docs/simple-solver.md)，案例组织见
[Case 结构](docs/case-structure.md)，两种 Solver 组织模式见
[Solver Programming Model](docs/solver-programming-model.md)，数值验证结果见
[验证说明](docs/validation.md)；多 Reynolds 数、网格无关性、格式与 MPI 对照见
[Ghia 方腔验证报告](docs/reports/cavity-ghia-validation.md)及其
[PDF 版本](docs/reports/cavity-ghia-validation.pdf)。
