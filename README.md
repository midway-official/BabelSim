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
fvm/fvc   分别表达隐式离散项与显式场运算
RunTime   管理时间、离散、求解、并行同步与收敛
Physics   如何组合方程与算子解决具体物理问题
```

`Equation` 的数学表达与内部离散 LDU 系统明确分离。Heat 与 SIMPLE 的 Physics 源码不接触
MPI、halo、CSR/LDU、Eigen 或 Field 底层存储；这些由 Runtime 自动处理。

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
- 不可压 SIMPLE；仅保留 `MomentumInterpolation` 与 `PressureCorrection` 两个
  具有独立数值语义的 CFD 专用算子；
- OpenFOAM 风格的轻量 `fvm/fvc`：数学表达式只在 `solve()` 时离散，不复制大矩阵；
- `heatFoam` 瞬态热传导 Solver，支持常数或 Field 系数；
- 原生 case/mesh/field 文件、通用并行结果写出与独立 VTK/Tecplot 后处理。

## 构建与运行

依赖：C++17 编译器、Eigen 3、MPI-3 实现和 GNU Make。默认编译器为 `mpic++`，
保证串行测试与 MPI 程序使用相同 ABI。

```bash
make -j

# 串行腔体
build/babelsim-solve -case cases/cavity

# MPI 通道流
mpirun -np 4 build/babelsim-solve -case cases/poiseuille

# 瞬态热传导：串行与 MPI 使用同一 Solver 源码
mpirun -np 2 build/babelsim-solve -case cases/heat

# 独立读取网格和并行结果，输出 ParaView/Tecplot 文件
build/babelsim-post -case cases/poiseuille -format vtk tecplot
```

`-time <名称>` 可为同一案例保存多个结果时刻。例如：

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
├── physics/incompressible.bs  # 密度、黏度等物性
├── numerics/simple.bs         # 算子方法、SIMPLE 与线性求解控制
├── output.bs                  # 结果目录与时刻名
└── results/<time>/rank-0000/  # 运行生成，不纳入 Git
```

每个 MPI rank 仅写出 owned cell 的 `U.csv`、`p.csv` 与 `metadata.bs`；ghost cell
不会写出。`babelsim-post` 按 global ID 检查完整性并合并为原始六面体网格的 VTK
`UNSTRUCTURED_GRID` 或 Tecplot `FEBRICK` 文件。

库代码需要从已存在的全局网格分区时仍可使用 `decompose()`；启动器和文件型并行程序
应使用 `readDistributedMesh(path, parallel)`。该接口在 rank 0 读取网格、广播尺寸和
patch 元数据、点对点发送局部顶点；接收方只持有本地 owned+ghost 几何和索引映射。

## 测试与验证

```bash
make test                 # 几何、算子、case/field IO、专用算子、通用输出
make test-mpi             # MPI 网格、halo、算子、装配、线性求解与 SIMPLE
make test-mpi-heat        # 1/2 rank 热传导场比较
make test-mpi-poiseuille  # 1/2/4 rank 案例启动器、结果比较、后处理
make validate-cavity      # 收敛的 Re=100 Ghia 腔体比较
make validate-poiseuille  # 收敛的 Poiseuille 解析解比较
```

非正交离散的标准入口是 `interpolation corrected`、`gradient leastSquares`、
`diffusion corrected`；网格角度较大时可改用 `diffusion limitedCorrected`，并通过
`nonOrthogonalCorrections` 设置压力修正右端项的显式迭代次数。串行测试还包含
三维扭曲腔体，MPI 测试包含分区界面上的修正通用算子和两个专用 CFD 算子。

详细设计与数据结构见 [架构说明](docs/architecture.md)，新增物理模型/求解器的流程
见 [Solver 开发指南](docs/solver-development.md)，数学表达规则见
[fvm/fvc 说明](docs/fvm-fvc.md)，两个内置 Solver 的对照说明见
[热传导](docs/heat-solver.md)、[SIMPLE](docs/simple-solver.md)，案例组织见
[Case 结构](docs/case-structure.md)，数值验证结果见 [验证说明](docs/validation.md)。
