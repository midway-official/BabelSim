# BabelSim 架构

## 目标与边界

BabelSim 是面向有限体积法（FVM）的轻量 CFD/PDE 框架。它借鉴 OpenFOAM 的一个核心
思想：求解器源码应首先表达物理方程和数值算法，而不是 MPI 通信、稀疏矩阵或内存布局。
它不复制 OpenFOAM 的运行时注册表、复杂模板层级和历史兼容结构。

当前框架以三维结构化六面体网格为唯一网格模型；二维是 `nz=1` 的退化情形。所有几何、
Field、算子和并行路径均复用同一实现。

```text
Case / 应用
    │  读取案例、选择求解器、并行结果写出
    ▼
Solver / Algorithm
    │  Heat：时间循环与热方程
    │  SIMPLE：UEqn、面动量插值、pEqn、修正、连续性
    ▼
Public Solver API
    │  Mesh、Field、Boundary、RunTime、fvm、fvc、solve
    ▼
FVM Backend（内部）
    │  离散方程、非正交修正、装配、线性后端
    ▼
Runtime（内部实现）
    │  FieldHistory、halo、全局归约、串行/分布式选择
    ▼
MPI / 内存 / I/O 基础设施
```

依赖方向只能向下。特别地：Physics/Solver 不依赖 `ParallelContext`、`HaloExchange`、
`SparseAssembly`、LDU、Eigen 或 MPI；代数与运行时不知道温度、压力、速度或 SIMPLE。

## OpenFOAM 设计考察

本次重构先对 OpenFOAM 的 `laplacianFoam`、`scalarTransportFoam` 和 `simpleFoam` 做了源码级
对照。标量输运 Solver 用 `runTime.loop()`、`fvm::ddt/div/laplacian` 和 `solve` 表达方程，
把网格、Field、边界、离散格式和线性控制留给 Case/框架；`simpleFoam` 则把
`createFields`、`UEqn.H`、`pEqn.H` 和 SIMPLE 外循环分开，压力方程中显式保留 `rAU`、
预测通量、非正交修正和速度/通量修正。OpenFOAM 的核心经验是“代码顺序对应数学顺序”，
而不是必须复制它的模板和对象注册机制。

BabelSim 采用相同的层次和表达习惯：`fvm` 是隐式方程贡献，`fvc` 是显式场计算；通用 FVM
提供基础算子，SIMPLE 仅保留 `MomentumInterpolation` 和 `PressureCorrection` 两个有独立
数值语义的专用组件。相关对照和原始资料见 [fvm/fvc 说明](fvm-fvc.md)、[热传导 Solver](heat-solver.md)
和 [SIMPLE Solver](simple-solver.md)。

## 当前模块审计

下表按真实 `#include`、构造调用和数据流整理，而不是按目录名称推断。依赖箭头均指向更底层
模块；Physics 只应使用 Public 列中的头文件。

| 模块 | 实际职责与层级 | 主要依赖 | 对外接口 | 内部接口/问题与动作 |
| --- | --- | --- | --- | --- |
| `src/apps` | 案例启动器、MPI 生命周期、结果输出（应用层） | `case`、IO、Physics、并行 | `babelsim-solve`、`babelsim-post` | 可以知道执行细节；不把这些细节传入 Physics |
| `src/io`、`case.h` | 案例字典、网格/Field/物性读取（案例层） | Mesh、Field、配置解析 | `readCase`、读取器 | 不参与离散；新增 Solver 增加小型读取器，不建 Manager |
| `src/physics`、`thermal.h`、`incompressible.h` | 热传导和 SIMPLE 物理/算法（物理层） | Mesh、Field、`fvm/fvc`、RunTime | `solveTransientHeat`、`SimpleSolver` | 不得依赖代数、并行或存储入口；已移除旧泄漏 |
| `fvm.h`、`fvc.h`、`methods.h` | 轻量数学表达和 FVM 方法（离散接口层） | Field、网格几何 | `fvm::*`、`fvc::*`、Methods | 不创建 LDU/通信；只保存描述符，求解时解释 |
| `src/discretization` | 梯度、插值、扩散、对流和方程装配（离散实现层） | Mesh、Field、内部离散方程 | Runtime 调用的算子核 | 可使用连续数组和非正交缓存，不应知道具体 Physics |
| `src/runtime`、`runtime.h` | 方程解释、Field 同步、时间历史、全局判据（执行数学层） | 离散、代数、Parallel | `RunTime::solve/evaluate/loop` | PImpl 隐藏后端；统一管理对象生命周期和 MPI 选择 |
| `src/algebra`、`assembly.h`、`linear_solver.h` | LDU/稀疏装配、Krylov 和预条件器（线性代数层） | Eigen、DiscreteEquation、Parallel | 仅框架级旧/内部接口 | 不知道温度、速度或 SIMPLE；不作为 Physics API |
| `src/parallel`、`parallel.h` | owned/ghost 映射、Halo、全局归约、并行写出（运行时基础设施） | MPI、Mesh、Field | 仅启动器/Runtime 使用 | 校验 MPI 生命周期和返回码；不向 Physics 暴露 |
| `src/core`、`mesh.h`、`field.h` | 三维结构化 Mesh/Field 数据和不变量（核心对象层） | 标量/向量值类型 | Mesh、Field、Boundary | 连续存储和容量由构造决定；内部受控修改，禁止任意 resize |

这一审计也解释了为什么没有新增 `PhysicsManager`、`EquationManager` 或 `SolverFactory`：当前
启动器只有明确的 Solver 分支，足以满足案例选择，同时避免无数学语义的注册层和循环依赖。

## 重构前后的关键变化

此前公开的 `Equation<T>` 实际是 LDU 系数容器，`SimpleSolver` 还公开保存了并行上下文、
halo、装配器、Eigen 向量和线性求解器。这使“方程”同时表示数学对象、离散对象和执行对象，
也迫使物理代码理解底层生命周期。

当前采用两个明确层次：

| 层次 | 代表对象 | 含义 | 不应知道的内容 |
| --- | --- | --- | --- |
| 数学表达 | `ScalarExpression`、`VectorExpression`、`ScalarEquationDefinition` | `lhs == rhs` 的轻量描述 | LDU、CSR、MPI、矩阵大小 |
| 离散方程 | Runtime 内部 `DiscreteEquation`/LDU | 某个 FVM 方法产生的局部代数系统 | 压力、温度、SIMPLE |

表达式只保存常数和 Field 引用；`operator+`、`operator-`、`operator==` 不创建大型 Field、
LDU 或 MPI 缓冲。仅在 `solve(...)` 时才同步必要 halo、执行非正交修正、装配并调用线性后端。

## Public Solver API

新的 PDE/Solver 开发者应主要使用下面的头文件和概念：

```text
babelsim/mesh.h          Mesh、patch、几何和拓扑
babelsim/field.h         ScalarField、VectorField、TensorField、边界条件
babelsim/fvc.h           显式场运算描述
babelsim/fvm.h           隐式方程项与 EquationDefinition
babelsim/runtime.h       RunTime、solve()、时间循环与收敛诊断
babelsim/thermal.h       热传导物性和最小 Solver
babelsim/incompressible.h  不可压缩物性、场和 SIMPLE 算法
```

典型隐式方程为：

```cpp
solve(
    runTime,
    fvm::ddt(rhoCp, T)
        == fvm::laplacian(k, T) + fvm::source(Q));
```

其中 `rhoCp` 与 `k` 可以是常数，也可以是 cell Field；Runtime 会为离散所需的面系数
完成同步和插值。`solve` 不要求调用者持有矩阵或线性求解器。

边界条件属于 Field：

```cpp
T.boundary("hot") = fixedValue(500.0);
T.boundary("wall") = zeroGradient();
U.boundary("side") = symmetry();
```

## 内部框架 API

下列对象属于 FVM、代数或运行时实现，不能成为 Physics Solver 的日常接口：

```text
Equation/LDU、SparseAssembly、PreparedLinearSolver、DistributedLinearSolver
ParallelContext、HaloExchange、MPI 通信器、owned/ghost 映射、Eigen 向量
```

它们仍然存在，因为需要提供高效分布式执行：局部 Mesh 仅保存 owned+ghost 单元，Field
连续存储于对应局部实体上；Runtime 在算子和矩阵向量乘前选择正确的 halo 同步宽度，在线性
残差、范数和 SIMPLE 外迭代判据处执行全局归约。它们不会向 Heat 或 SIMPLE 主算法泄漏。

`Field` 的容量仅由 `(Mesh, FieldLocation)` 构造时确定，外部不能 `resize`。Framework 内部
保留连续数据访问入口以服务 halo/计算热点；新 Solver 不应使用 `data()`、`mutableData()`、
`values()` 或 cell 循环实现 PDE。

## Mesh、Field 与边界

`Mesh` 描述空间的位置和连通性，而不携带物理含义。它预计算并连续保存：

```text
cell/face/vertex，owner/neighbour，cell-face/邻居关系，patch
cell/face 中心，体积及逆体积，面积向量及单位法向
正交系数、非正交向量、偏斜量、owner 插值权重
局部 owned/ghost 与 global-ID 映射
```

`Field<T>` 描述“空间中有什么”。`ScalarField`、`VectorField`、`TensorField` 通过同一连续
存储实现，位置由 `FieldLocation::Cell/Face/Vertex` 指定。压力 `p`、温度 `T`、速度 `U`
不是不同 Field 类型。

cell Field 绑定 patch 边界条件。已实现 `fixedValue/Dirichlet`、
`fixedGradient/Neumann`、`zeroGradient`、`inletOutlet`、`symmetry/Mirror`；`Wall`、`Inlet`
和 `Outlet` 是 Mesh patch 的物理角色，具体数学条件仍由 Field 指定。

## fvm 与 fvc

`fvm` 表示**隐式进入待求解方程的项**，例如 `fvm::ddt`、`fvm::div`、
`fvm::laplacian`、`fvm::source`。`fvc` 表示**直接计算为 Field 的显式量**，例如
`fvc::grad`、`fvc::div`、`fvc::flux`、`fvc::interpolate`、`fvc::laplacian`。

这个分离对应 OpenFOAM 的设计思路，但 BabelSim 的表达式描述是简单的引用列表而不是
模板化矩阵系统。详细规则见 [fvm-fvc.md](fvm-fvc.md)。

## SIMPLE 的位置

SIMPLE 不是通用 FVM 的一部分。通用层只定义梯度、散度、对流、扩散、插值、通量与时间项。
不可压缩算法层保留两个具有独立数值语义的组件：

- `MomentumInterpolation`：同位网格的 Rhie--Chow 风格面动量插值；
- `PressureCorrection`：由 `rAU`、预测通量与连续性形成压力修正，并一致修正 `p/U/phi`。

二者内部使用同一个 `RunTime`、`fvm/fvc` 与 Mesh/Field，不创建第二套矩阵、通信或数据布局。
其工作场在 SIMPLE 创建一次并跨外迭代复用，避免重复完整 Field 分配。

## 当前扩展边界

当前离散后端是 FVM。未来 FDM/FEM 可在“数学表达 → 离散方程”边界新增后端：复用 Mesh、
Field、Boundary、Expression、RunTime 的高层语义，并以各自的离散系统实现内部 assembly/solve。
在该能力真正实现前，BabelSim 不声称已经支持 FDM/FEM。

新增热传导、对流扩散或标量输运通常只需要新建一个物理源文件和 Case 字典，不需要改 MPI、
LDU、halo 或并行 I/O。新增压力速度耦合算法则应位于 `physics/incompressible`，并仅在具有
稳定独立数值意义时增加专用组件。
