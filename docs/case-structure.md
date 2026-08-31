# BabelSim Case 结构与启动器

## 设计原则

一个 Case 描述“要计算的问题”；一个 Solver 描述“如何用物理方程与算法求解”。两者不能
互相硬编码。BabelSim 借鉴 OpenFOAM 的 Case/solver 分离和 `fvSchemes`/`fvSolution`
职责分离，但只保留适合当前规模的七个文件入口。

```text
cases/<case-name>/
├── case.bs
├── mesh/<name>.mesh
├── fields/initial/<field>.field
├── physics/<model>.bs
├── numerics/methods.bs
├── numerics/solution.bs
├── control.bs
└── output.bs
```

运行后生成、且不进入 Git 的内容为：

```text
results/<time>/rank-0000/...   每个 rank 的 owned 数据
post/<time>.vtk                ParaView 文件
post/<time>.dat                Tecplot 文件
```

## case.bs

`case.bs` 是通用启动器读取的唯一固定入口：

```text
solver simpleFoam
mesh mesh/cavity.mesh
fields fields/initial
physics physics/simple.bs
methods numerics/methods.bs
solution numerics/solution.bs
control control.bs
output output.bs
```

启动命令始终相同：

```bash
build/babelsim-solve -case cases/cavity
mpirun -np 4 build/babelsim-solve -case cases/cavity -time mpi4
```

启动器读取 `solver` 后选择对应的 Case 运行入口。启动器只负责参数和 MPI 生命周期；每个
Physics 目录中的运行入口负责其 Case→Field→RunTime→输出流程，因此通用应用不需要知道
`T/U/p/C` 等具体 Field。Physics Solver 本身不需要知道进程数或通信器。

当前内置选择为：

| `solver` | 物理/算法 | 必需初值 |
| --- | --- | --- |
| `heatFoam` | 瞬态热传导 | `T.field` |
| `simpleFoam` | 稳态层流不可压缩 SIMPLE | `U.field`、`p.field` |
| `transportFoam` | 瞬态标量对流-扩散 | `C.field`、`U.field` |

新增 Solver 只需在启动器增加一个明确分支及其 Case 读取器；不要建立无边界的通用注册表或
Factory 层。

## 原生网格文件

原生 `.mesh` 文件描述结构化三维六面体网格、范围或显式顶点以及 patch：

```text
BABELSIM_MESH 1
dimensions 32 1 1
geometry cartesian
bounds 0 0 0 1 1 1
patch xmin hot wall
patch xmax cold wall
patch ymin side symmetry
patch ymax side symmetry
patch zmin front symmetry
patch zmax back symmetry
end
```

`dimensions nx ny nz` 永远是三维尺寸。`nz=1` 不是单独二维网格，而是退化三维网格。
对于非正交结构化网格可使用显式顶点模式；Mesh 构造时从顶点生成 cell/face、面积向量、
体积、正交系数、非正交修正和偏斜量。

并行时只有 rank 0 解析完整文件，随后按 x 方向生成局部 owned+ghost Mesh 并发送所需几何。
其他 rank 不保存完整全局 Mesh；Field 也只保存相同局部实体上的连续数据。

## 初值和边界 Field

`fields/initial/T.field`、`U.field`、`p.field` 等文件指定初值及各 patch 的数学边界。
示意：

```text
field T cell scalar
internal uniform 300
boundary hot fixedValue 500
boundary cold fixedValue 300
boundary side symmetry
end
```

可用的数学边界包括 `fixedValue`、`fixedGradient`、`zeroGradient`、`inletOutlet`、
`symmetry`/`mirror`。patch 的 `wall/inlet/outlet/symmetry` 角色来自网格；Field 文件决定
该物理量在该 patch 上的具体数学约束。

## physics、methods、solution 与 control

`physics` 只放物性和物理源项。例如热传导：

```text
density 1
heatCapacity 1
conductivity 0.1
source 0
```

`methods.bs` 只放离散格式。例如：

```text
interpolation linear
gradient greenGauss
convection upwind
diffusion orthogonal
time euler
convection C upwind
```

每行既可写默认格式（两个 token），也可写“算子类别、Field 名、格式”（三个 token）。例如
`convection C upwind` 仅覆盖 `div(phi,C)`，`convection U central` 可独立覆盖动量方程。

`control.bs` 只放时间区间和步长：

```text
startTime 0
endTime 0.05
deltaT 0.01
```

`solution.bs` 放线性求解器和算法收敛控制：

```text
scalarSolver bicgstab ilut 1e-14 1e-10 1000
```

`simpleFoam` 的 `solution.bs` 还包括 `maxIterations`、`nonOrthogonalCorrections`、速度/压力
欠松弛、连续性/速度容差以及 `velocitySolver`/`pressureSolver`。物性与数值方式分开，
所以改变黏度不需要改线性求解控制，改变非正交扩散也不需要改 Solver 源码。

## 并行输出与后处理

每个 rank 只写 owned cell/face，对应全局 ID；ghost 数据不会写出。metadata 记录网格和
Field 类型。独立程序合并这些分区文件：

```bash
build/babelsim-post -case cases/poiseuille -time mpi4 -format vtk tecplot

# 自动发现可兼容的 results/<time>，生成 post/series.pvd
build/babelsim-post -case cases/poiseuille -time all -format vtk
```

它会检查分区全局 ID 的完整性，并产生 VTK `UNSTRUCTURED_GRID`（ParaView）与 Tecplot
`FEBRICK`。`all` 会跳过与当前网格不兼容的旧结果并报告原因，并在 VTK 输出时写出 `.pvd`
索引。因此求解器不包含串行聚集、可视化格式或 MPI I/O 细节。

## 与 OpenFOAM 的关系

OpenFOAM 常把初值放在 `0/`、物性放在 `constant/`、离散/线性控制放在 `system/`。
BabelSim 保留“初值、物性、数值控制分离”的思想，但用 `fields/`、`physics/`、`numerics/`
与单个 `case.bs` 保持小规模、明确且易于新用户理解。

相关背景可见 OpenFOAM 的 [Case 文件结构说明](https://www.openfoam.com/documentation/user-guide/2-openfoam-cases/2.1-file-structure-of-openfoam-cases)、
[数值格式说明](https://www.openfoam.com/documentation/user-guide/6-solving/6.2-numerical-schemes) 与
[求解/算法控制说明](https://www.openfoam.com/documentation/user-guide/6-solving/6.3-solution-and-algorithm-control)。
