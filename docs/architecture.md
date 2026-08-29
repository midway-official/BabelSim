# BabelSim 架构与数据结构

## 1. 设计目标与边界

BabelSim 的目标不是维护一个只支持二维 Navier-Stokes 的程序，而是提供可组合的
有限体积计算基础。当前第一个物理算法是不可压 Navier-Stokes/SIMPLE，但 Mesh、
Field、通用算子、装配、线性代数、MPI 通信和结果写出都不包含压力、速度或 SIMPLE
的专有语义。

TaihoCFD 提供了经过回归验证的 SIMPLE 状态转移、迎风对流、正交扩散极限、
Rhie-Chow 插值、压力修正和 MPI 分区思路。BabelSim 保留这些数值算法，重新组织
其数据归属：网格不保存 `U/p`，场不实现离散，线性求解器不认识 CFD 变量。

不采用 `Manager`、多层 Factory、运行时算子注册表或深继承。核心计算路径使用
普通结构体、连续数组、整数索引和直接函数调用。

## 2. 运行链路与目录

求解和后处理的链路如下：

```text
babelsim-solve -case cases/poiseuille
  │
  ├─ 读取 case.bs
  ├─ rank 0 读取全局 Mesh，按 x 向只分发各 rank 的 owned + ghost Mesh
  ├─ 每个 rank 直接在本地 Mesh 上创建初始 Field（不复制全局 Field）
  ├─ 按 solver 名称分派物理算法
  ├─ SimpleSolver 组合通用算子、方程、装配和线性求解器
  └─ 每个 rank 写出 owned field

babelsim-post -case cases/poiseuille -format vtk tecplot
  │
  ├─ 读取全局 Mesh
  ├─ 读取所有 rank 的 metadata.bs 与 Field CSV
  ├─ 按 global_id 检查无重复、无缺失
  └─ 输出 VTK UNSTRUCTURED_GRID / Tecplot FEBRICK
```

生产 C++ 实现全部位于 `src/`：

```text
src/
├── core/                    # 网格几何与拓扑
├── io/                      # case、mesh、field、result 读写
├── discretization/          # 通用算子与稀疏装配
├── algebra/                 # 串行/分布式线性求解
├── parallel/                # 并行上下文、halo、并行写出
├── physics/incompressible/  # SIMPLE 与两个 CFD 专用算子
└── apps/                    # 通用启动器与后处理器
```

公共 API 位于 `include/babelsim/`；测试位于 `tests/`；案例位于 `cases/`；验证与
比较工具位于 `tools/`。

## 3. 抽象层及职责

| 层 | 主要类型/文件 | 职责 | 不应承担的职责 |
|---|---|---|---|
| Mesh | `Mesh`、`mesh.h` | 描述空间、拓扑、几何和局部所有权 | 保存物理变量或方程 |
| Field | `Field<T>`、`field.h` | 保存指定位置的数据与该场的边界条件 | 选择离散格式或组装矩阵 |
| Operator | `operators.h/.cpp` | 对 Field 或 Equation 执行通用数学操作 | 读取 case、隐藏 MPI 全局场 |
| Method | `Methods`、`methods.h` | 选择梯度、对流、扩散、时间离散策略 | 承载物理变量 |
| Equation | `Equation<T>` | 保存面耦合 LDU 系数和源项 | 生成 Eigen 矩阵或求解 |
| Assembly | `SparseAssembly` | 将 owned 行映射到稀疏矩阵与 RHS | 解释压力、速度等物理意义 |
| Linear Solver | `PreparedLinearSolver`、`DistributedLinearSolver` | 求解 `Ax=b` 与报告残差 | 处理边界条件或 CFD 流程 |
| Physics | `SimpleSolver` | 组织方程与算法的调用次序 | 读取文件、解析命令行、写 VTK |
| Apps / IO | `babelsim-solve`、`babelsim-post` | 启动、配置读取、结果交换格式 | 把数值算法塞入 IO |

依赖只从上层组合下层，禁止反向依赖。例如 `operators.cpp` 不包含 case 读取器，
`linear_solver.cpp` 不包含不可压求解器，`mesh.cpp` 不包含 Field。

## 4. Mesh：统一三维、预计算几何与 MPI 映射

`Mesh` 是结构化六面体网格。逻辑尺寸为 `(nx, ny, nz)`；二维算例使用 `nz=1`，
前后面通常配置为 symmetry。所有二维和三维算子均使用同一套六个面连接关系，
没有 `if (dimension == 2)` 的分支实现。

### 4.1 连续数据布局

`Mesh` 采用索引数组而非指针图，主要数组如下：

| 数据 | 类型 | 含义 |
|---|---|---|
| `vertices` | `MeshStorage<Vec3>` | 顶点坐标，按结构化逻辑顺序连续存储 |
| `cell_centres`、`cell_volumes`、`cell_inverse_volumes` | `MeshStorage<Vec3/double>` | 单元几何及逆体积缓存 |
| `cell_faces[6]` | `MeshStorage<std::array<Index,6>>` | 单元六个面的局部索引 |
| `cell_neighbours[6]` | 同上 | 六个相邻单元；边界为 `invalid_index` |
| `face_vertices[4]` | `MeshStorage<std::array<Index,4>>` | 四边形面的顶点索引 |
| `face_owner`、`face_neighbour` | `MeshStorage<Index>` | 面的 owner/neighbour 单元 |
| `face_patch`、`patches` | `MeshStorage` | 边界面的 patch 归属 |
| `face_centres`、`face_area_vectors`、`face_normals`、`face_areas` | `MeshStorage` | 面中心、面积向量、单位法向、面积 |
| `face_orthogonal_coefficients` | `MeshStorage<double>` | 扩散隐式正交部分的系数 |
| `face_non_orthogonal`、`face_skewness` | `MeshStorage<Vec3>` | 非正交与偏斜显式修正 |
| `face_owner_weights` | `MeshStorage<double>` | owner 到面上的线性插值权重 |

这些量在 `Mesh::structured()` 或 `Mesh::cartesian()` 中只计算一次。热点算子通过
内联只读访问器或数组索引读取，不重新计算几何，也不进行虚函数分派。
数组使用 `MeshStorage<T>` 保留连续内存；`resize/clear/push_back` 仅能由 Mesh 的
构造和分区成员函数调用，避免外部改变索引布局。`Mesh::validate()` 在装配、halo
和分布式读取边界检查拓扑、几何、所有权及 global ID 一致性。

内部面面积向量从 owner 指向 neighbour；边界面面积向量从 owner 指向域外。扩散将
面积向量拆为正交隐式部分与非正交显式部分：

```text
Sf = ((Sf · d) / |d|²) d + k_nonorth
```

其中第一项形成矩阵系数，`k_nonorth · grad(phi)` 作为修正项。偏斜量用于把
owner 与 neighbour 连线上的插值值重构到真实面中心。

### 4.2 局部网格与所有权

当前 MPI 分区沿 x 方向。`decompose()` 适用于调用方已经持有全局 Mesh 的库代码；
应用启动器使用 `readDistributedMesh()`：仅 rank 0 调用 `readMeshFile()`，随后广播
尺寸/patch 元数据并逐 rank 发送局部顶点。接收 rank 直接重建本地几何并建立完整的
owned/ghost 映射，整个读取过程不会在非零 rank 保存全局 Mesh。局部网格保留两层
ghost cell；两层是为了让非正交扩散在第一层 ghost 单元使用已重构的梯度。

局部 `Mesh` 还保存：

- `owned_cells`：本 rank 真正拥有的局部 cell 索引；
- `owned_faces`：至少连接一个 owned cell 的局部面索引；代数装配、通量更新和压力修正只遍历这些面，避免扫描 ghost-only 面；
- `cell_owned_indices`：局部 cell 到 owned 行号的映射，ghost 为 `invalid_index`；
- `cell_global_ids`：局部 cell 到全局 cell ID 的映射；
- `global_dimensions`、`global_i_offset`、`owned_i_begin/end`、`ghost_layers`：
  局部布局元数据。

物理边界仍是普通 patch；分区接口标为 `PatchKind::Processor`。现阶段不实现任意
非结构化图分区、xyz 多向切分或 MPI-IO，这些属于不改变核心接口的后续优化。

## 5. Field 与边界条件

`Field<T>` 是模板容器，当前常用别名为：

```cpp
using ScalarField = Field<double>;
using VectorField = Field<Vec3>;
using TensorField = Field<Tensor3>;
```

Field 保存 `const Mesh*`、位置 `FieldLocation`、名称、连续 `std::vector<T> values_`
和 cell field 的 patch 边界条件。`values_` 的长度由 `(Mesh, FieldLocation)` 唯一
决定；公开接口不提供 `resize`，`mutableData()` 仅用于框架内部 halo 原地写入，
`validateStorage()` 在通信、算子和输出入口检查不变量。`Vec3` 与 `Tensor3` 分别是
紧凑的三/九个 `double` 的小对象，因此一个 Field 内按实体连续，当前是 AoS 值布局，
没有隐含指针或逐单元分配。

位置可以是 Cell、Face、Vertex。当前不可压求解器使用 cell-centered（单元中心）`U/p` 和
face-centered `phi`；通用 halo 已支持 cell/face scalar、vector、tensor。标准化结果
写出目前针对 cell field，因为只有 cell 拥有完整的全局 ID 映射。

MPI 启动时每个 rank 仅构造局部 Field。当前 `.field` 文件是 uniform 初值，因此各
rank 可独立解析同一个小文件；processor patch 自动采用 ZeroGradient，接口数据由
`HaloExchange` 同步。cell ghost 由相邻 owned cell 覆盖，重复的 face field 由低 rank
owner 单向发布到高 rank，保证同步是幂等的。若以后加入非 uniform 场文件，应只按 `cell_global_ids` 分发或
读取本地 owned 数据，不能重新创建全局 Field。

边界条件属于 Field 而不是 Solver：FixedValue/Dirichlet、FixedGradient/Neumann、
ZeroGradient、InletOutlet、Symmetry/Mirror。向量 symmetry 去除法向分量；标量
symmetry 等价于零法向梯度。patch 的 `Wall/Inlet/Outlet/Symmetry` 仅是几何角色，
每个物理场仍自行选择数学边界条件。

## 6. 通用算子、离散方法、方程与装配

`Methods` 用小枚举表达离散策略，而非为每种格式建立继承树：

- `InterpolationMethod::Linear`；
- `InterpolationMethod::Corrected`（交点线性值加偏斜重构）；
- `GradientMethod::GreenGauss`、`LeastSquares`；
- `ConvectionMethod::Upwind`、`Central`；
- `DiffusionMethod::Orthogonal`、`Corrected`、`LimitedCorrected`；
- `TimeMethod::Steady`、`Euler`、`BDF2`。

通用算子包括 `gradient`、`interpolate`、`reconstruct`、`flux`、`divergence`、
`laplacian`、`addConvection`、`addDiffusion` 和 `addTimeDerivative`。它们只读取
局部 Mesh/Field/Method；调用方在读取 ghost 值之前显式调用 `HaloExchange`：

```cpp
HaloExchange halo(mesh, parallel);
halo.exchange(phi);                          // 同步 cell 输入场
gradient(phi, grad_phi, methods.gradient);   // 在局部网格上计算
halo.exchange(face_flux);                     // 后续若需使用接口面通量
addConvection(equation, face_flux, phi, methods.convection);
addDiffusion(equation, gamma, phi,
             methods.gradient, methods.diffusion);
```

`Corrected` 面值先在 owner 与 neighbour 连线和真实面中心的交点处做线性插值，再用
交点梯度乘以 `face_skewness` 重构到真实面中心。扩散和压力修正使用同一个
`integratedNormalGradient()`：面积向量先分解为正交隐式部分和非正交显式部分，
`LimitedCorrected` 将显式交叉扩散限制在正交通量幅值以内，以提高大角度网格的稳定性。
Green--Gauss 梯度也执行一次显式偏斜面值修正；Least--Squares 梯度直接由相邻
单元中心位移构造三维法方程。边界重构会施加固定值、固定法向梯度、对称投影和
入口出口流向条件，因此这些接口对 `nz=1` 和完整三维网格使用同一实现。

`Equation<T>` 对每个局部 cell/face 保存 `diagonal`、`upper`、`lower` 和 `source`。
这样通用算子可以同时处理 owned 与 ghost 邻接关系。`SparseAssembly` 只生成 owned
行，并在构造时预计算 Eigen 稀疏矩阵的非零位置；每次 `update()` 只改值，不重复
triplet 分配、排序或插入。

串行 `PreparedLinearSolver` 复用稀疏结构和预条件器；MPI `DistributedLinearSolver`
执行 owned block 的矩阵乘法，只交换矩阵接口所需的第一层 ghost 值，再加跨分区 face
系数。跨分区系数仅缓存接口面快照，不在每个外迭代复制完整 LDU 数组。所有
点积、范数和连续性统计使用 `MPI_Allreduce`。SIMPLE 还对健康状态、内层线性收敛
状态和外迭代物理收敛条件做全局归约，所有 rank 因而在同一个外迭代上继续或停止；
局部预条件器的 `MaxIterations` 不会再导致某个分区单独多做外迭代。

## 7. 不可压算法与专用算子

`SimpleSolver` 只负责以下调用顺序：

```text
组装动量方程
→ 求解三个速度分量
→ MomentumInterpolation::apply
→ PressureCorrection::assemble
→ 求解压力修正方程
→ PressureCorrection::apply
→ 连续性与场变化收敛判断
```

不将每一步拆成一次性小类，只保留两个具有独立 CFD 数值含义的专用组件：

- `MomentumInterpolation`：根据同位网格的 `U/p`、压力梯度、动量 mobility 和
  面几何执行 Rhie-Chow 型插值，生成压力稳定的面通量；
- `PressureCorrection`：组装压力修正 LDU 方程，并同步执行压力、cell velocity 与
  face flux 修正。

两者仍只依赖 Mesh、Field、Equation、ParallelContext 和通用边界/几何接口；不会
重新定义底层数据结构。压力出口和封闭域压力参考点处理保留在该算法组件中。

非正交扩散采用显式校正迭代：矩阵只保留紧凑的正交 LDU 系数，`SimpleSolver` 在
每个压力修正外迭代中先解正交系统，再按 `nonOrthogonalCorrections` 次数更新
压力修正梯度和右端项。矩阵、装配和通信接口不因该循环改变；因此专用算法仍然
可以替换为其他压力--速度耦合算法。

## 8. 原生文件与结果格式

`case.bs` 只描述入口和相对路径：

```text
solver simpleFoam
mesh mesh/poiseuille.mesh
fields fields/initial
physics physics/incompressible.bs
numerics numerics/simple.bs
output output.bs
```

`.mesh` 记录结构化尺寸、Cartesian bounds 或显式顶点，以及六个逻辑边的 patch。
`.field` 使用花括号，记录 Field 类型、位置、uniform 初值及以 patch 名称索引的
边界条件。物性与离散控制独立存放，因而新增物理模型可以复用同一网格和输出布局。

结果目录格式为：

```text
results/<time>/rank-0000/
├── U.csv
├── p.csv
└── metadata.bs
```

CSV 含 `global_id,x,y,z,value0,...`。后处理器不会访问求解器内存，只依靠原始 mesh、
metadata 和 rank CSV 恢复全局 cell 顺序并检查完整性。

## 9. 性能原则与当前范围

- 几何（含逆体积和单位法向）、稀疏结构、owned-face 列表、halo 通信平面、迭代向量和线性求解 workspace 均复用；
- Operator 热点无虚函数、无全局 Field、无不必要的大型临时对象；
- 正交网格跳过等价的偏斜重构 pass；分布式 Krylov 矩阵乘只交换第一层 ghost，场和梯度同步仍使用完整 halo；
- MPI 全局归约保持生命周期检查，但不在每次归约重复查询 rank/size；SIMPLE 的三个状态标志合并为一次整型归约；
- 结果阶段只写 owned cell，避免 ghost 重复和串行重组；
- case 解析仅发生在启动阶段，不进入迭代循环。

当前已验证结构化单节点/MPI 基础、三维扭曲非正交腔体和分区交界面上的全部主要
通用算子。下一阶段可增加 Advection-Diffusion、Heat 等
物理模型；它们应新增 `src/physics/<模型>/` 和必要的 IO/启动器分派，而不修改
Mesh、Field、HaloExchange、SparseAssembly、LinearSolver 或 `babelsim-post` 的核心。

## 10. 数值依据

实现遵循有限体积法中“正交隐式项 + 非正交显式项”的标准分解，以及偏斜面中心
梯度重构。可复核的主要资料包括：

- [OpenFOAM 数值格式说明：corrected 与 limited corrected snGrad](https://doc.cfd.direct/openfoam/user-guide-v8/fvschemes?s=2025)；
- [OpenFOAM correctedSnGrad 源码](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/finiteVolume/snGradSchemes/correctedSnGrad/correctedSnGrad.C)；
- [OpenFOAM snGradScheme 基类实现](https://github.com/OpenFOAM/OpenFOAM-dev/blob/master/src/finiteVolume/finiteVolume/snGradSchemes/snGradScheme/snGradScheme.C)；
- [Rhie--Chow 动量插值的数值研究](https://www.sciencedirect.com/science/article/pii/S0377042716301649)；
- [广义 Rhie--Chow 插值研究](https://www.sciencedirect.com/science/article/pii/S0021999113007523)；
- [基于网格偏斜度的梯度修正研究](https://doi.org/10.1063/5.0246823)。

这些资料用于确定离散公式和稳定性边界；BabelSim 当前仍限定为结构化六面体，
并不宣称已经覆盖任意非结构化、多面体或高于约 70 度非正交角的无条件稳定性。
