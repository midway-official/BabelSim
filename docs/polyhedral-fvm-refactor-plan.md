# BabelSim 通用 face-based polyhedral FVM 重构方案

日期：2026-09-22。源码审计基线：`c6c9dfa`。本文保留设计及验收方案；当前实现和实际验证状态见 [2026-09-22-status.md](reports/polyhedral-fvm/2026-09-22-status.md)。
配套执行指令：`polyhedral-fvm-refactor-prompts.md`。

## 1. 目标与硬约束

计算拓扑支持可变 face 顶点数、可变 cell face 数。每个计算 face 有且只有一个 owner，至多一个 neighbour。所有物理离散通过统一计算面累加守恒贡献，不依赖 tetra/hex/prism/polyhedron 类型分支。

用户要求：重构不能损失线性系统求解性能。将其落实为相同代数问题下的性能门禁，保留现有 Eigen/CSR、Krylov、预条件器和 MPI 通信重叠路径；将拓扑处理放在初始化，静态网格迭代中只更新系数和场。不能用更松容差、更少修正、更换预条件器或不收敛的短跑证明性能不变。

区分三个问题：

1. 相同矩阵、分区、初值、求解配置：验收代数后端吞吐和求解时间不退化。
2. 相同旧 hex 网格和离散配置：验收兼容迁移后的几何、矩阵、结果和整条流程性能。
3. 不同 polyhedral 网格：面数、nnz、条件数和迭代数可能改变，比较每 nnz 成本、误差及达成同等精度的时间。不能承诺任意不同网格总时间完全一致。

首轮交付静态、三维、有效闭合的面连接网格，包含 hex、tet、prism、pyramid、真正多边形面/多面体和显式切分的 coarse/fine 界面。支持拓扑上的任意面数，不等于接受自交、非流形或数值不可接受的任意几何。

首轮不扩展到运行时 AMR、ALE、VOF、overset、滑移/非匹配接口搜索、全新分布式 AMG、GUI 或完整 Gmsh/CGNS 导入。这些是后续项目。运行时只保留通用 face topology；旧 hex 算例通过一次性离线迁移工具转换为 v3，求解器不再保留旧 hex 入口。

## 2. 已核实的当前源码事实

| 位置 | 当前事实 | 重构动作 |
|---|---|---|
| `include/babelsim/mesh.h` | cell_vertices[8]、cell_faces[6]、cell_neighbours[6]、face_vertices[4]；已有 face_owner/neighbour | 改为连续变长连接与只读 range；明确 cellCount 不依赖 hex 顶点数组 |
| `src/core/mesh.cpp` | hex_faces 固定模板；quadGeometry；addTetrahedron 使用体积绝对值 | 固定模板留在导入适配层；通用多边形几何、有向体积和质量检查 |
| `src/io/mesh_reader.cpp` | 历史 reader 曾只读 BABELSIM_MESH 2 的 hex 格式 | 只读显式 faces/owner/neighbour 的 v3；v2 仅由离线迁移工具处理 |
| `src/discretization/operators.cpp` | 已有逐 face 算子；least-squares 等遍历 cell_faces | 保留面核；邻接改为变长 range，审计边界及重构权重 |
| `src/discretization/fvm_execution.cpp` | 有 cell_faces 访问 | 迁移 stencil，不改变 Solver 的方程流程语义 |
| `include/babelsim/discrete_equation.h` | diagonal/source 按 cell，upper/lower 按 face | 保留现有 LDU 系数契约 |
| `src/backend/eigen_assembly.cpp` | 图构建后缓存 coefficient positions；update 对非对角位置直接赋值 | 保留缓存；增加多个 face 对同一矩阵位置的汇总，防覆盖 |
| `src/algebra/distributed_solver.cpp` | 已有连续 CSR 及 pattern 复用；buildCoarseMapping 用 global_cells*6 图 | 保留迭代热路径；初始化图改为变长 CSR |
| `src/parallel/parallel_context.cpp` | 广播和局部分区由 8 顶点 hex 重建；quad face 匹配 | 直接分发/抽取原始面拓扑与 global IDs，不再重建 hex |
| `include/babelsim/result.h`、`src/io/result_reader.cpp` | 结果 provenance 保存 8 顶点 | 版本化拓扑标识及连接数据，保留旧结果读取 |
| `src/io/postprocess.cpp`、`src/parallel/parallel_writer.cpp` | 输出及核对假设 hex 顶点连接 | 输出真实 polyhedral faces；更新并行结果闭环 |

特别注意：当前 `mesh.cpp` 的正交扩散系数为 `(Sf·d)/(d·d)`，不是其他历史版本中的 `|Sf|²/(Sf·d)`。迁移必须以当前源码和实测基线为准；如需修改非正交数学方案，应独立提交并进行受控数值对照。

## 3. 可以参考 Fine/Marine 的哪些内容

ISIS-CFD 开发团队的公开论文《Towards industrial use of anisotropic adaptive mesh refinement in CFD》明确描述：未知量 cell-centered；离散 face-based；任意数量/形状的组成面；分裂面和 hanging nodes 按多个独立面处理。它支持本方案的计算拓扑设计，但没有公开证明其内部 C++ 布局、索引宽度或专有优化实现。

- 开发者论文：<https://internationalmeshingroundtable.com/assets/papers/2023/03-Wackers-compressed.pdf>，第 2 节。
- LHEEA 官方 ISIS-CFD 介绍：<https://lheea.ec-nantes.fr/research-impact/software-and-patents/isis-cfd?l=1>。
- 可公开核对的 owner/neighbour 多面体格式：<https://www.openfoam.com/documentation/user-guide/4-mesh-generation-and-conversion/4.1-mesh-description>。

下述数组、API、算法分层和性能门禁是针对 BabelSim 的设计，不冒充 Fine/Marine 的源码实现。

## 4. 数据模型与层次

建议概念分层，保持现有 Mesh 公共入口的调用边界，但将其唯一实现收口到通用面拓扑：

```text
通用显式面输入（旧 hex 算例先离线迁移）
              ↓  构造、校验、定向
MeshTopology：points + faces + owner/neighbour + patches
              ↓  初始化缓存
MeshGeometry：Cf, Sf, |Sf|, Cc, Vc, interpolation/nonorthogonal data
              ↓  face-based 离散
DiscreteEquation：diag[cell], upper[face], lower[face], source[cell]
              ↓  已缓存 face→coefficient 映射
现有压缩稀疏矩阵 / CSR + 现有线性求解器
```

建议连续数组：

```cpp
// 示意，须根据现有 MeshStorage 和 C++17 风格实现。
points[nPoints];
facePointOffsets[nFaces + 1];
facePointIds[sumFacePoints];
owner[nFaces];
neighbour[nFaces];                 // physical boundary = invalid_index
facePatch[nFaces];
cellFaceOffsets[nCells + 1];
cellFaceIds[sumCellFaces];
cellFaceSigns[sumCellFaces];        // owner +1 / neighbour -1
// 由连接关系生成的去重 cell 邻接图，用于分区、BFS、AMG。
cellNeighbourOffsets[nCells + 1];
cellNeighbourIds[sumUniqueNeighbours];
```

主数据是 points、face 顶点环、owner/neighbour 和 patch。cell_faces、cell_neighbours 为一次性派生缓存，不能要求调用者同时提供三套连接真相。cell_vertices 如确需输出时使用可派生 CSR；不得重新变为主拓扑。

提供 C++17 只读 range：`facePoints(f)`、`cellFaces(c)`、`cellNeighbours(c)`。热路径不返回临时 vector，不经虚函数，不做散列查找。构造器可使用 vector-of-vectors/map，finalize 后展平为连续数组。

局部 ID 保留 int32 以避免无依据地增大 SpMV 带宽；构建阶段检查 ID、连接总长和 MPI count 溢出。offset 使用适合总连接长度的明确类型；如扩展 64 位全局 ID，单独设计和测试，不能静默截断。

静态网格 finalize 后冻结拓扑和几何。符号图、稀疏位置、分区和 halo 映射按 mesh 生命周期复用。若未来增加拓扑/几何版本，明确哪些缓存随哪个版本失效；当前无需引入动态网格。

## 5. topology 不变量与计算子面

- 每面至少 3 个有效、互异顶点，组成有序简单闭环。
- owner 始终有效；内部面 neighbour 有效且与 owner 不同；全局内部面恰好关联两个 cell。
- 内部面顶点定向使 Sf 朝 owner 外侧、朝 neighbour；物理边界朝 owner 外侧。cell-face 符号引用与之相符。
- 物理边界恰好属于一个物理 patch；处理器耦合不应落入物理 BC 计算。
- cell 必须闭合、定向一致、正体积；拒绝自交、非流形、悬空面、重复面和非法索引。
- 校验 `sum(sign*Sf)` 的归一化闭合误差；同时检查边连接、面积、体积和面矩，不以一个总和替代完整有效性检查。
- 顶点环 identity 应规范化循环起点和反向，不能把任意排序后的顶点集当作多边形环。

粗面连接多个细 cell 时，输入/适配器把界面拆成共同覆盖的多个计算面：每个子面一个粗 cell owner 和一个细 cell neighbour，粗 cell 的 cellFaces 自然增长。不能让一个 face 带 neighbour 列表。首轮提供显式子面输入和测试，不声称已经实现自动几何求交。

子面之间无重叠/缝隙，合成面积向量和面矩应匹配母面；母面可保留 provenance，但不得与子面一起参与通量积分。细分产生的边点必须在相邻面环中一致，闭合检查不能把一个长边与两个短边误判成闭合。

同一对 cell 可能通过多个子面相邻：cell-face 关联保留所有面；分区/AMG 的 cell-neighbour 图去重；矩阵系数对这些面求和。这三种表示不能混为一谈。

## 6. 通用几何和离散数学

多边形面采用确定性的有效三角化，计算每个有向三角形面积向量、面积矩及体积贡献。凸、平面面是首批基线；简单凹多边形必须使用适用的三角化，不能一律以顶点或平均点作无条件扇形划分。

非平面面须定义唯一三角化表面及 warpage 指标。第一版对超过阈值的输入明确拒绝，或在导入阶段显式拆成平面计算面，并记录转换；不能让 owner/neighbour 各自选择不同三角化。保留旧扭曲 hex 的兼容几何需要专门 fixture 和数值对照，不能不加说明地改变旧算例。

cell 体积和一阶矩用闭合面三角形的有向四面体贡献累计。取靠近 cell 的参考点减小平移消减误差：

```text
Vt = sign(cell,face) * dot(a-r, cross(b-r,c-r)) / 6
V  = sum(Vt)
M  = sum(Vt * (r+a+b+c)/4)
Cc = M/V
```

不能对每个 Vt 取 abs；那会掩盖反向和凹体问题。面法向定向依赖拓扑闭合及有向体积校验，不能单凭 vertex-average 中心判断一般凹 cell。验证平移/缩放不变性、解析体积/质心、面积闭合及散度定理。

在已验证的几何上，连续方程面通量仍写为：

```text
residual[owner]     += flux(face)
residual[neighbour] -= flux(face)  // 内部面
```

存储积分形式还是体积归一化形式须沿用当前方程约定。验证质量守恒时用对应的体积权重，不能混用量纲。

扩散使用明确分解 `Sf = alpha*d + k`。对方程 `-div(Gamma grad(phi))=q`，内部面隐式块是 `Gamma_f*alpha * [[1,-1],[-1,1]]`，显式非正交通量按 owner/neighbour 相反符号入 RHS。迁移先保留当前 alpha；数值增强另作实验，尤其不能假设任意 polyhedron 上两点通量就准确。

least-squares 梯度、面插值、skewness 修正、限制器、边界法向距离、Rhie–Chow、压力修正、壁面距离/力积分均须审计。通用拓扑不自动保证二阶精度或非正交稳定。病态 LS 系统需诊断/拒绝或明确回退，不得制造有限但错误的结果。离散适用性检查如 `Sf·d>0` 与拓扑有效性检查分开报告。

Solver 仍按：读取场和参数 → 物性 → 清空系统 → 逐项组装 → BC/约束 → 求解 → 修正场 → 收敛检查。保留当前 equ/math API，避免网格重构夹带 DSL 重写。

## 7. 不退化的矩阵装配与线性求解

保留 diagonal/upper/lower 和现有 Eigen compressed sparse / interior-boundary CSR。Krylov 内核只接触代数数组和预建 halo，不读 face 顶点、cell-face CSR 或几何。

SparseAssembly 初始化时建立：

- `diagonalPosition[cell]`。
- 单面耦合的 `upperPosition[face]` / `lowerPosition[face]`，继续直接写。
- 重复耦合位置的连续 group offsets / face IDs，每组确定性求和后写一次。

更新时单面路径不增加每 face 的动态容器或分支；重复组每次从零起和，避免跨方程累积。清零/覆盖范围明确，包含系数从非零变成零的情况。线程化若将来引入，按矩阵位置分组避免竞争，不在本次强制增加原子操作。

MPI remote couplings 和 AMG coarse assembly 同样审计多面求和。halo 以唯一 ghost cell 为单位发送，不能按子面重复发送同一未知量。保留当前收发重叠策略、残差定义、预条件器更新频率。

AMG 的 `6*nCells` 粗化输入改为唯一邻接 CSR：可先发布 owned row 的 degree，再构建全局 offsets 并交换 packed neighbour IDs，或用等价的 packed rows 协议。明确这是保留现有全局图策略的兼容改造，不宣称具有新的分布式扩展能力。不得简单换成 `maxDegree*nCells`。在每个 mesh/partition 初始化时完成，不放到每个方程/迭代。

旧 hex 的邻接顺序、分区结果和 AMG 聚合可能影响迭代数。兼容测试保留或显式映射原编号；若为通用性改变顺序，分别比较冻结分区/冻结矩阵的代数性能与完整流程的实际收敛性能，不用前者掩盖后者退化。

## 8. MPI 与文件格式

MPI 直接从全局 face topology 抽取局部点、面、owned/ghost cell，保留 globalCellId/globalFaceId 及相对全局方向。ghost 深度继续由算子需求决定。

拥有 owner 和 neighbour/ghost 的分区接口仍是内部耦合面。局部外层 ghost 的截断边界需标记处理器/halo closure 并保留远端映射，不可当作真实物理边界；对 owned cell 生效的 stencil 必须齐全。允许局部为有效 owner 重定向，但需同步面环、符号和通量映射。区分 global face、local face copy 和负责输出/统计的 rank，禁止重复全局积分。

v3 主数据：header/version、points、faces(vertex count + ordered IDs + owner + neighbour + patch)、patches(name/kind)。cell count 和 cell_faces 从 owner/neighbour 导出。校验内部面/边界面覆盖，允许合法空 patch，限制和诊断超大 count。v2 只由 `tools/convert_hex_mesh_v2_to_v3.py` 一次性转换，运行时 reader 拒绝 v2；新格式 round-trip 必须保持拓扑及 patch 语义。

输出支持 VTK polyhedron 的真实 faces/faceoffsets 等合法连接，具体 writer 采用当前 VTK 官方格式验证；不能只写 cell 顶点列表并标成 hex。结果 mesh provenance 版本化并包含足以检测连接/patch 变化的信息；旧结果按旧契约读取。串行和 MPI 输出均检查 cell IDs、字段长度、几何和后处理读取。

## 9. 性能门禁与证据

本次仅制定门禁，没有测量结果。P0 先建立可复跑的旧版 release 基线；后续每阶段复测受影响部分。

冻结条件：同一机器、编译器、优化参数、Eigen/MPI 版本、rank 数、绑定、线程数、输入、行列/分区映射、初值、容差、预条件器、停止准则、输出开关。明确记录 CSR_SPMV/ASYNC_HALO 等实际 build 开关。

两条基准必需同时存在：

1. 冻结 A/b/x0 的代数 replay，保留必要的远程耦合、分区及 AMG 图元数据；避免重新离散的数值差异干扰判断。SPD 扩散、非对称对流扩散、压力矩阵使用当前合法的 CG/IC、BiCGSTAB/ILUT、AMG 组合。
2. 旧 hex 完整路径，从网格读入、几何/图初始化、组装到收敛，单独计算初始化、重复求解和总成本。

记录：cells、faces、nnz、每行 nnz 分布、重复耦合数量、ghost/halo 字节、峰值内存；geometry/setup、symbolic pattern、numeric update、preconditioner setup、SpMV、preconditioner apply、halo、reductions、solve 总时间；迭代数和真实 residual。

至少一次预热和 7 组交错 A/B 重复，数据不足时补样；报告中位数、离散度和比值置信区间。建议以 2% 为性能回归告警线，而不是允许损失的预算：若稳定复现超线退化必须修复；若噪声覆盖边界则标记证据不足，不宣称“零损失”。任何稳定可归因于重构的退化都要解释和消除，除非用户另行接受。

固定迭代数只衡量吞吐，必须同时有达到相同真实残差的 time-to-solution。MPI 测量用 rank 最大时间而不是只看 rank 0。不同网格同时报告 nnz 和每 nnz 时间，不能只按 cell 数比较。测试规模包含缓存内小矩阵和超出 LLC 的大矩阵；资源允许时采用现有 100k/1m case 配置，不能无条件启动超资源作业。

## 10. 阶段与验收

| 阶段 | 交付 | 出口门禁 |
|---|---|---|
| P0 审计/基线 | 固定假设清单、release 基线、代数 replay/报告 | 基线可复跑，矩阵/残差有证据 |
| P1 通用拓扑 | CSR/range、builder、显式 subfaces、离线旧算例迁移 | 多边形面、>6 面 cell、坏网格；旧入口不存在；完整构建可用 |
| P2 通用几何 | 有向几何、定向/闭合/质量检查 | 解析体积/质心、平移缩放、散度定理、旧 hex 差异说明 |
| P3 组装/代数桥接 | 多面汇总、AMG CSR 图、缓存 | 独立 dense 矩阵 oracle、A*x、重复 update、无迭代路径退化 |
| P4 离散/物理接入 | grad/div/laplacian、BC、SIMPLE/RANS 审计迁移 | 线性再现、MMS 网格序列、守恒、压力修正一致性 |
| P5 MPI | 通用图分区、局部面抽取、halo | 1/2/4 ranks 正确性；8 ranks 资源允许；分区子面守恒与性能 |
| P6 IO/输出 | v3、离线 v2→v3 工具、结果版本、polyhedron 输出 | round-trip、迁移后的旧算例/结果、串行并行后处理闭环 |
| P7 综合验收 | 数值/性能报告、文档与旧依赖清理 | 旧 hex 性能门禁、poly 数值门禁、准确列明未验证项 |

按依赖推进，每阶段允许为完整构建作最小调用点迁移，但不得借此提前宣布后续阶段已验证。P1/P2 在同一 feature branch 内渐进实现；已有用例必须可运行，真正 polyhedral 求解直到 P3/P4 出口才算接通。分支建议 `codex/polyhedral-fvm`；不要默认 commit/push。

关键测试样本：单 tet/prism/pyramid/cube、五边形柱体、切角立方体、明确凸多面体、多个 face 连接同一对 cell、1:2/1:4 粗细接口；非法非流形/自交/零面积/负体积/缺失或重叠子面。凹面及扭曲面根据声明的支持范围分别验证接受或有诊断地拒绝。

数值测试：常数场保持、解析线性梯度、全局/局部守恒、解析扩散/对流扩散 MMS（至少三层、预先声明范数和阶数目标）、压力 nullspace/reference、Poiseuille、cavity。SIMPLE/RANS 跑通或短 pilot 仅说明集成，不等于湍流物理验证。不能修改物性、松弛因子或容差来隐藏拓扑错误。

当前已有入口：`make test-architecture`、`make test`、`make test-workflow`、`make test-mpi`、`make test-simple-parallel`、`make test-rans`，以及 `tools/benchmark_backend.py` / `tools/summarize_performance.py`。先检查参数和依赖，再运行对应阶段所需检查；未完成的长任务不能报告通过。

完成定义：通用面拓扑 + 通用几何 + 正确代数汇总 + 有证据的离散精度/守恒 + MPI + IO + 旧版线性性能门禁全部成立。最终报告分别列“源码实现、已完成测试、性能证据、物理验证、剩余限制”，提供 HEAD、工作区和实际修改文件。
