# BabelSim polyhedral FVM 可执行 prompts

配套方案：`docs/polyhedral-fvm-refactor-plan.md`。以下是交给实现代理的指令，并不表示这些阶段已经执行。

## 总控 prompt（一次启动整个重构）

```text
在 /home/midway/BabelSim 执行通用 face-based polyhedral FVM 重构。
先阅读适用的 AGENTS.md、docs/polyhedral-fvm-refactor-plan.md 和本文件。
按 P0→P7 顺序实际实现、验证和修复，不停留在设计或仅增加接口。

必须满足：
1. face 顶点数及 cell 面数可变；每个计算 face 唯一 owner、至多一个 neighbour。
2. 同一通用拓扑支持迁移后的旧 hex 算例和任意 polyhedral；运行时不得保留 `Mesh::unstructured` 或其他旧 hex 网格入口。旧文件只由一次性离线工具转换为 v3。
3. 以连续 CSR/SoA 和 C++17 只读 range 实现；热路径无 vector-of-vectors、哈希查找、动态分配或按单元类型分派。
4. 绝不牺牲现有线性系统求解性能。保留 LDU、Eigen compressed/CSR、CG/BiCGSTAB、IC/ILUT/AMG 及 MPI 通信重叠。预处理拓扑与 coefficient mapping，迭代只使用代数数组。
5. 对同一对 cell 的多个计算 face，正确汇总矩阵贡献；普通单面耦合保留直接赋值快速路径。AMG 的 6*nCells 邻接改为 CSR，不能改成 maxDegree padding。
6. MPI 从原始通用 face 拓扑分区，保留 global IDs、方向、ghost 深度和守恒；不能由 hex 顶点重建。
7. 保持当前 equ/math 和 Solver 使用流程；非正交数值方案如需调整，独立证明，不能夹带到拓扑迁移。
8. reader 只接受 v3 显式面输入；提供独立 v2 hex→v3 迁移工具并迁移仓库旧算例，完成结果与 polyhedron 输出闭环。
9. 先建立旧版 release 基线，再建立相同 A/b/x0、分区、配置下的代数 replay；同时验证旧 hex 全流程。不得用放松容差、减少迭代或更换预条件器换取性能。
10. 不实现 AMR/ALE/VOF/overset/GUI 等额外模块，不默认提交或推送，不覆盖用户改动，不启动子代理。

开始时记录 HEAD、工作区、编译/MPI 环境。方案审计基线是 c6c9dfa，若当前 HEAD 不同，以当前源码重核事实，不回滚。
每阶段完成必要测试，记录命令、退出码、误差、性能及未完成项；失败则先定位修复再过门禁。
若已有基线测试失败，保留旧版失败证据并与新回归分开；不能削弱测试。
维护 docs/reports/polyhedral-fvm/ 下的阶段记录和可复跑报告；实际性能数据保留原始日志。
遇到资源限制或真实外部阻塞，记录可复现原因，完成所有独立可进行工作；不把未运行当通过。
最终分别报告实现、数值验证、MPI、性能、物理验证边界，以及 HEAD/工作区/改动文件。
```

## P0：基线与性能合同

```text
执行 polyhedral FVM P0。读取总方案和总控约束。
检查当前 HEAD/工作区、AGENTS.md、Makefile、实际工具链和性能脚本参数。
按 mesh→geometry→operators→assembly→Krylov/AMG→MPI→IO 检索固定 4/6/8 假设，区分维数常量与拓扑常量，记录文件和符号。
从现有版本建立可复跑 release 基线。需要隔离时使用独立 build 目录或 worktree，保留用户改动，不把工作区改动丢失到旧 HEAD 基线中。
记录已有测试结果；优先小规模关键测试，再执行基线必须的现有回归。
复用 benchmark_backend.py/summarize_performance.py；增加缺失的冻结 A/b/x0 代数 replay 及必要分区/remote coupling/AMG graph 输入，证明 replay 与原求解一致。
覆盖 SPD、非对称和压力问题，当前合法预条件器组合、串行及 MPI，缓存内与带宽受限规模。
记录 mesh/nnz/配置/hash、矩阵更新、预条件器 setup/apply、SpMV、halo/reduction、solve、真实残差和迭代数。
制定预热、交错 A/B 重复、统计区间和资源上限。2% 只作回归告警，不声称可接受 2% 损失。
产物：源码审计表、基线原始数据、复跑命令和报告模板。P0 不修改物理离散/求解算法。
```

## P1：通用 face topology

```text
执行 P1，沿用 P0 基线及总控约束。
在现有 Mesh/MeshAccess 边界内实现 facePointOffsets/Ids、owner/neighbour、cellFaceOffsets/Ids/signs、去重 cellNeighbour CSR、patch 映射和 C++17 只读 range。
由 faces/owner/neighbour 唯一派生 cell 关联，显式 cellCount，不再依赖 cell_vertices[8]。
增加显式面 builder；删除原 `Mesh::unstructured(hex...)` 及固定 hex 存储。所有旧算例先用离线转换工具写成 v3，再由唯一 reader 读入通用存储。
构建器可用临时动态容器，finalize 后只保留连续数据；检查 count/offset/ID/MPI 类型溢出。
实现连接、环方向、patch 覆盖、重复/非流形/闭合基础校验；保留后续几何检查入口。
最小迁移调用点以维持整个项目构建和旧 hex 测试。不得留下 array<6> 形式的运行时兼容影子拓扑。
测试 tet/prism/pyramid/hex、五边形面、>6 面 cell、同对 cell 多面、1:2 与 1:4 显式子面及非法输入。子面分裂时同步相邻面边点。
以计数、双向关联和方向等独立预期验收，不能只测试新 API 调新 API 的自洽性。
报告当前仅拓扑接通的能力，不宣称完成通用求解。
```

## P2：通用几何

```text
执行 P2，读取 P1 接口及总方案。
实现简单多边形的确定性三角化、Sf/Cf、cell 有向体积及一阶矩。不能把所有面当 quad，不能对每个四面体体积取 abs。
为平面凹面、非平面/warped 面、凹 cell 明确支持域和诊断；不支持的输入明确拒绝。超阈值 warped face 若选择切分，需要显式记录转换及一致 owner/neighbour。
分开拓扑有效性和离散适用性；校验闭合、定向、正体积、面矩、nonorthogonality、skewness 及 Sf·d。
预处理几何和权重；静态迭代不重算多边形几何。
保持当前扩散 alpha 及旧 hex 几何兼容，不擅自以历史版本替换当前公式。若有必要差异，用独立实验说明。
验证解析体积/质心、三角化一致性、平移缩放、散度定理、子面覆盖/面积矩和旧 hex 基线；容差按尺度归一化。
记录几何 setup 时间/内存及可支持的网格质量范围。
```

## P3：组装和代数桥接

```text
执行 P3，优先处理 SparseAssembly::update 的多面覆盖风险和 buildCoarseMapping 的 6 邻居假设。
保留 DiscreteEquation 的 diag/upper/lower/source 及现有求解器。
初始化时建立普通 face 直接 coefficient-position 映射和重复 cell-pair 的连续汇总分组。每次 update 覆盖/清零正确，禁止跨次累积；不对普通网格强制使用通用哈希装配。
检查对称/非对称 upper/lower、零系数、owner/neighbour 翻转、MPI remote coupling 及 AMG coarse assembly 的汇总。
将 AMG 图输入替换为变长去重 CSR，packed owned rows 交换及 global offsets 校验；缓存到拓扑/分区生命周期。保持现有聚合策略和预条件器流程，评估邻接顺序改变的影响。
用手工小矩阵或独立 dense oracle 核对多面 A/b、A*x、重复 update，覆盖一对 cell 两个及多个子面。
运行冻结 A/b/x0 代数 replay，检查迭代数、真实 residual、SpMV、预条件器 apply、halo/reductions 和总求解；检查 steady hot path 无新增分配或几何/拓扑访问。
任何可复现性能退化定位修复；证据不足就标记未验收，不能放松门禁。
```

## P4：离散和物理流程

```text
执行 P4，保持总方案的积分/体积归一化和通量符号约定。
迁移 operators.cpp、field_math.cpp、fvm_execution.cpp 等 stencil/邻接访问；统一逐 face 装配。审计梯度、插值、div/laplacian、LS 条件数、skewness/非正交修正、限制器、BC。
检查 SIMPLE 动量、压力修正、Rhie–Chow 质量通量修正、连续性诊断、压力参考/nullspace；检查 RANS 壁距及边界/力积分涉及的几何假设。
保留 equ/math DSL、物性/组装/BC/求解/修正顺序，避免给不同 cell 类型写独立数值路径。
测试常数保持、解析线性场、内部面抵消/整体守恒、至少三层网格 MMS、粗细子面界面。预先声明误差范数和预期阶数，失败分析几何、重构、边界和线性容差。
跑旧 hex heat/transport/Poiseuille/cavity 及通用网格 laminar SIMPLE，用相同输入比较；RANS 仅报告完成的集成或物理验证层级。
若需改非正交数学方案，单独提交级别记录和 A/B 对照，不能调松弛/容差来遮蔽错误。
复测 numeric assembly 和旧 hex 完整求解性能。
```

## P5：MPI 分区和 halo

```text
执行 P5，将 parallel_context.cpp 中广播、图分区和局部网格提取改为通用面拓扑。
发送 count+offsets+packed IDs；不再传每 cell 8 个顶点，不再调用 hex 构造器重建局部 mesh。
保留 globalCellId/globalFaceId、owner/neighbour、面方向/局部翻转映射、owned/ghost 和算子所需 ghost 深度。
处理器耦合与真实物理边界明确区分，外层 ghost closure 不影响 owned stencil。global flux/积分只计一次，局部贡献采用一致符号。
halo 以唯一 ghost cell 为单位，子面不导致重复未知量通信；保持异步通信与现有 CSR 热路径。
在 1/2/4 ranks，资源允许再到 8 ranks 验证 hex、mixed/poly、跨 rank 粗细界面和同对 cell 多面：场误差、真实 residual、质量守恒、global IDs/输出去重及重复运行。
分别比较冻结分区和实际分区的性能，记录 max-rank 时间、通信字节、ghost 比率、AMG setup 与收敛时间；不能以位级一致作为一般 MPI 唯一标准。
```

## P6：网格文件、结果与可视化

```text
执行 P6。增加 BABELSIM_MESH 3：points、faces(有序变长顶点/owner/neighbour/patch)、patches(name/kind)；cell count 与 cellFaces 派生。
reader 只接受 v3；用独立 `tools/convert_hex_mesh_v2_to_v3.py` 迁移旧算例后删除运行时 v2/hex 入口。严格校验索引/长度/版本/patch 和 count 溢出。
更新 result.h、result_reader.cpp、parallel_writer.cpp、postprocess.cpp 的 8 顶点 provenance 假设，版本化 mesh fingerprint/连接及结果元数据，兼容旧结果读取。
根据官方 VTK 文档输出合法 polyhedron 连接、faces/faceoffsets；用实际读取器验证 cell 数/拓扑/体积/字段和 MPI 合并，不只检查 XML 能解析。
测试 v3 round-trip、v2 旧例、错误输入、混合单元、>4 顶点 face、>6 面 cell、子面及串行/MPI 结果闭环。
固定数组只可存在于一次性迁移脚本内部；不把它们带入运行时 Mesh，不自动扩展到完整 Gmsh/CGNS 导入或修改所有历史结果文件。
```

## P7：综合验收和交付

```text
执行 P7。审计所有生产路径剩余的固定拓扑假设，区分合法的 3D 向量维数、预分配 hint、旧格式 adapter 和真正硬限制。
完成适用的 make test-architecture/test/test-workflow/test-mpi/test-simple-parallel/test-rans 及新增通用网格测试；按实际目标逐项记录退出状态，不能把启动视为完成。
执行 P0 的冻结代数与旧 hex 全流程 A/B 基准，预热后至少 7 组交错重复，报告中位数、离散度、比值区间、迭代数和真实残差，性能异常先定位修复。
通用网格报告 cells/faces/nnz、每 nnz 成本、守恒/MMS 阶数和收敛，不对不同网格无条件宣称总时间相同。
核对多面聚合、AMG 变长图、MPI 符号/halo 和输出读取的专项证据。
更新架构/用户文档、格式说明、已支持/拒绝的几何质量范围及后续扩展边界。
交付 docs/reports/polyhedral-fvm/ 中完整报告、复跑命令、原始证据索引和性能对比。
最终明确区分源码实现、自动测试、数值验证、性能验证及物理验证，列所有未完成项；输出 HEAD、git diff --stat、git status，不默认 commit/push。
只有全部出口满足才能宣布本次重构完成。
```
