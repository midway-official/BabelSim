# 网格分区与负载分析

[开发指南](README.md) · [性能测量](performance.md) · [安装说明](../user-guide/installation.md)

BabelSim 从全局网格图生成 MPI 所有权，再构造 owned 单元与完整所需幽灵层。图分区属于通用并行层，求解器作者不在 Physics 内直接管理 ownership 或 MPI halo。

## 当前构建选择

```bash
make USE_METIS=0
```

新构建禁用可选 METIS，使用确定性图生长回退。默认 USE_METIS=auto 找到库则启用；USE_METIS=1 要求库存在。改变该选项应使用新 BUILD 或清理旧构建，避免旧对象沿用前次选项。

METIS 安装可以由系统提供，不属于必需的 PETSc/Eigen/Hypre 源码集合。框架若启用了 METIS，外部静态库使用者链接时也需带该库；减少运行依赖的交付可以使用回退版本。

## 图与数据契约

图顶点为单元，无向边为内部共享面。邻接不依赖全局编号顺序，不能假设相邻编号等于几何相邻。owned 单元唯一，分区保留全局 ID、物理 patch 和面方向；跨分区内部面仍是几何内部面。

当前分配目标的精确单元配额为 `N/P + (rank < N%P)`。METIS 路径同时考虑单元数和面关联工作量，并在需要时恢复精确配额。连通图可请求连续分区；不连通域与严格配额不能保证任意图都同时满足所有连续性条件。

ghostLayers 默认为 3，当前最小支持值也是 3。图模型的一层割边或通信指标不能替代求解所需完整 stencil。工作 rank 构造局部网格，root 解析与分发时的临时全局数据不代表持续保留同样的内存布局。

## 应检查的指标

| 指标 | 作用与限制 |
| --- | --- |
| owned 单元数 | 检查实际所有权配额 |
| 三层 ghost 总数与最差 rank | 评估完整 halo 的重复单元与最差进程 |
| 局部面数与面关联数 | 几何/装配工作量的估计，不是全部浮点操作数 |
| 分区邻接数、割边 | 图通信结构，不等于 KSP/AMG 所有通信 |
| 连通块与配额修正迁移 | 描述连续性及精确配额代价 |
| 初始化时间和峰值内存 | 包括图、分发、局部构造及运行库 |
| 实际求解时间、迭代数 | 用生产案例验证图选择的最终效果 |

更少 ghost 或割边不保证 AMG 更快，Python 图模型 RSS 也不等于真实求解器内存。性能报告应采用相同物理输入、容差和工作量。

## 可选图研究工具

tools/study_mesh_partition.py 用实际网格拓扑比较图生长、几何排序与 METIS 候选，独立于生产求解器。它需要额外 Python 科学计算环境与 ABI 匹配的 METIS；这些是研究工具依赖，不是 make 构建 BabelSim 的必需环境。

```bash
python3 tools/study_mesh_partition.py --help
```

准备该工具依赖后，可用小型原生几何检查生成流程：

```bash
python3 tools/study_mesh_partition.py \
  --graph-root work/partitionGraphs --output work/partitionEvidence \
  --scales fixture --generate --fixtures
```

该工具的实际选项与依赖以当前源码为准。图生成/合法性、模型划分时间、生产求解器计算时间和物理验证分开记录。旧研究的环境、数字和冻结源码已归档到 archive/docs，不作为交付软件的普遍性能保证。

## 修改并行层后的验证

运行架构与局部分区/halo 契约，再进行生产案例同一时刻 1/2/4 rank 字段比较。非正交与可变面数网格、编号打乱、不连通区域和跨分区面通量需要专项覆盖。详见[验证指南](verification.md)。
