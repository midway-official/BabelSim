# 网格图划分与负载均衡模型

这项实验先在 Python 模型中比较通用分区策略，再选取候选接入 Parallel 层。
模型不调用物理求解器。其通信和存储指标来自真实网格拓扑，划分时间不能替代
SIMPLE/PISO 的实际计算时间，也不能预测 AMG 的迭代数或各层通信开销。

模型选择结果见 [2026-10-07 实验报告](../reports/partition-model-2026-10-07.md)。
工程规模 Voronoi 基底另用完整 16,000 多边形的实际面环验证，避免仅靠小型几何检查。
加密采样为 `0.002 + 0.996*x^3`，保持与镜像边界的数值可分辨性。

## 网格范围与契约

每个顶点代表一个控制体，每条无向边代表一个内部共享面。CSR 无自环、无重复邻接，
每条边双向出现，物理边界面单独计数。ID 可以打乱，几何中心与邻接不依赖 ID 顺序。
所有权唯一，owned 单元数保持当前框架的精确配额：`N/P + (rank < N%P)`。

实验包括 13 类可实现的 polyhedral 网格：六面体立方域、薄板、长通道、刚体旋转的长通道、
带圆柱孔洞的域、分支通道、沿流向分级加密的通道、扭曲棋盘三角棱柱、每立方体六四面体的
一致细分、变面数 Voronoi 多边形挤出、局部密集的 Voronoi 挤出、两个不连通域，以及
打乱单元编号的立方域。

生成器也能输出小型 `BABELSIM_MESH 3` 文件。模型邻接与生成的实际面环独立对齐；
框架 `readMeshFile()`/`Mesh::validate()` 检查正体积、面闭合、面方向及连接。
这两类检查验证生成器及其拓扑映射。工程规模图由相同构造生成，并逐图检查 CSR 与物理面数。
三层 halo 还与独立的广度优先参考计算比较，避免把一层指标当作框架需要的完整 stencil。

## 候选及依据

- `graph-growth`：当前 C++ 实现的确定性散布种子和配额轮流生长策略的 Python 模型。
- `rcb`：沿当前区域最长坐标方向递归切分，精确分配各子树配额。
- `rib`：根据单元中心的主惯性方向递归切分，可处理任意分区数量。
- `morton`：逐坐标归一化后，以 Morton 次序分配精确配额。
- `morton-isotropic`：使用相同的全局尺度归一化，保留物理长宽比，再按 Morton 次序分配配额。
- `metis-cut`：多级图划分，以割边为目标，只约束单元数。
- `metis-cut-work`：同时约束单元数和面关联数，以割边为目标。
- `metis-volume-work`：使用相同工作量约束，以一层通信量为目标，并减少分区邻接度。
- `metis-volume-contig`：在连通图上进一步请求连续分区；不连通图采用前一策略。

[METIS 手册](https://karypis.github.io/glaros/files/sw/metis/manual.pdf)描述了割边、通信量、
多约束负载和连续分区选项。模型直接调用匹配公开头文件的本机 METIS 库，固定每个候选的
`NCUTS=2`，分别考察种子 0、7、29。单元约束容差为 1.001，面关联工作量约束为 1.05。
其结果再恢复精确单元配额，单独记录修正时间及迁移数。

[Zoltan RCB 文档](https://sandialabs.github.io/Zoltan/ug_html/ug_alg_rcb.html)和
[RIB 文档](https://sandialabs.github.io/Zoltan/ug_html/ug_alg_rib.html)提供几何切分的参考定义。
这里的 RCB/RIB 是 Python 对照实现，未调用 Zoltan。多级图方法的研究背景见
[Karypis 图划分项目](https://papers.karypis.org/glaros/projects/gp.html)。

恢复精确配额时，优先在 surplus/deficit 分区之间迁移使割边增加最少的边界单元。
当二者不邻接时，通过分区图传递配额；优先用至多 64 个已展开节点的局部路径证明迁移后原分区的 owned 邻居仍相互连通。
用割边代价优先队列选择候选，迁移后仅更新该单元及邻居的代价。
同一所有权状态下缓存每个单元的连通判断；所有权变化后失效。找不到局部证据时优先尝试其他候选；只有必要迁移才完整检查。记录确实切断连接的强制迁移。
这个预算限制候选筛选的工作量，不是 halo 宽度，也没有截断求解器的 stencil。
不连通域和局部配额约束可能需要不连续的所有权，因此显式记录每分区连通块数与强制迁移，
不承诺任意图都能同时满足精确配额及连续性。

## 指标与选择

主要比较三层 ghost 总数及最大 rank 的 ghost 数，最大 rank 的局部面数及负载偏差，
邻接 rank 数、面关联工作量偏差、划分时间，以及单元和面的重复存储倍率。
一层通信量和割边单独报告。模型邻居数指 ghost cell 的 owner rank；面场及 AMG 的邻居另行测量。候选的几何平均比值说明整体表现，最差比值检查退化；
完整逐图数据保留，以免平均值遮盖某一网格族的问题。

`ownedFaceWorkImbalance = max(sum_owned(cellFaceCount)) / mean(sum_owned(cellFaceCount))`。
面关联数是通用装配和几何循环工作量的估计，不是浮点操作数，也不表示各 AMG 层真实负载。
`cellStorageDuplication = sum(localCells)/N`；面存储倍率由每个局部子网格所附着的真实面计数。
模型的进程 RSS 含 Python、共享图和内存分配器，实际求解器内存还包含场、方程、PETSc 和 AMG。

所有任务逐项运行，并使用与求解器基准相同的排他锁。
图生成、读取、验证、Numba 编译、指标计算及绘图均在划分计时之外。
图生长及配额修正使用 Numba 编译的 Python 内核，几何对照使用 NumPy，METIS 使用其原生库。
这减少 Python 解释器开销的干扰，但仍应把模型时间理解为实现候选的筛选证据。
最终求解器性能必须用相同案例、精度、固定工作量及完整计算窗口再验证。

## 重现

依赖为 Python 3.10+、NumPy、SciPy、Numba、psutil、Matplotlib，以及 ABI 匹配的 METIS
公开头和共享库。本次使用 NumPy 2.2.6、SciPy 1.15.3、Numba 0.62.1、METIS 5.1.0。
本机系统 SciPy 1.8 与 NumPy 2.2 不匹配，实验依赖放在独立目录，未改动现有求解环境。

```bash
python3 tools/study_mesh_partition.py \
  --graph-root /tmp/babelsim-partition-model/graphs \
  --output /tmp/babelsim-partition-model/broad \
  --scales small engineering --ranks 4 12 --seeds 0 \
  --generate --run

python3 tools/study_mesh_partition.py \
  --graph-root /tmp/babelsim-partition-model/graphs \
  --output /tmp/babelsim-partition-model/repeated \
  --scales small engineering --ranks 4 12 --seeds 0 7 29 \
  --algorithms graph-growth metis-cut-work metis-volume-work metis-volume-contig \
  --run
```

原始每次记录、图参数与哈希、源码快照、运行版本、CSV 汇总和图表保存在 output。
`--resume` 要求代码和研究配置相同；算法或修正逻辑变化使用新证据目录。
小型原生几何文件可使用 `--scales fixture --generate --fixtures` 生成。
