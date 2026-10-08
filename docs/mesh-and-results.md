# 网格、场与结果文件

本页按当前 [mesh reader](../src/io/mesh_reader.cpp)、[field reader](../src/io/field_reader.cpp)、
[Case](../src/io/case.cpp) 和 [后处理实现](../src/io/postprocess.cpp) 说明文件契约。

## 原生网格

运行时只接受 `BABELSIM_MESH 3`。三维几何显式存储顶点及有向面；单元由面所属的
owner/neighbour 构成，不按固定八顶点或六面假设读取。文件结构为：

```text
BABELSIM_MESH 3
vertices <N>
<x> <y> <z>
...
faces <F>
face <vertex-count> <owner> <neighbour> <patch> <vertex-id> ...
...
patches <P>
patch <name> <kind>
...
end
```

索引从 0 开始。内部面有有效 neighbour、无物理 patch；边界面 neighbour 为 `-1`，
patch 指向 patch 表。顶点环的面积矢量朝 owner 外侧，内部面朝 neighbour。
`kind` 可为 `generic`、`wall`、`inlet`、`outlet`、`symmetry`（兼容 `mirror`）、`processor`。
案例文件应描述物理边界；MPI 分区生成 processor patch。

reader 拒绝非有限坐标、错误版本、缺失段落和尾部额外 token；Mesh 构造进一步检查
面/单元连接、方向、闭合和正体积等几何约束。实际文件示例见
[heat.mesh](../cases/heat/mesh/heat.mesh)，构造 API 见 [mesh.h](../include/babelsim/mesh.h)。
网格文件不支持把 C++ 测试 fixture 的容错行为作为输入规范。

## Case 字典与路径

`case.bs` 必须给出 `solver`、`mesh`、`fields`、`physics`、`methods`、`solution`、
`control`、`output`，`ghostLayers` 可选、默认 3，要求不小于 3。推荐文件路径相对
Case 根目录，reader 按根目录拼接；它并未统一禁止输入文件使用绝对路径。
结果 `directory` 必须相对 Case，不能包含 `.`/`..` 路径段。

`physics` 管物性，`methods` 管离散，`solution` 管线性与算法数值控制，`control` 管
`startTime`/`endTime`/`deltaT`，`output` 管输出。重复或未知键以及未消费的已配置参数
会被拒绝。方法和线性键见 [数值配置](numerical-configuration.md)。

## 初场和重启

初场使用带括号的 `field <name> { ... }`，必需 `type`、`location cell`、`internal` 和
全部物理 patch 的边界。标量、向量、张量分别为 1、3、9 个值，实例见
[T.field](../cases/heat/fields/initial/T.field) 和 [U.field](../cases/poiseuille/fields/initial/U.field)。

`internal uniform <value>` 设置均匀值；`internal file <path>` 加载非均匀值。路径相对
`.field` 文件，可读取 `globalCellId` 加分量的普通文本、单个结果 CSV，或合并一个
时间目录的 `rank-*/<field>.csv`。每个全局单元必须恰好出现一次，值需有限。
MPI processor patch 由 halo 提供值，不需要为它写独立输入条目。

重启需给每个未知场选择同一物理时刻的结果，并把 `startTime` 设为该时刻。时间历史
不在初场文件中；BDF2 重启后的首步用 Euler，再积累下一层历史。`-time` 仅命名输出，
不自动加载旧结果。使用新的运行名，避免混入旧分区结果。

## 输出与时间目录

```text
directory results
timeName final
writeInterval 10
writeFields U p
```

`writeFields` 可选白名单，`excludeFields` 可选黑名单。加载的 cell 场默认选择输出；
程序创建的场需由 Solver 的 `problem.output(field)` 或白名单选择。

| 运行方式 | 序列目录 | 最近写出的快照 |
| --- | --- | --- |
| 不带 `-time` | `results/<物理时间>/rank-XXXX/` | `results/<timeName>/rank-XXXX/` |
| `-time demo` | `results/demo/<物理时间>/rank-XXXX/` | `results/demo/rank-XXXX/` |

`Case::write()` 同时写当前时刻和最近快照；`finish()` 结束 Case 并保留最终快照。
它们负责 I/O，不判物理收敛。内置求解器只写接受的步；失败后以前接受步的结果可能仍在。
SIMPLE 达到稳态判据后写出一次。物理时间目录使用数值格式，例如 `0.01`、`0.05`。

每 rank 只写 owned 单元，包含 `metadata.bs`、`mesh.geometry` 和字段 CSV；CSV 以
`global_id,x,y,z,value0` 开头，向量/张量追加分量。元数据、分区数、global ID 覆盖和
网格一致性用于读取校验，不能把 ghost 重复写入当作完整覆盖。

## 后处理与对照

```bash
build-petsc/babelsim-post -case cases/heat -time demo/all -format vtk
build-petsc/babelsim-post -case cases/heat -time demo/latest -format vtk
build-petsc/babelsim-post -case cases/poiseuille -time mpi4 -format vtk tecplot
```

`all` 处理未命名序列，`demo/all` 处理命名序列并生成 `post/demo/series.pvd`；`latest`
选择最后保存的物理时刻。指定运行名处理最近快照。VTK 使用真实 polyhedron 拓扑，
Tecplot FEBRICK 限于兼容八顶点六面体结果。

串并行比较选择同一时刻，命令示例：

```bash
python3 tools/compare_parallel_results.py \
  cases/poiseuille/results/serial cases/poiseuille/results/mpi4 \
  --atol 5e-6 --rtol 5e-6
```

比较通过表示给定输入及容差下的结果一致，不能代替解析解、物理基准或网格/时间步研究。
