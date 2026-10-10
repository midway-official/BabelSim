# 文件格式参考

[配置参考](README.md) · [网格与边界手册](../user-guide/mesh-and-boundaries.md) · [结果手册](../user-guide/results-and-restart.md)

## BABELSIM_MESH 3

语法按空白分词，结构顺序固定：

```text
BABELSIM_MESH 3
vertices <N>
<x> <y> <z>
... 共 N 行
faces <F>
face <Nv> <owner> <neighbour> <patch> <v0> ... <vNv-1>
... 共 F 行
patches <P>
patch <name> <kind>
... 共 P 行
end
```

上面的尖括号为说明，不是文件中的字面文本。可直接运行的完整文件见[one-cell.mesh](../examples/one-cell.mesh)。

顶点、单元、patch 编号从 0 开始。内部面 `neighbour >= 0` 且 `patch=-1`；边界面 `neighbour=-1` 且 patch 指向有效边界。面环按 owner 外法向排列，内部面朝向 neighbour。单元由面连接定义，不另写 cell 行。

patch kind 可用 `generic`、`wall`、`inlet`、`outlet`、`symmetry`、`processor`，`mirror` 可作为 symmetry 别名。用户物理网格通常不手工构造 processor patch，它由运行分区产生。名称必须适合在场字典中精确引用。

文件要求有限顶点坐标、合法面连接、有效编号和明确 `end`，结尾不能包含额外未知数据。几何构造还检查封闭性、方向与正体积。格式 1/2 或其他工具的网格不是当前运行时输入格式。

## 场文件

```text
field C {
    type scalar
    location cell
    internal uniform (0)
    boundary {
        inlet  { type fixedValue value (1) }
        outlet { type zeroGradient }
        walls  { type symmetry }
    }
}
```

`type` 支持 scalar/vector/tensor，分量数为 1/3/9；tensor 行优先为 xx、xy、xz、yx、yy、yz、zx、zy、zz。加载位置为 cell。字段名、类型、位置、internal 和全部物理 patch 必须与求解器和网格一致。场文件路径为 case 的 fields 目录加 `<name>.field`。

boundary 接受 `fixedValue`、`fixedGradient`、`zeroGradient`、`symmetry`、`inletOutlet`。需要 value 的条件必须提供正确分量数；zeroGradient/symmetry 不提供 value。兼容别名为 fixed_value/dirichlet、fixed_gradient/neumann、zero_gradient、mirror。原生参考统一使用主名称。

## internal file

```text
internal file <relative-or-absolute-path>
```

相对路径以当前场文件目录为基准，可指向：

1. 空白分隔文本，每行 `global_id component0 ...`。
2. CSV，包含 `global_id,x,y,z,value0...` 表头。
3. 一个结果快照目录，由其中 `rank-*/<field-name>.csv` 合并读取。

例如标量文本：

```text
0 0.25
1 0.75
```

CSV：

```csv
global_id,x,y,z,value0
0,0.25,0.5,0.5,0.25
1,0.75,0.5,0.5,0.75
```

这些两行示例只能用于确有两个对应全局单元的网格。所有全局单元必须恰好出现一次，分量数必须正确。用户负责确保全局编号对应相同网格；坐标列不意味着提供自动网格插值或任意网格迁移。

## 单元结果 CSV

每个 rank 的标量头为 `global_id,x,y,z,value0`；向量追加 value1/value2；张量追加至 value8。坐标为单元中心，全局编号用于跨进程合并。只有拥有者写出，幽灵单元不重复写。

结果文件保存浮点场数据，不附带单位换算。定量分析应依据算例定义解释场含义。中间压力校正、PETSc 状态与 BDF2 全历史不在这些文件中。

## metadata.bs 与 mesh.geometry

每个 `rank-0000` 等目录保存自己的 metadata 与 geometry。metadata 记录结果格式版本、快照名称、rank、进程总数、全局/拥有单元数和输出场描述；不记录完整输入字典、求解器版本或全部算法状态。mesh.geometry 按拥有单元保存全局编号及顶点坐标，用于网格来源核对；完整网格连接由匹配的 case 网格提供。它们是与 CSV 配套的结果文件，应一起复制；不要手工修改元数据以掩盖不一致输入。

后处理会核对本次 case 网格、已保存几何、rank 集合、字段列与全局编号。原始输入网格修改后，旧结果不能继续以新几何导出。需要保留结果时，应同时保留匹配的 case 输入。

## 导出格式

VTK `.vtu` 保存单元几何与字段；`.pvd` 记录物理时间及对应 VTU 文件。一般多面体使用 VTK polyhedron 表达。Tecplot `.dat` 当前使用 FEBRICK，限制为适合的八顶点六面体连接。导出不进行重采样或恢复未保存的时刻。
