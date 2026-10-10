# 03 首次运行：热传导

[用户手册](README.md) · 上一章：[安装与构建](installation.md) · 下一章：[算例管理](case-management.md)

本教程完成一个 128 单元的热传导计算：左侧温度为 1，右侧为 0，其余边界不传热；初始温度为 0，计算时间从 0 到 0.05，时间步 0.01。参数采用一致单位，温度数值可视为归一化量。

## 1. 创建独立算例

先完成 `make`，然后执行：

```bash
python3 docs/examples/create_case.py --template heat --output work/heat
```

脚本生成网格、完整字典和初始场，不依赖原有结果。若目录已存在，脚本拒绝覆盖；需要重复教程时使用新的输出目录。生成的 `work/heat/case.bs` 为：

```text
solver heat
mesh mesh/domain.mesh
fields fields/initial
physics physics/model.bs
methods numerics/methods.bs
solution numerics/solution.bs
control control.bs
output output.bs
```

主文件只负责选择求解器及各输入路径。`physics/model.bs` 为：

```text
density 1
heatCapacity 1
conductivity 0.1
source 0
```

`control.bs` 为：

```text
startTime 0
endTime 0.05
deltaT 0.01
```

这些物性对应 `density * heatCapacity * dT/dt - div(conductivity * grad(T)) = source`。导热率增加会使相同时间内热量传播得更快；`source` 是单位体积热源。

## 2. 调用求解器

```bash
build/babelsim-solve -case work/heat -time serial
```

`-case` 指向包含 `case.bs` 的目录；`solver heat` 选择实际求解算法；`-time serial` 设置本次结果名称。它不会改变物理起始时间，也不会加载名为 serial 的重启场。

成功时输出 5 条以 `heat time=` 开头的诊断记录，进程退出码为 0。在终端可以紧接着执行 `echo $?` 检查退出码。不同机器的残差末位可能不同，不应将某个完整日志字符串作为数值标准。

## 3. 查看保存的文件

```text
work/heat/results/serial/
  rank-0000/
    metadata.bs
    mesh.geometry
    T.csv                  最近一次显式写出的场
  0.01/                    第一个保存时刻
  0.02/
  0.03/
  0.04/
  0.05/                    本教程最后时刻
```

每个时刻的 rank 目录均有结果元数据、几何文件和 `T.csv`。CSV 头为 `global_id,x,y,z,value0`，包括单元全局编号、中心坐标和温度值。这里保存单元温度；固定边界值不是作为额外单元写出，因此最热单元一般低于边界温度 1。

## 4. 导出时间序列

```bash
build/babelsim-post -case work/heat -time serial/all -format vtk
```

在 ParaView 中打开 `work/heat/post/serial/series.pvd`，选择温度字段 `T` 并播放时间序列。导出程序使用实际网格连接，而非仅将 CSV 转成点云。没有安装 ParaView 也可以先检查 `.vtu` 和 `.pvd` 已生成，或直接用 CSV 做定量分析。

只导出最后时刻和 Tecplot 文件：

```bash
build/babelsim-post -case work/heat -time serial/latest -format vtk tecplot
```

对应输出 `work/heat/post/serial/0.05.vtu` 和 `0.05.dat`。本教程采用八顶点六面体，适合当前 Tecplot FEBRICK 导出。

## 5. 使用两个 MPI 进程

```bash
mpirun -np 2 build/babelsim-solve -case work/heat -time mpi2
python3 tools/compare_parallel_results.py \
  work/heat/results/serial/0.05 work/heat/results/mpi2/0.05 \
  --atol 1e-8 --rtol 1e-8
```

MPI 结果分为 `rank-0000/`、`rank-0001/`，每个全局单元只由拥有它的进程写一次。比较工具按全局编号比较场值，允许浮点归约与 Krylov 路径引入的小误差。比较必须使用相同输入与相同物理时刻。

## 6. 修改与重复运行

将 `physics/model.bs` 中的 `conductivity` 改为 `0.2`，然后用新运行名保存：

```bash
build/babelsim-solve -case work/heat -time higherConductivity
```

保留参数修改记录，再比较 `serial/0.05` 和 `higherConductivity/0.05`。运行名本身不保存旧字典；正式工作应为不同设置复制独立输入目录或使用版本管理。

下一步可按[详细操作示例](tutorials.md)完成输运和流动计算，或者在[外部求解器教程](../developer-guide/extending-solvers.md)中编写并调用自己的热传导求解器。
