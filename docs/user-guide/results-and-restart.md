# 09 结果处理与继续计算

[用户手册](README.md) · 上一章：[运行与并行计算](running-and-parallel.md) · 下一章：[详细操作示例](tutorials.md)

## 结果组成

每个写出快照包含 `rank-*/metadata.bs`、`rank-*/mesh.geometry` 和 `rank-*/<场名>.csv`，rank 名采用 `rank-0000` 的四位补零形式。CSV 只写每个进程拥有的单元，不重复幽灵单元。运行名目录保存最近一次显式写出的快照，其数值时间子目录保存序列。

上述命名序列使用显式求解 `-time baseline`。省略该参数时，数值时间目录直接位于 `results/0.01` 等位置，而最近快照写入 output 的 `timeName`（默认 `results/final`）；仅设置 timeName 不会把数值序列再嵌入同名目录。

`output.bs` 示例：

```text
directory results
timeName baseline
writeInterval 10
writeFields U p
```

`writeInterval` 按时间步编号计数。正式瞬态求解器在接受且满足间隔的步以及最终接受步写出；程序不会自动保存所有中间压力校正。默认输出已从文件加载的单元场；程序创建的派生场必须由求解器选择输出或列入 `writeFields`。完整字段说明见[字典参考](../reference/dictionaries.md)。

## 选择结果时刻

| 后处理 `-time` | 数据选择 |
| --- | --- |
| `baseline` | 命名运行目录中的最近写出快照 |
| `baseline/0.03` | 命名运行的指定物理时刻 |
| `baseline/latest` | 命名运行的最后一个已保存数值时刻 |
| `baseline/all` | 命名运行所有已保存数值时刻 |
| `latest`、`all` | 未命名序列根目录的最后/全部数值时刻 |
| 省略 | 使用 `output.bs` 的 `timeName` 选择 |

`latest` 代表最新已保存时刻，不代表求解器必然成功到达预定结束时间。应核对元数据与日志。

## 导出几何场

```bash
build/babelsim-post -case work/heat -time baseline/all -format vtk
build/babelsim-post -case work/heat -time baseline/0.03 -format vtk tecplot
```

时间序列写入 `post/baseline/`，包括每个时刻的 `.vtu` 与 `series.pvd`。单个指定时刻写 `0.03.vtu`/`0.03.dat`；选择 `baseline` 本身写 `post/baseline.vtu` 等文件。单快照导出不生成新的序列 PVD。

VTK 使用实际顶点与多面体连接，可在 ParaView 中显示体网格与单元场。当前 Tecplot 输出是 FEBRICK，适用于可按其规则表达的八顶点六面体；一般多面体应使用 VTK。后处理会核对元数据、进程目录、全局编号及网格几何，计算后更换网格会导致拒绝导出。

`tools/merge_parallel_csv.py` 面向特定旧 CSV 形式，输出点云；它不替代当前 `babelsim-post` 的几何导出流程。

## 从保存场继续计算

继续计算需要建立新输入算例，明确起始时刻、每个必要场和边界条件。求解命令的 `-time` 不会自动重启。

例如从热传导 `0.03` 继续到 `0.05`：

```bash
python3 docs/examples/create_case.py --template heat --output work/heatRestart
```

编辑 `work/heatRestart/fields/initial/T.field`，保留边界，将初值行改为：

```text
internal file ../../../heat/results/serial/0.03
```

路径从 `fields/initial/` 出发，指向另一个算例的时刻目录。编辑新算例 `control.bs`：

```text
startTime 0.03
endTime 0.05
deltaT 0.01
```

然后运行并比较：

```bash
mpirun -np 2 build/babelsim-solve -case work/heatRestart -time continued
python3 tools/compare_parallel_results.py \
  work/heat/results/serial/0.05 work/heatRestart/results/continued/0.05 \
  --atol 1e-8 --rtol 1e-8
```

本例使用 Euler、相同网格和物性，允许从不同进程数保存的全局场重分区。输入场读取按全局编号映射，使用者必须保证编号对应同一几何单元；仅能读入场文件不证明其网格版本正确。

## 重启范围

流动继续计算需要一致时刻的 `U`、`p` 以及所启用模型的输运场。不要只恢复压力或只恢复速度后声称接续了原有状态。`wallDistance` 等几何输入仍需与网格一致。

当前保存的是场快照，不是完整算法检查点。BDF2 的双时间层历史不由单个结果快照恢复，重新启动会从 Euler 步建立新历史；因此不会保证与不中断的 BDF2 逐步结果完全一致。预条件器、内层迭代历史和全部内部面通量也不作为检查点交付。需要严格完整重启时，应开发相应存储并验证。

## 定量后处理

选择具体数值时刻，按全局编号合并字段，再计算与问题对应的物理量。速度流线不能代替涡量、Q 或旋涡强度；小质量误差不能代替压力/速度耦合误差。频谱分析必须记录采样间隔和可用时长。报告图像应给出物理时刻、单位、色标定义和使用的网格版本。
