# Case：问题配置、运行与结果

Case 描述“问题是什么”，Solver 描述“方程和算法是什么”。BabelSim 保留初值、物性、
数值控制分离的思想，不照搬 OpenFOAM 的全部字典结构。

```text
cases/heat/
├── case.bs
├── mesh/heat.mesh
├── fields/initial/T.field
├── physics/thermal.bs
├── numerics/methods.bs
├── numerics/solution.bs
├── control.bs
└── output.bs
```

## 入口与名称

```text
solver heat
mesh mesh/heat.mesh
fields fields/initial
physics physics/thermal.bs
methods numerics/methods.bs
solution numerics/solution.bs
control control.bs
output output.bs
```

现有选择是 heat（热传导）、simple（稳态层流 SIMPLE）、transport（对流扩散）。
不再接受 BabelSim 旧的 heatFoam/simpleFoam/transportFoam 名称。
它们不代表 OpenFOAM 程序或输入格式。

普通 Solver 由 Case 取得命名场、物性和算法参数；不用实现自己的 reader。
启动器只初始化运行、构造 Case、调用显式选择表并处理成功/失败。

## 网格和场

网格仍使用 BABELSIM_MESH 1，存尺寸、笛卡尔边界或显式顶点、patch。
统一三维结构化六面体，二维是 nz=1。各 patch 名必须唯一。

实际可读取的场语法是花括号形式，不是旧文档中的简写：

```text
field T
{
    type scalar
    location cell
    internal uniform (0)
    boundary
    {
        hot   { type fixedValue value (1) }
        cold  { type fixedValue value (0) }
        lower { type symmetry }
        upper { type symmetry }
        front { type symmetry }
        back  { type symmetry }
    }
}
```

每个物理 patch 都必须配置。标量值也写括号；vector 是三个分量，tensor 是行优先九个分量。
目前读取支持 cell 场和 uniform 初值；面通量由数学算子生成，不读取 face 初值文件。
支持 fixedValue、fixedGradient、zeroGradient、inletOutlet、symmetry/mirror。

## 物理与数值分开

`physics/thermal.bs`：

```text
density 1
heatCapacity 1
conductivity 0.1
source 0
```

`numerics/methods.bs`：

```text
interpolation linear
gradient greenGauss
convection upwind
diffusion orthogonal
time euler
convection C upwind
```

三列条目覆盖某个 Field 的默认方法。非正交问题使用 corrected 插值/扩散、适当梯度；
Solver 不手写几何修正。

`numerics/solution.bs`：

```text
scalarSolver bicgstab ilut 1e-14 1e-10 1000
```

依次是方法、预条件器、绝对容差、相对容差、最大迭代数。
vectorSolver 配置矢量方程；未给出时使用 RuntimeControl 默认值。
SIMPLE 的压力使用 scalarSolver、速度使用 vectorSolver；旧 velocitySolver/pressureSolver
需改为这两个通用键。SIMPLE 自身的松弛、最大外迭代和容差仍在此文件，由算法读取。

`control.bs`：

```text
startTime 0
endTime 0.05
deltaT 0.01
```

Euler 的最后一步可缩短到 endTime，不越过终点。BDF2 首步 Euler、后续等步长；
不支持以非整数步数终止的 BDF2。稳态 simple 使用算法外循环，不把迭代次数冒充物理时间。

## 自动输出

```text
directory results
timeName final
writeInterval 2
```

writeInterval 是正整数步数，省略时为 1。上例时间设置输出 0.02、0.04、0.05，
末时刻不因未落在间隔上而丢失。timeName 只是最终结果别名，不是伪造的物理时间。

```text
results/
├── 0.02/rank-0000/T.csv + metadata.bs
├── 0.04/rank-0000/...
├── 0.05/rank-0000/...
└── final/rank-0000/...      最终状态的兼容入口
```

每个 rank 只写 owned cell。输入场自动输出，中间数学场和面通量不自动输出。
metadata 的 time 记录真实物理时间。当前不支持 checkpoint/restart 或计算中动态改变输出场集合。

为独立实验指定名称：

```bash
mpirun -np 4 build/babelsim-solve -case cases/heat -time mpi4
```

此时序列放在 results/mpi4/<物理时间>/，最终状态仍可从 results/mpi4/rank-*/ 读取。
不同进程数、网格和参数的实验请使用不同名称或 output.directory，避免旧时间目录混入新实验；
当前不会自动清理已有结果，也不保证同一路径并发写入安全。
不同 rank 数量的旧结果会在写入前被拒绝；请换一个标签。标签和 timeName 不可为数字、all 或 latest，
以免覆盖物理时间目录或与后处理选择冲突。-time 是运行标签，不是重启时刻。

## 独立后处理

```bash
# 默认序列
build/babelsim-post -case cases/heat -time all -format vtk tecplot
# 命名运行的序列
build/babelsim-post -case cases/heat -time mpi4/all -format vtk tecplot
# 最晚物理时间，按数值而不是字典序选择
build/babelsim-post -case cases/heat -time mpi4/latest -format vtk
# 只处理最终别名
build/babelsim-post -case cases/heat -time mpi4 -format vtk tecplot
```

-format vtk 现在生成 XML UnstructuredGrid 的 .vtu，PVD 引用 .vtu；不再用 legacy .vtk
作为时间集合的子文件。ParaView 打开 post/series.pvd 或 post/mpi4/series.pvd。
Tecplot 仍输出 FEBRICK .dat，每个时刻独立文件。

all/latest 只扫描数值目录，忽略 final 和实验标签；all 遇到缺 rank、缺全局 ID、
网格不匹配或元数据时间不一致会报错，不能静默跳过坏时间步并生成看似完整的动画。
旧手工保存的非数值标签仍可用 -time <标签> 单独处理。

这条工作流由 make test-workflow 验证，包括 1/2/4 进程、时间间隔、短末步、
数值排序和可用时的真实 ParaView PVDReader。
