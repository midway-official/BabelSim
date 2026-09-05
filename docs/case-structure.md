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

现有选择是 heat（瞬态热传导）、transport（瞬态对流扩散）、simple（稳态层流 SIMPLE）和
transientSimple（瞬态层流 SIMPLE）。
不再接受 BabelSim 旧的 heatFoam/simpleFoam/transportFoam 名称。
它们不代表 OpenFOAM 程序或输入格式。

普通 Solver 由 Case 取得命名场、物性和算法参数；不用实现自己的 reader。
启动器只调用 runApplication，由框架初始化运行、构造 Case、选择已注册的 Solver 并处理成功/失败。
每个 Solver 在自己的源文件注册名称和函数；不修改启动器或 Case reader。

### solver 名称如何对应 C++ 函数

`case.bs` 的 `solver heat` 只保存字符串 `heat`，不会自动查找名为 `runHeat` 的函数。
对应关系在求解器自己的 `src/physics/heat/main.cpp` 中声明：

```cpp
const SolverRegistration heat("heat", runHeat);
```

这一行放在 `namespace babelsim` 中、`runHeat` 函数外，并包含 `babelsim/application.h`。
SIMPLE 与 Transport 分别只注册自己。通用启动器不再维护对应表：

```cpp
#include "babelsim/application.h"
int main(int argc, char* argv[]) {
    return babelsim::runApplication(argc, argv);
}
```

`readCase()` 把配置名称交给 `Case::solver()`；`runApplication()` 比较已注册名称，
找到 `heat` 后调用对应的 `runHeat(problem)`。函数必须在构建时链接到可执行程序；
这不是动态插件、文件名查找或自动命名规则。没有注册、未知名称、重名或空函数会报错。
注册对象只记录名称与函数，启动前不调用 MPI；实际检查与分派仍由框架完成。
外部 Solver 使用相同的一行注册和通用 main，无需维护另一套表。

注册只负责把 Case 中的名称连接到求解函数，不会把具体 Solver 变成公共类。
`include/babelsim/` 不提供 Heat、Transport 或 SIMPLE 头文件；Case 用户只选择名称，
Solver 作者只依赖公共 Framework API。内置求解器的实现全部留在各自的 `src/physics/` 目录。

### physics 条目与 physics() 接口

`problem.physics()` 返回 `case.bs` 的 `physics` 条目所指向的参数字典。
例如 `physics physics/thermal.bs` 对应：

```cpp
const double rho = problem.physics().positive("density");
const double k = problem.physics().nonnegative("conductivity");
```

它与 `problem.solution()` 的命名规则一致：接口直接对应配置条目。路径可以更换，
无需修改读取参数的 Solver。字典在 Case 构造时读取，调用接口时复用同一对象，
仍保留数值检查、稳定引用和未使用参数检查。
旧 `Case::properties()` 已改为 `Case::physics()`，不保留别名；外部 Solver 需更新调用并重新编译。
`case.bs` 和各物理参数文件的格式、键名不变，`Parameters` 通用字典类型也不变。

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
convection linearUpwind
diffusion orthogonal
time euler
convection C upwind
```

三列条目覆盖某个 Field 的默认方法。对流可选一阶 `upwind`、二阶
`linearUpwind` 和 `central`；`linearUpwind` 使用迎风单元梯度重构面值，隐式部分仍保留
一阶迎风基底。非正交问题使用 corrected 插值/扩散、适当梯度；Solver 不手写几何修正。

`numerics/solution.bs`：

```text
scalarSolver bicgstab ilut 1e-14 1e-10 1000
vectorSolver bicgstab ilut 1e-12 1e-8 1000
```

每行在配置名后依次填写方法、预条件器、绝对容差、相对容差、最大迭代数。
scalarSolver 配置标量方程，vectorSolver 配置矢量方程，两个条目均为必填，不能省略。
缺少任意一个都会在 Case 创建时报告 solution.bs 路径和缺失键名，不再回退到默认值。
即使 Heat/Transport 当前只求标量，也必须提供两项；配置矢量求解器不会创建或求解额外矢量方程。
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
瞬态 transientSimple 使用 `problem.loop()` 推进物理时间，并在每个时间步内部完成 SIMPLE
压力速度校正；其 `methods.bs` 必须选择 `euler` 或 `bdf2`，不能选择 `steady`。

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
metadata 的 time 记录真实物理时间。当前不支持 checkpoint/restart。
可用 problem.output(derived) 选择派生 cell 场输出，或 output(input,false) 关闭输入场输出；
场必须属于当前 Case，face 场暂不支持保存。各进程必须以相同逻辑修改输出选择，建议在声明阶段完成。

Case::validate() 只校验参数；start() 和首次 loop() 才关闭声明阶段。
SIMPLE 构造不再隐式关闭声明，便于组合热/流体等算法。输入和中间场都在开始计算前创建，
scalar/vector/tensor 中间场可指定初值；面场支持 scalar/vector/tensor 三种值类型。
非均匀初值可在 Solver 中用 Field::evaluate(位置函数) 定义，但文件 reader 仍只接受 uniform。

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
