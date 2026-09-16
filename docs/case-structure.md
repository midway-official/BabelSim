# BabelSim Case 配置与场生命周期

Case 描述问题数据和运行控制；Solver 描述 PDE 和算法。配置解析、场所有权、时间元数据和结果
写出由 Case/运行时提供，Physics 只用类型安全的只读接口。

## 1. 目录

~~~text
cases/mySolver/
├── case.bs
├── mesh/mesh.mesh
├── fields/initial/T.field
├── physics/physics.bs
├── numerics/methods.bs
├── numerics/solution.bs
├── control.bs
└── output.bs
~~~

case.bs 必须为每个路径提供一个值：

~~~text
solver mySolver
mesh mesh/mesh.mesh
fields fields/initial
physics physics/physics.bs
methods numerics/methods.bs
solution numerics/solution.bs
control control.bs
output output.bs
~~~

路径相对于案例根目录。未知键、重复键、绝对输出路径和不完整 case 会在构造时拒绝。
ghostLayers（若提供）必须是至少 3 的整数。

## 2. 配置职责

### physics/physics.bs

放物性、模型选择和模型常数，例如：

~~~text
density 1.0
dynamicViscosity 0.001
turbulenceModel none
~~~

### numerics/methods.bs

离散格式只在这里定义，并由 Case 构造时读取一次：

~~~text
interpolation linear
gradient leastSquares
convection upwind
diffusion corrected
time euler
gradient T greenGauss
diffusion p orthogonal
~~~

前五项是默认格式；带字段名的行是覆盖。运行时通过 problem.methods() 提供只读的
Methods，不允许 Solver 自己重新读取或改变它。

### numerics/solution.bs

这里放线性系统、算法迭代和模型数值控制：

~~~text
scalarSolver bicgstab ilut 1e-14 1e-10 1000
vectorSolver bicgstab ilut 1e-12 1e-8 1000
maxIterations 1000
velocityRelaxation 0.7
pressureRelaxation 0.3
continuityTolerance 1e-8
~~~

单个未知量可以覆盖默认线性配置：

~~~text
scalarSolver.T cg incompleteCholesky 1e-13 2e-9 321
vectorSolver.U bicgstab ilut 1e-13 1e-9 500
~~~

通过 readLinearControl(problem, field) 读取；调用者不接触 token。

### control.bs

~~~text
startTime 0
endTime 1
deltaT 0.01
~~~

TimeStepper 校验有限值、正时间步和结束时间。BDF2 要求等步长且时间区间包含整数步。

### output.bs

~~~text
directory results
timeName final
writeInterval 5
writeFields T U derived
excludeFields U
~~~

writeFields 和 excludeFields 是可选列表。默认策略是写出文件加载的 cell 场；程序创建的
cell 场必须列入 writeFields 或由 problem.output(field) 显式选入。Face 场当前用于离散中间
量，不能直接写到结果格式。字段名称在实际 write 时确认已经声明。

## 3. 参数读取和有效配置

~~~cpp
const auto& physics = problem.physics();
const double rho = physics.positive("density");
const double nu = physics.nonnegative("dynamicViscosity", 0.0);
const std::string model = physics.word("turbulenceModel", "none");

const auto& solution = problem.solution();
const bool clip = solution.boolean("clipTurbulence", false);
const int maxIter = solution.integer("maxIterations", 1000, 1, 1000000);
~~~

Parameters 统一处理缺失值、默认值、类型、有限值和范围。positive、nonnegative、fraction
和 integer(min,max) 在读取点完成校验。inspect(key) 返回 configured、consumed 和行号，
用于调试生效配置。problem.validate() 会调用 physics/solution 的 requireAllUsed；它只校验，
不推进时间，不写结果，不关闭声明阶段。

## 4. 场的三种明确操作

~~~cpp
auto& T = problem.scalarField("T");             // 读取 fields/initial/T.field
auto& U = problem.vectorField("U");             // 读取 U.field
auto& phi = problem.createFaceField("phi");     // 创建面场，不读文件
auto& k = problem.createScalarField("k", 0.0);  // 创建 cell 场
auto& same = problem.existingScalarField("k");  // 查找已声明场
~~~

命名场由 Case 绑定当前 Mesh 和边界。加载 API 只加载 cell 文件；create API 只创建程序
场。第一次 create 的初值生效，重复 create 或 load/create 混用会报错，不会静默忽略。
声明阶段结束后不能增加 Case 场；局部 math 返回值可继续使用。

压力修正等齐次边界场使用：

~~~cpp
auto pPrime = field::homogeneousLike(p);
~~~

它为给定 cell 场生成同网格、同位置、对应齐次边界类型的场，不修改原场。函数名应由
调用者给出清楚的物理含义。

## 5. 时间和输出生命周期

~~~cpp
auto time = time::start(problem);
auto history = time::history(T);
problem.validate();

while (time.value() < time.end()) {
    time.advance();
    history.save(T, time.dt());
    // 组装、求解、收敛判断
    if (time.step() % readWriteInterval(problem) == 0 || time.finished())
        write(problem, time);
}
~~~

TimeStepper 只负责推进和末步截断；History 由 Solver 明确保存；运行时不判断物理收敛、不
打印。稳态 Solver 使用自己的 for 循环，不调用 Case::loop。

## 6. 结果结构

每个写出时间目录包含按 rank 分区的 cell 场文件和 metadata。使用不同进程数重跑同一结果
目录会被拒绝，应用应通过 run name 隔离实验。Case::write() 只写已经选择的字段，不宣称
求解成功；SolverResult 的成功与否由 application 处理。
