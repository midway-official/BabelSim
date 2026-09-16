# 用 BabelSim DSL 开发物理求解器

本文是当前 BabelSim Physics SDK 的唯一入门指南。它描述已经存在于
`include/babelsim/` 和 `src/physics/` 的接口；代码片段按当前源码组织，可以直接作为新
Solver 的起点。底层 MPI、halo、稀疏矩阵、Eigen、CUDA 和网格遍历均由运行时与离散后端
实现，Physics 只表达场、方程和算法。

## 1. 先建立正确的心智模型

一个 Solver 是普通的 C++ 函数：

```cpp
babelsim::SolverResult runMySolver(babelsim::Case& problem);
```

函数中按教材顺序写出：

1. 从 `Case` 加载初始场；
2. 从分类字典读取物性和数值控制；
3. 创建程序拥有的中间场和绑定未知量的 `equ::Equation`；
4. 在循环中保存历史、清空方程、逐项离散；
5. 求解并显式修正其它场或面通量；
6. 用 `diagnostics` 观察当前方程和场的状态；
7. 用通用 `monitor::Reporter` 报告观测值；
8. 返回 `SolverResult`，由应用层决定进程退出码。

Solver 不创建解析器、Manager、Factory 或公共基类，也不调用 MPI 或打印库。稳态 SIMPLE、
瞬态 SIMPLE、Heat、Transport 仍是四个独立的 Solver；两个 SIMPLE 的代码故意分别放在
自己的 `main.cpp`，Rhie–Chow 组合也只写在对应的 `main.cpp`。

## 2. Case 与配置

### 2.1 目录和 `case.bs`

```text
cases/mySolver/
├── case.bs
├── mesh/mesh.mesh
├── fields/initial/T.field
├── physics/physics.bs
├── numerics/methods.bs
├── numerics/solution.bs
├── control.bs
└── output.bs
```

`case.bs` 只声明路径和 Solver 名称：

```text
solver mySolver
mesh mesh/mesh.mesh
fields fields/initial
physics physics/physics.bs
methods numerics/methods.bs
solution numerics/solution.bs
control control.bs
output output.bs
```

Case 构造时读取网格、`physics`、`solution`、`methods` 和输出控制。方法字典只解析一次；
`problem.methods()` 是只读生效配置，不应在 Solver 中再次调用解析器。

### 2.2 配置分类

| 文件 | 内容 | 典型键 |
| --- | --- | --- |
| `physics/*.bs` | 模型和物性 | `density`、`dynamicViscosity`、`turbulenceModel` |
| `numerics/methods.bs` | 空间、时间离散方法 | `gradient`、`convection`、`diffusion`、`time` |
| `numerics/solution.bs` | 线性系统和算法迭代控制 | `scalarSolver`、`maxIterations`、`velocityRelaxation` |
| `control.bs` | 物理时间区间 | `startTime`、`endTime`、`deltaT` |
| `output.bs` | 结果目录、时间名和字段筛选 | `directory`、`writeInterval`、`writeFields` |

同一个键只在所属字典消费。拼写错误、重复键、缺少必填项和未消费项都会在
`problem.validate()` 或 Case 构造时报告，Physics 不接触原始 token。

### 2.3 类型安全的参数读取

`Parameters` 提供带校验的接口：

```cpp
const auto& physical = problem.physics();
const double rho = physical.positive("density");
const double cp  = physical.positive("heatCapacity");
const double q   = physical.number("source", 0.0);
const bool clipped = problem.solution().boolean("clipTurbulence", false);
const int maxIterations =
    problem.solution().integer("maxIterations", 1000, 1, 1000000);
const std::string model =
    physical.word("turbulenceModel", "none");
```

`positive`、`nonnegative`、`fraction`、整数范围和有限值检查在参数层完成。不要在
Physics 中调用 `entry().tokens` 或把字符串转换成数字。需要诊断生效状态时使用：

```cpp
const ParameterInfo info = problem.solution().inspect("maxIterations");
// info.configured、info.consumed、info.line
```

在声明阶段完成所有场、物性和设置读取后调用一次 `problem.validate()`。它只校验，不推进
时间、不写文件，也不会提前禁止继续声明场。第一次进入 `time::start` 或
`setTime` 时 Case 才锁定声明阶段。

## 3. 场的加载、创建和查找

Case 是当前问题的所有者，因此命名场使用 `problem.*Field`。这不是隐藏的求解器，而是
保证场和同一个 Mesh、边界和输出生命周期绑定。

### 3.1 三种动作必须写清楚

```cpp
// 读取 fields/initial/T.field；同名重复读取返回同一个稳定引用
auto& T = problem.scalarField("T");
auto& U = problem.vectorField("U");

// 由程序创建，不读取文件；初值在第一次 create 调用时生效
auto& phi = problem.createFaceScalarField("phi");
auto& correction = problem.createScalarField("correction", 0.0);

// 访问之前已经声明的场，不加载、不创建第二份对象
auto& correctionAgain = problem.existingScalarField("correction");
```

- `scalarField/vectorField/tensorField(name)` 只表示从初始场目录加载单元场。
- `createScalarField/createVectorField/createTensorField` 只创建程序场；重复 create 会报错，
  不会静默忽略新的初始化值。
- `createFaceScalarField/createFaceVectorField/createFaceTensorField` 创建面场，不能从结果输出
  配置直接写出。
- `existing*Field` 只查找已经声明的对象，用于模型之间明确共享一个场。
- 声明阶段结束后不能创建新 Case 场；局部派生量仍可用值语义的 `math::*` 返回值创建。

压力修正场的边界必须明确表达其含义：

```cpp
auto pPrime = field::homogeneousLike(p, "pPrime");
```

它复制 `p` 的网格和位置，把定值边界变成齐次定值、定梯度边界变成零梯度；名称比
`correction` 更准确的 API 仍由调用者给出。这个操作只准备边界条件，不执行压力修正。

### 3.2 输出选择

文件加载的单元场默认写出。程序创建的场默认不写出，可在 `output.bs` 中声明：

```text
directory results
timeName final
writeInterval 5
writeFields T U wallShear
excludeFields U
```

`writeFields` 和 `excludeFields` 中的名称在实际 `write()` 时验证必须是已声明的单元场。
求解器仍可对已拥有的场显式选择：

```cpp
problem.output(wallShear);        // 加入后续写出
problem.output(U, false);         // 关闭某个场
```

显式选择只改变输出列表，不立即写文件。默认输出策略和配置筛选因此与方程代码分离。

## 4. math 与 equ：一个名字，一种语义

### 4.1 `math` 是已经计算的场数学

`math` 不维护生命周期，也不保存历史。常用操作：

```cpp
const auto gradT = math::grad(T);                // cell vector field
const auto gradU = math::grad(U);                // cell tensor field
const auto faceT = math::interpolate(T);         // face scalar field
const auto phi = math::flux(U);                  // face flux field
const auto divPhi = math::div(phi);              // explicit cell divergence
const auto divStress = math::div(stress);        // tensor -> vector divergence
const auto diffusionFlux = math::flux(k, T);     // explicit face diffusion flux
const auto lapT = math::laplacian(T);            // explicit cell Laplacian
```

`math::interpolate` 永远表示插值，不偷偷执行 Rhie–Chow。SIMPLE 的 Rhie–Chow 由求解器
main 直接用 `math::flux`、`math::interpolate`、`math::grad` 和面区域更新组合出来。
`math::div(phi)` 是显式场散度；它不会从某个 Equation 猜测未知量。

场代数也属于 `math`：

```cpp
const auto stress = muEff *
    (math::transpose(gradU) -
     (2.0 / 3.0) * math::isotropic(math::trace(gradU)));
U -= rAU * math::grad(pPrime);
```

`Field` 自身的 `fill`、赋值、`+=`、`-=` 和缩放只修改目标场，不负责通信或方程
装配。

### 4.2 `equ` 是立即加入 Equation 的离散项

```cpp
auto temperatureEquation = equ::createEquation(T);
temperatureEquation.reset();
equ::ddt(temperatureEquation, rho * cp, history);
equ::laplacian(temperatureEquation, k, -1);
equ::source(temperatureEquation, Q);
const SolveResult result =
    equ::solve(temperatureEquation, linearControl);
```

`createEquation(T)` 明确创建“绑定未知量 T 的离散方程”。它不会从字段名推断 PDE；
后续每个 `equ::*` 调用立即把一项加入同一个线性系统。约定如下：

- `equ::ddt(eq, capacity, history)` 加入容量乘时间导数；
- `equ::div(eq, phi)` 对绑定未知量加入隐式 \\(\\nabla\\cdot(\\phi x)\\)；
- `equ::laplacian(eq, gamma, multiplier)` 加入
  `multiplier * div(gamma * grad(x))`；`-1` 是通常的扩散左端符号；
- `equ::reaction` 加入隐式反应对角项；
- `equ::source` 只加入已知右端，不进行隐式线性化；
- `equ::relax(eq, previous, alpha)` 使用当前 BabelSim 的对角和右端松弛约定；
- `equ::faceFlux(eq, solution)` 返回与该 Equation 实际离散项一致的标量面通量。

`equ::apply`、`equ::residual` 和 `diagnostics::relativeResidual` 观察已组装的系统，
不会再次组装。Equation 的 `reset`、`diagonal`、`rhs`、`referenceIfUnanchored` 是状态或
结构操作，属于 Equation 对象本身。需要 SIMPLE 响应系数时，几何体积和对角系数保持
独立：`const auto rAU = geometry::cellVolumes(problem.mesh()) / equation.diagonal();`。

## 4.1 几何量也是可组合的场

`geometry` 只把已经加载网格中的几何数组 materialize 成值场；它不读取文件、不注册输出，
也不包含物理算法。静态网格的几何量应在循环外取得：

```cpp
const auto V  = geometry::cellVolumes(problem.mesh());
const auto Sf = geometry::faceAreaVectors(problem.mesh());
const auto Af = geometry::faceAreas(problem.mesh());
const auto C  = geometry::cellCentres(problem.mesh());
```

`V` 是控制体积，`Sf` 是与通量和散度使用同一拓扑方向的面面积向量，`Af` 是面面积，
`C` 是单元中心；还可使用 `faceCentres` 和 `faceUnitNormals`。它们返回普通 Field，
可以直接参与 `math::interpolate`、场代数和 `math::dot`。

## 5. 时间、历史和循环

只有瞬态算法需要时间服务：

```cpp
auto time = time::start(problem);
auto T_old = time::history(T);
problem.validate();

while (time.value() < time.end()) {
    time.advance();               // 截断最后一步并更新 Case 时间元数据
    T_old.save(T, time.dt());     // 每个物理步一次；内迭代不能调用
    // 组装、求解、检查
    if (time.step() % writeInterval == 0 || time.finished())
        write(problem, time);
}
```

`TimeStepper::advance()` 是唯一推进时间的调用；它处理最后一步和浮点边界。运行时不
替 Solver 写历史、不判断物理收敛、不打印。Euler/BDF2 的选择来自一次性加载的
`problem.methods().time`，`equ::ddt(eq, capacity, history)` 使用显式历史，并在
BDF2 首步自动使用 Euler。稳态 Solver 不应创建时间循环；它检查方法为
`TimeMethod::Steady` 后直接进入自己的迭代循环。

Solver 使用显式 `TimeStepper` 或自己的迭代循环；Case 不提供隐藏的 `problem.loop()`，
因此代码中清楚可见何时保存历史、何时组装和何时写出。

## 6. 线性求解、诊断和通用监视器

线性控制按未知场读取，签名和行为一致：

```cpp
auto& T = problem.scalarField("T");
const auto linear = readLinearControl(problem, T);

auto& U = problem.vectorField("U");
const auto vectorLinear = readLinearControl(problem, U);
```

`solution.bs` 的默认项是：

```text
scalarSolver bicgstab ilut 1e-14 1e-10 1000
vectorSolver bicgstab ilut 1e-12 1e-8 1000
```

对单个场可覆盖：

```text
scalarSolver.T cg incompleteCholesky 1e-13 2e-9 321
```

Solver 只看到 `LinearSolverConfig` 和 `SolveResult`，看不到 CSR 或预条件器对象。
`SolveResult::converged()` 是该线性系统是否达到后端容差；Physics 自己还要判断外层
算法的物理残差、场变化和守恒。

诊断只观察，不修改：

```cpp
const double r = diagnostics::relativeResidual(eq, T);
const double dT = diagnostics::relativeChange(T, T_previous);
const auto mass = diagnostics::fluxBalance(phi);
const bool finite = diagnostics::all(std::isfinite(r) && std::isfinite(dT));
```

监视器是通用的名称—数值报告器：

```cpp
const monitor::Reporter report("mySolver");
report.iteration(iter + 1, maxIterations, converged, {
    {"residual", r}, {"change", dT}, {"mass", mass.relative},
    {"converged", converged}
});
```

Reporter 只负责格式、周期和输出进程选择；它不理解 SIMPLE、RANS 或任何收敛准则。收敛
布尔值由 Solver 根据自己的算法计算，并作为普通 metric 传给 Reporter。Physics 返回：

```cpp
return SolverResult::completed();
return SolverResult::notConverged();
return SolverResult::numericalFailure();
```

应用层才把状态映射成命令行退出码。

## 7. 最小 Heat Solver

下面的完整结构与仓库中的 Heat 相同，故意没有公共 Heat 基类：

```cpp
SolverResult runHeat(Case& problem) {
    const monitor::Reporter report("heat");
    auto& T = problem.scalarField("T");

    const auto& physics = problem.physics();
    const double rho = physics.positive("density");
    const double cp  = physics.positive("heatCapacity");
    const double k   = physics.nonnegative("conductivity");
    const double Q   = physics.number("source", 0.0);

    const auto linear = readLinearControl(problem, T);
    const int writeInterval = readWriteInterval(problem);
    auto time = time::start(problem);
    auto history = time::history(T);
    auto equation = equ::createEquation(T);
    problem.validate();

    while (time.value() < time.end()) {
        time.advance();
        history.save(T, time.dt());

        equation.reset();
        equ::ddt(equation, rho * cp, history);
        equ::laplacian(equation, k, -1);
        equ::source(equation, Q);

        const auto solved = equ::solve(equation, linear);
        report.record({{"time", time.value()},
                       {"residual", solved.relative_residual}});
        if (!solved.converged())
            return SolverResult{solved.status};

        if (time.step() % writeInterval == 0 || time.finished())
            write(problem, time);
    }
    return SolverResult::completed();
}
```

代码逐项对应 \\(\\rho c_p \\partial_t T = \\nabla\\cdot(k\\nabla T)+Q\\)。场加载和
文件输出由 Case 提供，时间推进和历史由 time 提供，方程项由 equ 提供，算法仍由这个
Solver 函数掌握。

## 8. 新 Transport Solver 的最小结构

```cpp
auto& C = problem.scalarField("C");
auto& U = problem.vectorField("U");
auto& phi = problem.createFaceScalarField("phi");
phi = math::flux(U);

const auto& physics = problem.physics();
const double storage = physics.positive("storage");
const double D = physics.nonnegative("diffusivity");
const double source = physics.number("source", 0.0);

auto time = time::start(problem);
auto history = time::history(C);
auto equation = equ::createEquation(C);
problem.validate();

while (time.value() < time.end()) {
    time.advance();
    history.save(C, time.dt());

    equation.reset();
    equ::ddt(equation, storage, history);
    equ::div(equation, phi);
    equ::laplacian(equation, D, -1);
    equ::source(equation, source);

    const auto solved = equ::solve(equation, readLinearControl(problem, C));
    if (!solved.converged()) return SolverResult{solved.status};
}
```

这里的 `phi` 是程序创建的面通量，不会误读同名初始文件；如果需要输出派生场，再调用
`problem.output(field)` 或在 `output.bs` 中列出它。

## 9. SIMPLE 的组织规则

稳态和瞬态 SIMPLE 不共享实现。二者在各自 main 中都应按以下顺序展开，但变量、Equation、
循环和修正代码分别维护：

**稳态 SIMPLE：一个外层循环**

```cpp
for (int iter = 0; iter < maxIterations; ++iter) {
    const VectorField U_previous = U;

    momentumEquation.reset();
    equ::div(momentumEquation, phi, rho);
    equ::laplacian(momentumEquation, muEff, -1);
    equ::source(momentumEquation, -math::grad(p));

    const double rU = diagnostics::relativeResidual(momentumEquation, U);
    // relax -> response -> velocity solve

    // 在本 main 中显式组合 Rhie–Chow：
    // math::grad / math::interpolate / math::flux / math::add/subtract

    pressureCorrectionEquation.reset();
    equ::laplacian(pressureCorrectionEquation, rAU, -1);
    equ::source(pressureCorrectionEquation, -math::div(phiH));
    pressureCorrectionEquation.referenceIfUnanchored(0.0);
    // pressure solve -> p/U/phi correction

    // optional RANS transport, diagnostics and convergence decision
}
```

允许压力方程内部保留显式 non-orthogonal 修正子循环；那是数值修正循环，不改变稳态
只有一个主要外迭代的结构。

**瞬态 SIMPLE：两层主要循环**

```cpp
auto time = time::start(problem);
auto UHistory = time::history(U);
while (time.value() < time.end()) {
    time.advance();
    UHistory.save(U, time.dt());

    for (int iter = 0; iter < maxIterations; ++iter) {
        // ddt + momentum assembly
        // residual before relaxation
        // velocity solve
        // this main's explicit Rhie–Chow
        // pressure correction and optional non-orthogonal sub-loop
        // turbulence, diagnostics and convergence
        if (converged) break;
    }
    if (!converged) return SolverResult::notConverged();
}
```

不要把上面的流程替换成 `simple.solve()` 或 `problem.loop()` 黑盒。SIMPLE 专用
Rhie–Chow、压力/速度/通量修正、松弛约定和收敛判定必须在对应 main 中可见。

## 10. 独立的 RANS 模型

RANS 模型不是 SIMPLE 的隐藏公式，也不共享其它模型的闭合实现。模型接口只暴露物理对象：

```cpp
auto turbulence = rans::load(problem, U, phi);
const auto& muEff = turbulence.effectiveViscosity();

equ::laplacian(momentumEquation, muEff, -1);
equ::source(momentumEquation,
            math::div(turbulence.deviatoricStressRemainder(U)));

const auto transport = turbulence.solveTransport();
```

Spalart–Allmaras、k–omega、k–epsilon 各自拥有需要的场、输运方程、边界、历史和闭合
常数。SIMPLE 只知道有效黏度、应力余项和带名字的输运诊断，不知道某个模型内部如何
计算源项。模型的 `turbulenceRelaxation`、`turbulenceTolerance` 等数值控制位于
`solution.bs`，模型常数位于 `physics.bs`。

## 11. 注册和编译

每个 Solver 在自己的 `main.cpp` 注册：

```cpp
namespace babelsim {
SolverResult runMySolver(Case&);
const SolverRegistration mySolver("mySolver", runMySolver);
}
```

应用入口只调用通用 `runApplication(argc, argv)`。Physics 源码只包含公共头和自己目录
的私有头；不得包含 `src/internal`、其它 Solver 目录、MPI、Eigen、CSR 或原始 Field
访问。Makefile 把该文件加入对应构建目标即可，不需要修改集中式 Solver 表。

## 12. 开发检查清单

- 场加载是否使用 `scalarField/vectorField/tensorField`，程序场是否使用 `create*Field`？
- 复用已声明场是否使用 `existing*Field`，并且没有依赖调用顺序的静默初始化？
- 是否在所有读取/创建/边界定义后调用一次 `problem.validate()`？
- 方法是否只通过 `problem.methods()` 查询，线性控制是否用具体未知场读取？
- `math` 是否只做显式场数学，`equ` 是否只向绑定 Equation 加离散项？
- `equ::laplacian` 的系数和 `-1` 符号是否在调用点清楚可见？
- 时间历史是否每个物理步只保存一次，内迭代是否没有推进时间？
- 是否在重新组装前检查 `diagnostics::relativeResidual`，避免无意义的重复组装？
- 收敛判定是否仍在 Solver 中，Reporter 是否只是通用报告？
- 输出字段是否通过 `output.bs` 选择，而不是为了写出而修改方程代码？
- SIMPLE 的 Rhie–Chow 是否仍只位于对应 main，稳态/瞬态是否没有共享算法实现？
- 是否运行 `make -j4 test`、工作流和需要的 MPI/外部 Solver 检查？

相关参考：

- [DSL API 语义](procedural-dsl.md)
- [Case 配置](case-structure.md)
- [架构与维护边界](architecture.md)
- [验证与回归](validation.md)
