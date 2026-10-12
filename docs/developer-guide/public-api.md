# 02 公开 C++ 接口

[开发指南](README.md) · 上一章：[架构](architecture.md) · 下一章：[新增与调用求解器](extending-solvers.md)

## 从数学到代码

BabelSim 区分已知场数学与未知量方程。math 立即计算并返回独立 Field；equ 立即向已绑定未知场的 Equation 加入离散贡献。普通 C++ 负责物理时间、非线性与耦合循环。

```cpp
auto equation = equ::createEquation(problem, spec.equation("temperature"), T);
equation.reset();
equ::ddt(equation, rho * cp, history);
equ::laplacian(equation, conductivity, -1, "diffusion");
equ::source(equation, source);
const auto result = equ::solve(equation);
```

这段代码装配 `rho*cp*dT/dt - div(k*grad(T)) = source`。它所在的完整可编译程序见[external_heat.cpp](../examples/external_heat.cpp)。

## 场的声明、初始化与使用

求解器中的命名场统一由 Case 管理。加载的未知场、程序创建的物性场、通量和压力修正场使用同一种 Field，不需要选择局部或全局存储，也不需要手动注册。所有 `load/create/fields` 入口都返回稳定引用，统一使用 `auto&` 接收。

```cpp
auto& U = problem.load.vectorField("U");
auto& p = problem.load.scalarField("p");
auto& phi = problem.create.faceScalarField("phi");
auto& pPrime = problem.create.homogeneousCorrection(p, "pPrime");
auto& conductivity = problem.create.scalarField("conductivity", 0.1);

// 边界和数值操作直接作用于这个场。
conductivity.fill(0.2);
conductivity.evaluate([](const Vec3& x) { return 0.1 + 0.01 * x.x; });
pPrime.fill(0.0);

// 在其他模型或函数中取得同一个场，不再次读文件或创建。
auto& samePressure = problem.fields.scalarField("p");
problem.output(pPrime); // 选择结果输出，无需再注册
```

| 入口 | 语义 |
| --- | --- |
| `problem.load.scalarField/vectorField/tensorField(name)` | 从配置的 fields 目录加载单元值及边界；重复加载返回同一对象 |
| `problem.create.scalarField/vectorField/tensorField(name, initial)` | 不读文件，创建单元场；初值默认零，默认边界为 zeroGradient |
| `problem.create.faceScalarField/faceVectorField/faceTensorField(name, initial)` | 创建面场；初值默认零，可显式给定 |
| `problem.create.homogeneousCorrection(original, name)` | 创建零单元场；fixedValue 对应 fixedValue(0)，symmetry 保留，其余边界对应 zeroGradient |
| `problem.fields.scalarField/vectorField/tensorField(name)` | 取得已经声明的单元场，不区分其初值来源 |
| `problem.fields.faceScalarField/faceVectorField/faceTensorField(name)` | 取得已经声明的面场 |
| `problem.read.config(description)` | 按成员绑定描述读取、校验参数，返回可补全的普通 C++ 配置结构体 |
| `problem.read.parameters(description)` | 动态工具入口，按描述返回 ParameterValues；通常求解器使用 config |
| `problem.read.physics()/solution()` | 访问 Case 已加载的参数字典 |
| `problem.read.methods()/timeControl()/outputControl()` | 访问已加载的数值方法、时间和输出配置 |

名称在一个 Case 内唯一，创建只能执行一次。重复创建会报错并提示用 `problem.fields` 查找；读取尚未声明的场会报错；不能用 `load` 把程序创建场重新解释成输入文件。修正场也遵守同样规则，并要求原场位于当前 Case 的单元上。

`auto& pPrime = problem.create.homogeneousCorrection(...)` 中的 `&` 是统一的场访问规则。写成 `auto pPrime = ...` 会按 C++ 值语义复制场；项目架构检查会拒绝求解器和示例中这类直接按值接收命名场的写法。这是源码检查，不是对外部 C++ 编译器语义的改变。

## 场的数学操作与输出

```cpp
auto& total = problem.create.scalarField("total");
total = T + C;                        // 更新同一个场，保留它的名称和边界
const ScalarField previous = T;       // 显式保留一份独立数值快照
T.assignScaled(0.5, previous);
T.addScaled(0.5, C);
problem.output(total);                // 选择后续写出
```

Field 持有明确的 Mesh、位置、名称与边界语义。支持整场赋值、缩放、加法和乘系数；赋值要求布局兼容，并保留目标场的物理边界约束。`T + C`、梯度等数学表达式返回的计算结果会随作用域释放，不自动登记为命名场，因此循环中计算表达式不会持续扩充 Case 的场容器。`math::evaluate(operation, target)` 可以复用已有目标存储。

输入单元场遵循现有默认输出选择；程序创建的单元场（包括修正场）通过 `problem.output(field)` 或 `output.bs` 的 `writeFields` 选择输出。面场仍可正常参与数学运算，但当前结果文件格式只支持单元场。输出选择与场的生命周期独立。

公开接口不提供原始数组 `data()`、单元 `operator[]` 或 resize。通过整场函数实现自定义逐值操作，避免把分区/同步细节复制到 Physics。

## 配置与生命周期

`problem.read.physics()`、`problem.read.solution()` 返回只读 Parameters。通过 word、number、positive、nonnegative、fraction、integer、boolean 读取并消费配置。fraction 为 `(0,1]`；带默认值的 getter 仍验证范围。contains/inspect 不替代实际读取，不会让未使用键自动通过。

求解器优先使用 `const auto config = problem.read.config(spec)`。配置结构体只声明成员类型；`spec.bind(parameter, &Config::group, &Group::member)` 同时登记参数说明和成员绑定。默认值、范围、依赖默认值只写在描述中，循环访问 `config.fluid.density` 等普通成员。见[模块描述](module-descriptions.md)。

`equ::createEquation(problem, spec.equation("momentum"), U)` 复用描述中的未知场名称和离散项，并检查实际字段名称、位置和网格。已有直接传名称与离散项的入口仍供没有模块描述的自定义低层算法使用。

完成物理初值、模型和时间历史登记后，显式调用 `problem.finalizeInitialization()`。它检查已登记模块要求的字段、类型/位置、输入来源、修正场来源、方程未知量/离散项、边界能力和未消费配置，成功后关闭新的场、模块与方程配置声明。失败时仍可补齐声明后重试；重复完成调用无副作用。它不推进时间、保存历史、做预测或求解，也不能代替具体模型的物理初值检查。

`time::makeStepper(problem)` 只构造时钟；`advance()/setTime()/write()/finish()` 要求初始化已完成。write/finish 写快照，finish 结束 Case，不判断求解收敛。查询已有字段和配置不推进时间。

## math：显式运算

需要空间重构的函数应传完整 OperatorOptions，常用 `equation.options("diffusion")` 或已声明 operationOptions。不要依赖空选项获得全局默认方法。

| 调用 | 返回与定义 |
| --- | --- |
| grad(scalar, options) | 单元向量梯度 |
| grad(vector, options) | 单元张量梯度 |
| interpolate(cellField, options) | 面场插值 |
| normalGradient(scalar, options) | 外法向面梯度 |
| flux(U, options) | 积分速度面通量 `U_f · S_f` |
| flux(coefficient, scalar, options) | 数学扩散面通量 `coefficient * grad(scalar) · S_f` |
| div(faceScalarFlux) | 单元体积归一的通量散度 |
| div(phi, knownField, options) | 显式对流散度 |
| div(vector, options) | 向量散度 |
| div(tensor, options) | `(div T)_i = dT_ij/dx_j` |
| laplacian(coefficient, scalar, options) | 显式标量拉普拉斯运算 |
| babelsim::map(field, callable) | 对场值执行类型化函数并产生独立结果，位于 babelsim 名字空间 |

标量扩散 math::flux 的符号是正的数学 `k*grad(T)`；物理傅里叶热流使用负号。`math::evaluate(operation, output)` 可填已有场以复用空间；FaceRegion::Interior 是几何内部面集合，包括跨分区内部面，不包括物理边界。

具体重载以[math.h](../../include/babelsim/math.h)为准；标量拉普拉斯与向量/张量运算不能按名称推断为任意类型都支持。

## equ：积分方程

Equation 表示单元积分系统 `A*x=b`，当前显式实例化标量 double 与向量 Vec3。TensorField 的存在不表示已经提供张量未知量线性方程。

| 操作 | 贡献/契约 |
| --- | --- |
| createEquation(Case, name, unknown, terms) | 读取已命名配置并绑定未知场 |
| createEquation(unknown, EquationControl) | 使用已解析的完整控制，不读取字典 |
| reset | 清空项，保留分配、布局与未知量绑定 |
| copy | 独立复制装配系数，仍绑定相同未知场 |
| div(eq, phi, scale, term) | 左端 `scale*div(phi*unknown)`；phi 已乘面面积 |
| laplacian(eq, coefficient, multiplier, term) | 左端 `multiplier*div(coefficient*grad(unknown))` |
| reaction(eq, coefficient) | 左端局部线性反应项，自动积分 |
| source(eq, value[, scale]) | 右端体积源，自动积分 |
| ddt(eq, capacity, history) | 非保守 `capacity*dUnknown/dt`，非 `d(capacity*unknown)/dt` |
| diagonal / rhs | 已积分的对角系数与右端 |
| scale / add / addDiagonal / addRhs | 操作已经积分的系数，不再乘单元体积 |
| apply / residual | 当前装配系统的 `A*x` / `b-A*x` |
| relax | 对角除以 alpha，RHS 加对应旧值贡献 |
| reference / referenceIfUnanchored | 标量参考值/检测常数不确定模式后的参考处理；不验证 Neumann 相容性 |
| solve(eq) | 求解绑定未知量并返回 SolveResult |
| faceFlux(eq, scalar) | 使用装配时冻结的实际扩散项、修正与边界数据 |

所有装配输入在调用时冻结。更改一个系数字段不会自动重新装配旧 Equation。div 使用的 phi 同时成为唯一边界流向上下文；不能在未 reset 的同一方程中换成另一通量而期待覆盖。

常用扩散左端选择 `multiplier=-1`。equ::faceFlux 包含装配 multiplier 的符号，和重新调用 math::flux 不一定等价。守恒压力修正使用实际装配通量，避免复制简化公式丢失非正交修正。

## 时间与历史

```cpp
auto clock = time::makeStepper(problem);
auto history = time::history(T);
while (!clock.finished()) {
    clock.advance();
    history.save(T, clock.dt());
    // 同一物理步的全部装配/求解/校正使用这份 history
}
```

关联场可通过 `time::HistoryGroup` 一起保存并检查同步；`time::weights/predict` 提供公共时间工具，`time::beginStep` 提供显式接受与拒绝。用法与边界见[时间推进与关联历史](time-integration.md)。

TimeStepper 只在 advance 时改变时间，查询 value/dt/step 无副作用。History 只有 save 才推进，必须每物理步一次，不能放进非线性内循环。BDF2 第一次保存历史不足时使用 Euler，此后使用两层历史；显式重构修正还需要求解器正确提供当前迭代/时间预测。

内置 heat 与公开热示例包含相同的 BDF2 场外推，以便新时刻的空间延迟修正一致。自定义耦合示例限定 Euler，使外层时间与内层耦合的契约清楚。

## 诊断、状态与应用入口

`diagnostics::relativeResidual(eq, field)` 按当前已装配系统观察残差，不重新装配；与线性停止阈值使用的缩放不同。`diagnostics::relativeChange` 观察整场变化；其他通量与范数接口见[solver.h](../../include/babelsim/solver.h)。并行归约由框架执行，不需在求解器中增加 MPI_Allreduce。

monitor::Reporter 的 record/iteration 接收名称与值，按进程/间隔输出；它不执行收敛判断。数值失败应保留 `SolveResult.status`，例如 `return SolverResult{result.status}`；外层耗尽返回 notConverged，完成返回 completed。

在文件作用域声明 `SolverRegistration("customHeat", runCustomHeat)`，名称字符串需保持有效且不能重复。main 使用 `runApplication(argc, argv, errorHandler)`；错误回调用于清楚报告配置失败。入口管理应用与 MPI 生命周期，不替求解器自动 finish/write。
