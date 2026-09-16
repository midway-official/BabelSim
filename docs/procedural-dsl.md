# BabelSim 过程式 PDE DSL 参考

本文是当前 Physics SDK 的语义速查；完整的新求解器示例见 [DSL 求解器开发指南](dsl-solver-guide.md)。
BabelSim 使用 C++ embedded DSL，不定义另一门脚本语言；普通 C++ 的 if、for、while 和 function
就是算法控制流。

## 名字空间职责

| API | 语义 | 典型接口 |
| --- | --- | --- |
| Case / Field | 网格上下文、命名场、边界、输出和时间元数据 | scalarField、createFaceScalarField、write |
| math | 已经计算的场数学 | grad、div、interpolate、flux、laplacian |
| equ | 对绑定未知量立即加入离散项 | createEquation、ddt、div、laplacian、source |
| time | 时间推进和显式历史 | start、TimeStepper::advance、History::save |
| diagnostics | 只读观测 | relativeResidual、relativeChange、fluxBalance |
| monitor | 通用报告格式 | Reporter::record、Reporter::iteration |

## 场生命周期

~~~cpp
auto& T = problem.scalarField("T");             // 加载 fields/initial/T.field
auto& U = problem.vectorField("U");             // 加载 fields/initial/U.field
auto& phi = problem.createFaceScalarField("phi");     // 创建面场，不读文件
auto& rAU = problem.createScalarField("rAU");   // 创建 cell 中间场
auto& same = problem.existingScalarField("rAU");// 查找已声明对象
~~~

scalarField/vectorField/tensorField 只加载文件中的 cell 场；create*Field 只创建程序拥有的
场；createFace*Field 只创建面场；existing*Field 只查找同一对象。重复 create、把已创建场
重新加载、或声明阶段之后创建都会抛错，第一次 create 的初值不会静默忽略。临时数学量使用
math 返回值，跨时间步状态使用 time::History。

## math 与 equ

math::div(phi) 是给定面通量的显式散度；equ::div(equation, phi) 是向 equation 绑定的
未知量加入隐式散度。math::interpolate 永远只是插值，不偷偷执行 Rhie–Chow 或压力修正。

~~~cpp
auto equation = equ::createEquation(T);
equation.reset();
equ::ddt(equation, rho * cp, history);
equ::laplacian(equation, k, -1);
equ::source(equation, Q);
const auto solved = equ::solve(equation, linear);
~~~

laplacian(equation, gamma, multiplier) 的定义是 multiplier * div(gamma * grad(unknown))；
扩散左端通常使用 -1。source 只加入已知右端。Equation 的 reset、diagonal、rhs、
referenceIfUnanchored 是对象状态操作；需要 SIMPLE 响应系数时显式写
`geometry::cellVolumes(mesh) / equation.diagonal()`；faceFlux 从实际离散项返回一致面通量。

## geometry 几何场

```cpp
const auto V = geometry::cellVolumes(problem.mesh());
const auto Sf = geometry::faceAreaVectors(problem.mesh());
const auto Af = geometry::faceAreas(problem.mesh());
```

这些值场只反映网格几何，不读取配置、不写结果、不实现物理算法。`Sf` 的方向与面通量
和散度算子一致，因而可以和 `math::dot`、`math::interpolate` 直接组合。

## 时间、配置和诊断

~~~cpp
auto time = time::start(problem);
auto history = time::history(T);
while (time.value() < time.end()) {
    time.advance();
    history.save(T, time.dt());
}
~~~

methods 在 Case 构造时只读取一次，problem.methods() 只读查询；BDF2 首步由 ddt 用 Euler
启动。physics 放物性和模型，solution 放线性/迭代控制，control 放时间，output 放写出筛选。
Parameters 提供 word、boolean、number、positive、nonnegative、fraction、integer、inspect；
Physics 不访问 entry().tokens。readLinearControl(problem, field) 支持按未知场的覆盖。

diagnostics 只观察；Reporter 只格式化、按周期和进程选择报告。收敛由具体 Solver 判断并
返回 SolverResult，application 才映射退出码。
