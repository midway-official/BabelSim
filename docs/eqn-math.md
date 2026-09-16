# math / equ 数学与离散语义

BabelSim 保留两个名字空间，因为“已计算的数学场”和“加入待求解离散系统的项”是不同
动作。Physics 开发者只需要理解返回类型和绑定未知量，不需要知道后端矩阵布局。

## math：显式场数学

~~~cpp
const auto gradP = math::grad(p);                    // cell vector
const auto gradU = math::grad(U);                    // cell tensor
const auto phi = math::flux(U);                      // face scalar flux
const auto divPhi = math::div(phi);                 // cell scalar divergence
const auto divStress = math::div(stress);            // cell vector divergence
const auto f = math::interpolate(p);                 // face scalar interpolation
const auto diffusionFlux = math::flux(k, T);         // face diffusion flux
const auto lap = math::laplacian(T);                 // explicit laplacian
~~~

math 操作立即求值并返回 Field。math::div 的参数类型决定它是面通量、矢量场或张量场
的散度。math::interpolate 永远只是插值，没有 Rhie–Chow、压力修正或隐藏的 CFD 算法。
场代数、grad、div、dot、cross、trace、transpose、isotropic 等保持纯数学语义；文档只列出
当前实现已经提供的操作。

Field 的 fill、赋值、+=、-=、缩放是目标场修改。它们不负责历史、通信、输出或 Equation
装配。需要跨时间步保存时使用 time::History；需要 Case 命名场时使用
create*Field。

## equ：绑定未知量的有限体积离散

~~~cpp
auto equation = equ::createEquation(T);  // 明确绑定未知量 T
equation.reset();                         // 清空项，保留结构和绑定
equ::ddt(equation, rho * cp, history);
equ::div(equation, phi);                  // 隐式 div(phi*T)
equ::laplacian(equation, k, -1);          // -div(k*grad(T)) 左端
equ::source(equation, Q);                 // 已知 RHS
const auto result = equ::solve(equation, T, linear);
~~~

Equation 不从字段名猜 PDE。每个 equ 调用立即向同一系统加入一项，并冻结该调用所需的
输入。定义：

- ddt(eq, capacity, history)：capacity 乘时间导数；
- div(eq, phi, scale)：scale 乘隐式 div(phi 乘绑定未知量)；
- laplacian(eq, coefficient, multiplier)：multiplier 乘
  div(coefficient 乘 grad(unknown))；
- reaction：隐式反应对角；
- source：已知右端，不自动线性化；
- relax：按当前 BabelSim 约定修改对角和 RHS；
- faceFlux：从该 Equation 的实际离散项得到一致面通量。

laplacian 的符号必须在调用点可见。通常扩散放到左端时直接使用 -1；不要使用没有
数学含义的额外布尔或字符串参数。

## Equation 状态和只读观察

| API | 含义 |
| --- | --- |
| reset | 清空当前贡献，保留未知量和工作区 |
| copy | 复制独立的已组装系统 |
| diagonal | 积分后的 aP |
| volumeScaledInverseDiagonal | V/aP，用于 SIMPLE 响应 |
| rhs | 积分后的 b |
| reference / referenceIfUnanchored | 标量规范约束 |
| apply / residual | 对当前系统做代数应用 |
| diagnostics::relativeResidual | 读取 b-Ax 的归一化残差 |

这些接口不会自动再次组装。SIMPLE 可以在组装后直接读取 relativeResidual，再进行松弛
和求解，因此不需要为了诊断再次调用 momentum assembly。

## 时间离散

~~~cpp
auto time = time::start(problem);
auto history = time::history(T);
time.advance();
history.save(T, time.dt());
equ::ddt(equation, rho * cp, history);
~~~

methods.bs 在 Case 构造期间读取一次；equ::ddt 从 History 和 problem.methods 的 time
设置选择 Euler/BDF2。BDF2 的首步自动用 Euler。History 不在内层迭代中推进。

## SIMPLE 专用数学保持在 Physics

Rhie–Chow 由 SIMPLE main 组合 math::grad、math::interpolate、math::flux、math::add 和
math::subtract。non-orthogonal pressure correction 也是该 main 的显式子循环。
不要给 math::interpolate 增加特殊 overload，也不要把 SIMPLE 流程包装为 simple.solve()。
