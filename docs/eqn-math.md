# eqn / math：数学表达与离散语义

## 为什么要区分 eqn 和 math

同一个连续算子在有限体积程序中有两种不同用途：

1. 构造待求解方程的贡献，包括未知量的隐式项和已知源项；
2. 使用当前场直接计算为一个已知 Field。

混淆两者会让 Solver 不清楚何时生成矩阵、何时只是计算场。BabelSim 因而采用与
OpenFOAM 相同的**语义分离**，但使用表达用途的名称，不沿用它的命名或复杂矩阵模板体系。

```text
eqn::operator   轻量方程项 → solve → FVM 数值前端 → 离散方程 → 计算后端
math::operator  轻量求值描述 → math::evaluate → FVM 数值前端 → 已计算的 Field
```

两类描述都只保存 Field 引用和常数。不会在表达式构造、`+`、`-` 或 `==` 时复制大型场，
不会分配 LDU/CSR，也不会执行 MPI。

名称不改变求值时机：`math::grad(p)` 是数学量的轻量描述，不会偷偷分配结果 Field。
独立计算时写 `math::evaluate(math::grad(p), gradP)`；用作方程的已知项时可写
`-math::grad(p)`，由 `solve()` 在执行该方程时求值。数学量与方程项的区别是用途，
不是“创建表达式时是否立刻执行”。这种约定保留工作场复用和原有并行同步契约。

## 命名迁移与文件边界

| 旧公开名称/文件 | 当前名称/文件 |
| --- | --- |
| `fvm::`、`babelsim/fvm.h` | `eqn::`、`babelsim/eqn.h` |
| `fvc::`、`babelsim/fvc.h` | `math::`、`babelsim/math.h` |
| `src/discretization/fvm_expression.cpp` | `src/discretization/equation_expression.cpp` |
| `tests/fvc_runtime_test.cpp`、`tests/parallel_fvc_test.cpp` | `tests/math_runtime_test.cpp`、`tests/parallel_math_test.cpp` |
| `docs/fvm-fvc.md` | `docs/eqn-math.md` |

不保留旧命名空间别名或转发头。外部 Solver 同步修改 include 和限定名后重新编译；
原 `solver.h` 聚合入口不变，Case 文件无需更改。
方程描述的 `FvmTermKind/ScalarFvmTerm/VectorFvmTerm` 改为
`EquationTermKind/ScalarEquationTerm/VectorEquationTerm`，仍是同样的轻量项，不是矩阵。

内部 `FvmExecution` 与 `fvm_execution.h/.cpp` 保留：FVM 在这里明确表示有限体积数值前端，
不是公开操作的名字。它负责把轻量数学描述解释为已有算子和 `DiscreteEquation`，但不再
包含 Eigen、稀疏装配、Halo 或 MPI。`ComputeBackend` 接收离散方程，负责同步、归约、
装配与线性求解。本次不宣称实现了 FEM/FDM 或动态后端插件。

## eqn 方程贡献

除显式体源外，以下输运项隐式离散到未知量方程。

| API | 数学意义 | 未知量位置 |
| --- | --- | --- |
| `eqn::ddt(T)` | \(\partial T/\partial t\) | cell scalar |
| `eqn::ddt(rhoCp, T)` | \(\rho c_p\partial T/\partial t\) | cell scalar 系数与场 |
| `eqn::ddt(rho, U)` | \(\rho\partial U/\partial t\) | cell scalar 系数与 vector 场 |
| `eqn::div(phi, T/U)` | \(\nabla\cdot(\phi T/U)\) | face scalar 通量、cell 场 |
| `eqn::div(rho,phi,U)` | \(\nabla\cdot(\rho\phi U)\) | 常数缩放的 face 通量、cell vector 场 |
| `eqn::laplacian(k, T/U)` | \(\nabla\cdot(k\nabla T/U)\) | 常数、cell 或 face 系数 |
| `eqn::Sp(a,T/U)` | 隐式线性项 \(aT/aU\)，按有符号系数加入对角 \(aV\) | 常数或 cell scalar 系数，cell scalar/vector 未知量 |
| `eqn::source(Q)` | 已知体源 \(Q\) | scalar/Vec3 常数或 cell scalar/vector Field |
| `eqn::source(a,F)` | 已知体源 \(aF\) | 常数系数与 cell scalar/vector Field；不生成临时场 |

例如瞬态对流扩散方程可写为：

```cpp
solve(
    eqn::ddt(rhoCp, concentration)
      + eqn::div(phi, concentration)
        == eqn::laplacian(diffusivity, concentration)
         + eqn::source(source));
```

这里 `rhoCp` 可预先作为 cell Field 建立；FVM 数值前端在 `solve` 时验证位置，经计算后端
同步 ghost，将 cell 扩散系数插值到 face，并根据 `Methods` 选择对流、梯度和非正交扩散处理。

`eqn::laplacian` 在表达式中位于右端时对应正扩散项。等价地可以将其取负放在左端；不支持
把正 Laplacian 直接放左端，因为这通常与正定扩散矩阵的符号约定相反，FVM 数值前端会报出清晰错误。

## 显式 math 量

`math` 不组装方程。调用者提供结果 Field，因此可以复用工作场，避免隐藏的大临时量：

```cpp
VectorField& grad_p = problem.vectorField("gradP", Vec3{});
ScalarField& phi = problem.faceField("phi");
ScalarField& div_phi = problem.scalarField("divPhi", 0.0);

math::evaluate(math::grad(p), grad_p);
math::evaluate(math::flux(U), phi);
math::evaluate(math::div(phi), div_phi);
```

当前显式描述包括：

| API | 输出 | 含义 |
| --- | --- | --- |
| `math::grad(T)` | cell vector | \(\nabla T\) |
| `math::grad(U)` | cell tensor | \(\nabla U\) |
| `math::normalGradient(T)` | face scalar | 当前扩散格式下的法向梯度（不是面积积分值） |
| `math::flux(k,T)` | face scalar | 扩散面通量；k 可为 cell 或 face scalar |
| `math::flux(k,math::reconstruct(T,gradT))` | face scalar | 同上，但复用给定的 cell 梯度，不重复重构 |
| `math::flux(U)` | face scalar | \(U_f\cdot S_f\)；cell 输入先插值，face 输入直接点积 |
| `math::div(phi)` | cell scalar | \(\nabla\cdot\phi\) |
| `math::div(U)` | cell scalar | \(\nabla\cdot U\) |
| `math::div(phi,T/U)` | cell scalar/vector | 显式对流散度 |
| `math::interpolate(cellField)` | face scalar/vector | 面中心重构 |
| `math::reconstruct(field,grad)` | face scalar/vector | 已知梯度的偏斜面重构 |
| `math::laplacian(k,T)` | cell scalar | 显式扩散散度 |
| `math::subtract(rAU, math::grad(p'), U)` | cell vector | 原位执行 \(U\leftarrow U-rAU\nabla p'\) |
| `math::subtract(math::flux(rAU,p'), phi)` | face scalar | 原位执行 \(\phi\leftarrow\phi-rAU_fS_f\cdot\nabla p'\) |
| `math::add(math::flux(faceVector),phi,region)` | face scalar | 按指定面区域累加面积积分通量 |

对流显式求值使用 `Methods::convection` 的 Upwind/LinearUpwind/Central 选择。
LinearUpwind 从迎风单元按梯度重构到面中心；隐式方程采用一阶迎风系数作为稳定基底，
二阶重构量以延迟修正加入右端，因此不会复制完整离散矩阵。Central 且选择 corrected
插值时会进行面中心偏斜修正。扩散、梯度和插值的非正交实现由 FVM 数值前端自动选择当前方法，
调用者不应自行 halo 或重构 ghost。

`math::subtract` 是面向修正算法的原位显式操作：FVM 数值前端复用梯度、面系数和面通量工作场，
再同步结果。它避免 SIMPLE Solver 为 `grad(p')`、`rAU_f` 和每个 face 的扩散通量编写
cell/face 循环；该 API 不生成隐式方程，也不暴露工作场的存储。

面通量的 `add/subtract` 可选 `math::FaceRegion::All`（默认）或 `Interior`。
后者只更新几何内部面，包含并行分区交界，物理边界保持原值。面区域不是 rank、patch
或 ghost 选择器。SIMPLE 用它组合 Rhie–Chow，但通用 FVM 不含 Rhie–Chow 公式。
给定重构梯度的输入必须对应当前场；调用者若改变 T，应先更新 gradT。框架负责同步
这两个输入，却不会把用户给定的数学量偷偷改成另一种重构。可传 face 系数以精确保留
算法选择的插值结果；传 cell 系数时仍按被微分场 T 的格式进行插值。
增量复用执行层面通量工作场，不分配新的整场数组，不产生临时矩阵；多个公开调用
分别履行同步契约，不承诺与融合单面核具有完全相同的通信次数。

## 公开同步契约

只包含 math.h 即可调用上述 evaluate/add/subtract。
所有参与进程必须按同样顺序调用；框架同步每个输入，计算，再同步输出。
调用后 cell 结果的 ghost 和分区共享面值已可被后续算子使用。返回并不表示每个 rank 持有全局 Field。

旧公开单面局部核 `integratedNormalGradient(T,grad,face)` 已被移除，因为逐面局部核无法
独立保证整场通信顺序。需要法向导数时声明一个 faceField，再调用
`math::evaluate(math::normalGradient(T), normal)`。
面积积分通量应使用 `math::evaluate(math::flux(k,T), flux)`，不要遗漏面积因子。

evaluate 不支持同位输入和结果使用同一 Field；例如 laplacian(T,T)、div(phi,U,U)
会拒绝。修正算法使用明确的 subtract 原位接口。
面扩散通量的系数不能与目标通量别名；错误网格、位置和未知面区域也会拒绝。
Field::assign/fill/evaluate 等点运算本身不通信，下一个 eqn/math 入口负责所需同步。

## 方程控制和响应

`solve(eq, relaxed(alpha))` 使用方程级欠松弛，alpha 在 (0,1]。
标量形式增加对角及对应旧值源项；矢量形式沿用 SIMPLE 的整体行缩放表示。
`solveWithResponse(eq,rAU,relaxed(alpha))` 额外输出当前缩放表示下的 V/aP，
没有返回矩阵或指针。注意它不是对原始未缩放物理源的导数：
矢量欠松弛已经把原始源乘 alpha；若按原始物理源解释增量，需计入该系数。
该约定保持现有 SIMPLE 的数学路径，不声称与其他软件的 rAU 缩放约定逐项相同。

`solve(eq, referenceValue(c))` 给具有常数零空间、满足相容条件的标量方程设参考值，
锚定固定全局单元，与进程数无关。它不是一般区域约束或边界条件替代；
矢量方程不能使用标量规范。公开 solve 返回统一 SolveResult，
三个矢量分量必须全部收敛，健康失败不能被合并成成功。

## 非正交与偏斜处理

对内部面，离散法将面积积分梯度分成紧凑的正交隐式部分和显式非正交部分：

\[
S_f\cdot\nabla\phi
= a_f(\phi_N-\phi_P)
+ \left(S_f-a_fd_{PN}\right)\cdot(\nabla\phi)_f.
\]

`DiffusionMethod::Corrected` 使用完整显式修正；`LimitedCorrected` 限制非正交修正；
`Orthogonal` 只使用第一项。真实面中心不在 cell-centre 连线上的偏斜由 corrected
interpolation/Green--Gauss reconstruction 处理。梯度可选择 Green--Gauss 或 Least-Squares。

所有这些细节属于 Method/FVM 层。Solver 只在 Case 的 `numerics` 文件选择方法。

## 时间层与内迭代

同一时间步可以重复求解一个 Field，而不推进旧时间层。RunTime 在下一次时间推进时通知 FVM 执行层更新历史，
因此非线性修正和双场耦合复用同一个 `solve()`。BDF2 在无第二个历史层时以 Euler 启动；
后续仍使用均匀步长公式。当前不支持自适应 BDF2。

普通 Solver 通过 `Case::loop()` 取得这套时间语义和自动输出，不构造 RunTime。
`case.h`、`solver.h` 的公开头文件闭包不包含 MPI、Eigen 或内部运行类。

## 求解时发生的事情

`solve(equation)` 的固定顺序为：

```text
FVM 数值前端检查表达式与 Field 位置
→ 经 ComputeBackend 同步该算子所需的局部/ghost 输入
→ FVM 按 Methods 形成 owned 行的 DiscreteEquation 并施加边界条件
→ ComputeBackend 装配其线性代数表示
→ 默认后端选择串行或分布式线性求解并用全局残差判断收敛
→ ComputeBackend 写回 unknown 并同步 ghost
```

这个顺序解释了为什么 Solver 不应调用 `data()`、`HaloExchange` 或离散方程入口：
这些会绕开 Runtime 的一致性与生命周期管理。

接口调用以整场同步、一次全局归约或一套离散方程为单位，不在 cell/face 热循环内进行虚调用。
因此默认 Eigen/MPI 实现保持原有连续存储、几何缓存、稀疏模式复用和通信路径。维护者可以在
构建时替换计算后端，而 Solver 与上述 FVM 数值代码无需修改。

## 隐式与显式的选择

一般规则：对未知输运 Field 的稳定主导项使用 `eqn`，用于算法修正、诊断或已知场的量使用
`math`。例如 SIMPLE 的动量方程用 `eqn::div`、`eqn::laplacian`，但压力梯度是
`math::grad(p)`；压力修正中的 `div(phi)` 与 `interpolate(rAU)` 是显式量。

这也是 OpenFOAM `UEqn.H`/`pEqn.H` 的设计思想。BabelSim 只保留这一清晰分界，而没有复制
OpenFOAM 的 `tmp<>`、`fvMatrix` 或对象注册机制。


## 场的边界约束与派生边界迹

默认 cell Field 表示可求解未知量，值更新不修改其 patch 约束。历史场的 `assign`
只复制当前数值，主未知量的 `fill/evaluate` 也不会覆盖固定边界值。

已知系数、点值函数和其他派生量可在声明时调用 `useCalculatedBoundary()`：

```cpp
auto& coefficient = problem.scalarField("coefficient", 0.0);
coefficient.useCalculatedBoundary();
coefficient.evaluate(T, [](double t) { return 1.0 + t*t; });
// 单元值与每个物理边界面的迹都按同一个函数更新。
solve(eqn::ddt(T) + eqn::Sp(2.0, T) == eqn::laplacian(coefficient, T));
```

`fill/evaluate/assign/assignScaled/addScaled/assignProduct/addProduct` 都遵守这个显式模式，
支持原位点值代数；例如 g=f² 的边界迹为 (trace f)²，而不是自动继承 f 的边界类型。
模式不依赖字段名称、单位、物性或方程类别。派生迹和通量上下文均为值快照，不持有源场
或 lambda 的长期引用；输入改变后，应重新执行对应派生计算。

calculated Field 不能作为 `solve` 的未知量，也不能用 `boundary/setBoundary` 访问或
设置 patch 约束。它的每面迹可通过 `math::interpolate` 使用；无需向 Solver 暴露 face 索引。
`math::grad` 的 calculated 输出由输入固定值、法向梯度或对称约束修正边界梯度。
散度、对流散度和 Laplacian 的 calculated 输出采用相邻单元常值外推，因为仅给输入值
约束不能唯一确定其高阶导数边界。该外推是明确的零阶边界近似，不承诺边界导数的高阶精度。
执行层先同步结果单元值，再更新边界迹，保证分区上复合算子的输入一致。

`inletOutlet` 是依赖面通量符号的混合约束：外向通量小于 0 时 fixedValue，否则 zeroGradient。
方程中的对流项在装配任何项之前提供上下文，使梯度、扩散、对流共享同一面条件；同一
方程不接受相互不同的对流通量对象。独立数学运算可调用 `field.setBoundaryFlux(phi)`；
没有通量上下文时明确报错。`math::flux(U)` 按目标中的前次通量判定本次边界，然后更新
U 的上下文，因此迭代调用应复用目标；初次目标通常为零通量。

## 通用张量与隐式线性项

Tensor3 使用行 i 列 j 表示分量。`transpose/trace` 和加减、标量乘法不带物理语义。
`math::evaluate(math::div(TensorField), vector)` 计算
\((\operatorname{div}T)_i=\sum_j\partial_j T_{ij}\)，每个面使用 \(T_f S_f\)，
按行复用现有插值与梯度重构格式；内部面对相邻单元贡献等大反向。

`eqn::Sp(a,F)` 表示隐式乘法算子 aF，放在左侧时增加 aV 对角，移到右侧则改变符号。
它不自动判定生产或破坏，不自动裁剪负系数。怎样线性化非线性反应项由 Solver 决定；
对于正的衰减率 a，`Sp(a,F)` 在左侧通常有利于隐式稳定性。

## 残差和成功状态

线性求解器和 MPI 后端共用以下定义（全局 owned 行、未预条件 L2 范数）：

\[
s=\max(\|b-Ax_0\|_2,\|b\|_2),\quad
\|b-Ax\|_2\le\max(\mathrm{atol},\mathrm{rtol}\,s).
\]

只有 s=0 时归一化分母使用 1，此时初始解已满足零残差。最终状态重新计算真残差，
不因舍入停滞、迭代递推残差较小或绝对内积较小而放宽用户容差。
有限性失败不能标为收敛；剩余预算允许时串行/MPI 都可从真残差重启 Krylov。
当前双精度平方和范数仍有极端上溢/下溢范围限制，不等于任意指数范围求解保证。

`diagnostics::residual(equation)` 重装配**原方程**，不施加欠松弛或参考点惩罚，不求解，
不推进历史。它返回 `EquationResidual{norm,scale}`，其中
\(norm=\|Ax-b\|_2\)，\(scale=\|Ax\|_2+\|b\|_2\)。
`relative()` 在非零尺度下为二者之比；零尺度返回 norm，非有限输入返回无穷大。
`converged(atol,rtol)` 可组合绝对/相对阈值。该尺度用于当前非线性方程，不等同于
线性求解的初始残差尺度；二者各有明确含义，均不设置隐藏容差下限。

Solver 必须在更新系数、修正未知量或裁剪之后再检查原方程。`SolveResult::Converged`
只说明刚求解的线性系统满足容差，不能证明之后修改的场仍满足原方程。
