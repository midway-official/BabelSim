# fvm / fvc：数学表达与离散语义

## 为什么要区分 fvm 和 fvc

同一个连续算子在有限体积程序中有两种不同用途：

1. 作为未知量的隐式项进入线性方程；
2. 使用当前场直接计算为一个已知 Field。

混淆两者会让 Solver 不清楚何时生成矩阵、何时只是计算场。BabelSim 因而采用与
OpenFOAM 相同的**语义分离**，但不复制它复杂的矩阵模板体系。

```text
fvm::operator  轻量方程项 → solve → Runtime 内部离散方程 → 线性系统
fvc::operator  轻量求值描述 → fvc::evaluate → Runtime 内部计算 → 已计算的 Field
```

两类描述都只保存 Field 引用和常数。不会在表达式构造、`+`、`-` 或 `==` 时复制大型场，
不会分配 LDU/CSR，也不会执行 MPI。

## fvm 方程贡献

除显式体源外，以下输运项隐式离散到未知量方程。

| API | 数学意义 | 未知量位置 |
| --- | --- | --- |
| `fvm::ddt(T)` | \(\partial T/\partial t\) | cell scalar |
| `fvm::ddt(rhoCp, T)` | \(\rho c_p\partial T/\partial t\) | cell scalar 系数与场 |
| `fvm::ddt(rho, U)` | \(\rho\partial U/\partial t\) | cell scalar 系数与 vector 场 |
| `fvm::div(phi, T/U)` | \(\nabla\cdot(\phi T/U)\) | face scalar 通量、cell 场 |
| `fvm::div(rho,phi,U)` | \(\nabla\cdot(\rho\phi U)\) | 常数缩放的 face 通量、cell vector 场 |
| `fvm::laplacian(k, T/U)` | \(\nabla\cdot(k\nabla T/U)\) | 常数、cell 或 face 系数 |
| `fvm::source(Q)` | 已知体源 \(Q\) | 常数或 cell scalar Field |
| `fvm::source(a,C)` | 已知体源 \(aC\) | 常数系数与 cell scalar Field；不生成临时场 |

例如瞬态对流扩散方程可写为：

```cpp
solve(
    fvm::ddt(rhoCp, concentration)
      + fvm::div(phi, concentration)
        == fvm::laplacian(diffusivity, concentration)
         + fvm::source(source));
```

这里 `rhoCp` 可预先作为 cell Field 建立；FVM 后端在 `solve` 时验证位置、同步 ghost、
将 cell 扩散系数插值到 face，并根据 `Methods` 选择对流、梯度和非正交扩散处理。

`fvm::laplacian` 在表达式中位于右端时对应正扩散项。等价地可以将其取负放在左端；不支持
把正 Laplacian 直接放左端，因为这通常与正定扩散矩阵的符号约定相反，Runtime 会报出清晰错误。

## 显式 fvc 量

`fvc` 不组装方程。调用者提供结果 Field，因此可以复用工作场，避免隐藏的大临时量：

```cpp
VectorField& grad_p = problem.vectorField("gradP", Vec3{});
ScalarField& phi = problem.faceField("phi");
ScalarField& div_phi = problem.scalarField("divPhi", 0.0);

fvc::evaluate(fvc::grad(p), grad_p);
fvc::evaluate(fvc::flux(U), phi);
fvc::evaluate(fvc::div(phi), div_phi);
```

当前显式描述包括：

| API | 输出 | 含义 |
| --- | --- | --- |
| `fvc::grad(T)` | cell vector | \(\nabla T\) |
| `fvc::grad(U)` | cell tensor | \(\nabla U\) |
| `fvc::flux(U)` | face scalar | \(U_f\cdot S_f\) |
| `fvc::div(phi)` | cell scalar | \(\nabla\cdot\phi\) |
| `fvc::div(U)` | cell scalar | \(\nabla\cdot U\) |
| `fvc::div(phi,T/U)` | cell scalar/vector | 显式对流散度 |
| `fvc::interpolate(cellField)` | face scalar/vector | 面中心重构 |
| `fvc::reconstruct(field,grad)` | face scalar/vector | 已知梯度的偏斜面重构 |
| `fvc::laplacian(k,T)` | cell scalar | 显式扩散散度 |
| `fvc::subtract(rAU, fvc::grad(p'), U)` | cell vector | 原位执行 \(U\leftarrow U-rAU\nabla p'\) |
| `fvc::subtract(fvc::flux(rAU,p'), phi)` | face scalar | 原位执行 \(\phi\leftarrow\phi-rAU_fS_f\cdot\nabla p'\) |

对流显式求值使用 `Methods::convection` 的 Upwind/Central 选择；Central 且选择 corrected
插值时会进行面中心偏斜修正。扩散、梯度和插值的非正交实现由 Runtime 自动选择当前方法，
调用者不应自行 halo 或重构 ghost。

`fvc::subtract` 是面向修正算法的原位显式操作：Runtime 复用梯度、面系数和面通量工作场，
再同步结果。它避免 SIMPLE Solver 为 `grad(p')`、`rAU_f` 和每个 face 的扩散通量编写
cell/face 循环；该 API 不生成隐式方程，也不暴露工作场的存储。

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

所有这些细节属于 Method/FVM/Runtime 层。Solver 只在 Case 的 `numerics` 文件选择方法。

## 时间层与内迭代

同一时间步可以重复求解一个 Field，而不推进旧时间层。RunTime 只在下一次时间推进时更新历史，
因此非线性修正和双场耦合复用同一个 `solve()`。BDF2 在无第二个历史层时以 Euler 启动；
后续仍使用均匀步长公式。当前不支持自适应 BDF2。

普通 Solver 通过 `Case::loop()` 取得这套时间语义和自动输出，不构造 RunTime。
`case.h`、`solver.h` 的公开头文件闭包不包含 MPI、Eigen 或内部运行类。

## 求解时发生的事情

`solve(equation)` 的固定顺序为：

```text
检查表达式与 Field 位置
→ 同步该算子所需的局部/ghost 输入
→ 按 Methods 组装 owned 行
→ 施加 Field 边界条件
→ 选择串行或分布式线性后端
→ 用全局残差判断线性收敛
→ 写回 unknown 并同步 ghost
```

这个顺序解释了为什么 Solver 不应调用 `data()`、`HaloExchange` 或 `Equation::discrete()`：
这些会绕开 Runtime 的一致性与生命周期管理。

## 隐式与显式的选择

一般规则：对未知输运 Field 的稳定主导项使用 `fvm`，用于算法修正、诊断或已知场的量使用
`fvc`。例如 SIMPLE 的动量方程用 `fvm::div`、`fvm::laplacian`，但压力梯度是
`fvc::grad(p)`；压力修正中的 `div(phi)` 与 `interpolate(rAU)` 是显式量。

这也是 OpenFOAM `UEqn.H`/`pEqn.H` 的设计思想。BabelSim 只保留这一清晰分界，而没有复制
OpenFOAM 的 `tmp<>`、`fvMatrix` 或对象注册机制。
