# 层流不可压缩 SIMPLE Solver

`SimpleSolver` 实现稳态、常密度、层流不可压缩 Navier--Stokes 的压力速度耦合。它不是
通用 FVM 的一部分；通用层提供 `grad/div/flux/interpolate/laplacian`，SIMPLE 层只组织
这些运算与必要的 CFD 专用数值步骤。

## 方程与算法

\[
\nabla\cdot U = 0,
\qquad
\nabla\cdot(\rho\phi U) = -\nabla p + \nabla\cdot(\mu\nabla U),
\qquad \phi=U_f\cdot S_f.
\]

一次 BabelSim 外迭代的结构为：

```text
UEqn
  solve div(rho*phi,U) == -grad(p) + laplacian(mu,U)
  得到 U 和 rAU = V/aP
        ↓
动量插值
  用 U、p、rAU 与 grad(p) 构造抗棋盘格的预测面通量 phiHbyA
        ↓
pEqn / 压力修正
  interpolate(rAU) 与 div(phiHbyA) 形成压力修正方程
  非正交网格执行指定次数的显式修正循环
        ↓
correct p, U, phi
        ↓
FluxBalance
  全局连续性残差与速度变化量决定所有 rank 的共同停止时刻
```

这与 OpenFOAM `simpleFoam` 的 `UEqn.H → pEqn.H → SIMPLE loop` 思路相同。详细的
数据结构和矩阵操作被 Runtime 隐藏，但算法顺序没有被抽象掉。

`SimpleSolver::iterate()` 的主逻辑保持为以下五个步骤；私有工作区和具体非正交实现不进入
这里：

```cpp
result.velocity = solveMomentum();
const PressureStep pressure = solvePressure();
correctVelocity();
correctFlux();
checkContinuityAndConvergence();
```

## 代码中的动量方程

`src/physics/incompressible/momentum_interpolation.cpp` 中的动量预测器写为：

```cpp
solve(
    fvm::div(mass_flux, U) ==
        -fvc::grad(p) + fvm::laplacian(mu, U),
    momentum_control);
```

`mass_flux` 是私有工作区中由 Field 的 `assignScaled(rho, phi)` 形成的 \(\rho\phi\)。
`momentum_control` 仅在 SIMPLE 内部传递速度欠松弛和 `rAU` 输出位；Runtime 在组装后的
对角中计算 \(rAU=V/a_P\)，并完成方程求解。普通 Solver 不会看到这个控制对象。

该 Solver 源码不出现 MPI、通信器、rank、halo、矩阵、LDU、CSR、Eigen 或线性求解器。

## 两个保留的专用数值步骤

OpenFOAM 的压力方程代码也不会把所有内容简化为一条通用 Laplacian：它明确计算 `rAU`、
`HbyA`、预测通量和压力修正。BabelSim 只保留两个相同层级、具稳定数值意义的组件：

### 动量插值

同位网格上，简单插值得到的面通量可能产生压力棋盘格。该组件以

\[
\phi_f = U_f\cdot S_f
+ (rAU\nabla p)_f\cdot S_f
- rAU_f\,S_f\cdot\nabla p
\]

构造 Rhie--Chow 风格的预测通量。内部使用 `fvc::flux`、`fvc::interpolate` 与统一的
非正交 `integratedNormalGradient`；其工作场持久复用，不在每次外迭代分配完整 Field。

### 压力修正

它通过

\[
-\nabla\cdot(rAU\nabla p') = -\nabla\cdot\phi^*
\]

调用通用 `fvm::laplacian` 和 `fvc::div`。其专用职责仅为选择压力参考、控制非正交循环，
并施加

\[
p\leftarrow p+\alpha_p p',\qquad
U\leftarrow U-rAU\nabla p',\qquad
\phi\leftarrow\phi-rAU_fS_f\cdot\nabla p'.
\]

它们是 `SimpleSolver` 的私有实现，不向其他 PDE Solver 暴露长参数表，也没有独立的矩阵、
Field 存储或 MPI 路径。

## 使用与配置

应用启动器将 `physics/incompressible.bs` 读为 `FluidProperties`，将
`numerics/simple.bs` 分成：

```text
Methods：interpolation、gradient、convection、diffusion、time
SIMPLE：maxIterations、nonOrthogonalCorrections、速度/压力欠松弛、外层容差
Runtime：velocitySolver、pressureSolver
```

`SimpleSolver` 只接受 `RunTime`、`IncompressibleFields`、`FluidProperties` 与
`SimpleControl`。RunTime 的 MPI 绑定会自动根据进程是否已经初始化 MPI 决定串行或分布式
路径；SIMPLE 本身不接收 `ParallelContext`。

```cpp
RunTime run_time = RunTime::forMesh(
    mesh, simpleRunTimeControl(methods, control));
SimpleSolver simple(run_time, fields, fluid, control);

for (int i = 0; i < control.max_iterations; ++i) {
    const SimpleIterationResult result = simple.iterate();
    if (!result.healthy || result.converged) break;
}
```

应用代码可显示上述循环；`iterate()` 内部固定保留可读的 SIMPLE 步骤，避免把每个步骤拆为
缺乏数学意义的 Manager/Factory。

## 收敛状态的含义

`SimpleIterationResult` 刻意区分三件事：`healthy` 表示所有线性结果有限且没有数值失败；
`linear_converged` 要求三个动量分量和全部非正交压力修正线性求解都满足各自容差；
`converged` 则额外要求全局连续性和相对速度变化同时满足 SIMPLE 外迭代容差。因此不会出现
“外迭代 `converged=true`、但某个内部线性方程未收敛”的情况。`relative_pressure_correction`
仅是 \(\lVert p'\rVert/\max(\lVert p\rVert,\epsilon)\) 的诊断量，不参与停止判据。

## 与 OpenFOAM simpleFoam 的对照

| OpenFOAM 的组织 | BabelSim 的组织 |
| --- | --- |
| `createFields.H` | Case 读取后构造 `IncompressibleFields` |
| `simple.loop()` | 应用外层迭代与 `SimpleIterationResult` |
| `UEqn.H` | `fvm::div == -fvc::grad + fvm::laplacian` |
| `rAU/HbyA/phiHbyA` | SIMPLE 私有工作区与动量插值步骤 |
| `pEqn.H` | `solvePressure()`、`correctVelocity()`、`correctFlux()` |
| `correctNonOrthogonal()` | `nonOrthogonalCorrections` 循环 |
| 全局 continuity check | `diagnostics::fluxBalance` 与 `diagnostics::all` |

OpenFOAM 把 Rhie--Chow 和压力修正组织在 `pEqn.H` 的算法段，而不是完全变成通用
`fvm::laplacian`。BabelSim 同样保留算法语义，但将两个步骤收敛在 `SimpleSolver` 私有实现中，
并让其内部复用统一 FVM/Runtime 基础设施。

参考 [OpenFOAM simpleFoam 源码](https://github.com/OpenFOAM/OpenFOAM-6/blob/master/applications/solvers/incompressible/simpleFoam/simpleFoam.C)
和其压力修正组织示例 [pEqn.H](https://github.com/OpenFOAM/OpenFOAM-7/blob/master/applications/solvers/multiphase/interFoam/pEqn.H)。
