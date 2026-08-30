# 层流不可压缩 SIMPLE Solver

`SimpleSolver` 实现稳态、常密度、层流不可压缩 Navier--Stokes 的压力速度耦合。它不是
通用 FVM 的一部分；通用层提供 `grad/div/flux/interpolate/laplacian`，SIMPLE 层只组织
这些运算与两个必要的 CFD 专用步骤。

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
MomentumInterpolation
  用 U、p、rAU 与 grad(p) 构造抗棋盘格的预测面通量 phiHbyA
        ↓
pEqn / PressureCorrection
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

## 代码中的动量方程

`src/physics/incompressible/simple_solver.cpp` 中的动量预测器写为：

```cpp
solve(
    run_time,
    fvm::div(mass_flux, fields.velocity) ==
        -fvc::grad(fields.pressure) +
        fvm::laplacian(fluid.dynamic_viscosity, fields.velocity),
    momentum_control);
```

`mass_flux` 是 Runtime 高层 `scale(rho, phi, mass_flux)` 产生的 \(\rho\phi\)。
`momentum_control` 的唯一算法信息是速度欠松弛和输出 mobility Field；Runtime 在组装后的
对角中计算 \(rAU=V/a_P\)，并完成方程求解。

该 Solver 源码不出现 MPI、通信器、rank、halo、矩阵、LDU、CSR、Eigen 或线性求解器。

## 两个保留的专用算子

OpenFOAM 的压力方程代码也不会把所有内容简化为一条通用 Laplacian：它明确计算 `rAU`、
`HbyA`、预测通量和压力修正。BabelSim 只保留两个相同层级、具稳定数值意义的组件：

### MomentumInterpolation

同位网格上，简单插值得到的面通量可能产生压力棋盘格。该组件以

\[
\phi_f = U_f\cdot S_f
+ (rAU\nabla p)_f\cdot S_f
- rAU_f\,S_f\cdot\nabla p
\]

构造 Rhie--Chow 风格的预测通量。内部使用 `fvc::flux`、`fvc::interpolate` 与统一的
非正交 `integratedNormalGradient`；其工作场持久复用，不在每次外迭代分配完整 Field。

### PressureCorrection

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

没有独立的矩阵、Field 存储或 MPI 路径。

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

## 与 OpenFOAM simpleFoam 的对照

| OpenFOAM 的组织 | BabelSim 的组织 |
| --- | --- |
| `createFields.H` | Case 读取后构造 `IncompressibleFields` |
| `simple.loop()` | 应用外层迭代与 `SimpleIterationResult` |
| `UEqn.H` | `fvm::div == -fvc::grad + fvm::laplacian` |
| `rAU/HbyA/phiHbyA` | `mobility`、`MomentumInterpolation`、`face_flux` |
| `pEqn.H` | `PressureCorrection::solve/apply` |
| `correctNonOrthogonal()` | `nonOrthogonalCorrections` 循环 |
| 全局 continuity check | `RunTime::fluxBalance` 与全局 `RunTime::all` |

OpenFOAM 把 Rhie--Chow 和压力修正组织在 `pEqn.H` 的算法段，而不是完全变成通用
`fvm::laplacian`。BabelSim 同样保留算法语义，但只用两个组件，并让其内部复用统一 FVM/
Runtime 基础设施。

参考 [OpenFOAM simpleFoam 源码](https://github.com/OpenFOAM/OpenFOAM-6/blob/master/applications/solvers/incompressible/simpleFoam/simpleFoam.C)
和其压力修正组织示例 [pEqn.H](https://github.com/OpenFOAM/OpenFOAM-7/blob/master/applications/solvers/multiphase/interFoam/pEqn.H)。
