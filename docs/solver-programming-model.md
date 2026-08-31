# Solver Programming Model：方程驱动与算法驱动

BabelSim 正式采用两种 Solver 组织方式。它们不是两套 Framework，也不通过共同基类强行
统一；区别只在最上层如何组织计算，底层共同复用 `Field → fvm/fvc → Discretization →
Linear Algebra → RunTime → MPI`。

## Equation-driven Solver

单方程或弱耦合问题首先选择方程驱动模式。Solver 的中心是一个数学方程，时间循环只负责
重复求解它：

```cpp
while (run_time.loop()) {
    solve(
        fvm::ddt(rho_cp, T) ==
            fvm::laplacian(k, T) + fvm::source(Q));
}
```

热传导、Poisson、扩散、对流--扩散和组分输运都属于这一类。当前
`solveTransientHeat()` 和 `solveTransientScalarTransport()` 是黄金模板：它们只使用
Field、物性、边界、`fvm`、`solve()` 和时间循环，不创建矩阵，不操作 Field 存储，也不
接触 MPI。

## Algorithm-driven Solver

多个方程与修正步骤协同工作的算法选择算法驱动模式。Algorithm 组织 Equation 和
Correction，而不是建立另一套离散或线性求解基础设施。当前 SIMPLE 的一次外迭代为：

```cpp
result.velocity = solveMomentum();
const PressureEquationResult pressure = solvePressure();
correctVelocity();
correctFlux();
checkContinuityAndConvergence(result, pressure);
```

其中动量和压力仍是普通 Equation-driven 方程：

```cpp
solve(
    fvm::div(rho, phi, U) ==
        -fvc::grad(p) + fvm::laplacian(mu, U));

solve(
    -fvm::laplacian(rAU, pPrime) ==
        -fvm::source(divPhiHbyA));
```

SIMPLE 只额外组织 `phiHbyA` 动量插值、非正交压力修正、速度/通量修正和连续性诊断。
PIMPLE、Projection Method 或多物理场耦合以后也应遵守同一原则。

## SIMPLE 状态边界

| 类别 | 当前对象 | 生命周期与职责 |
| --- | --- | --- |
| 物理状态 | `U`、`p`、`phi` | 属于案例和最终结果 |
| 算法状态 | `pPrime`、`rAU`、`phiHbyA`、`UPrevious` | 跨一个外迭代的方程和修正步骤 |
| 数值工作区 | `gradP`、`rAUGradP`、`rAUGradPFace`、`rAUFace`、`divPhiHbyA` | 仅用于避免重复分配，不具有独立物理含义 |
| 执行状态 | halo、全局归约、矩阵、预条件器、时间历史 | 完全位于 RunTime/Parallel/Algebra 内部 |

这些对象不会作为长参数表传给专用算子。`SimpleSolver` 私有实现直接复用同一套 `fvm/fvc`
和 `solve()`；公开 Solver API 只保留物理场、物性、控制和迭代结果。

## Pressure Equation 的算法语义

当前实现明确遵循：

```text
Momentum Equation → rAU → phiHbyA
                  → div(phiHbyA)
                  → non-orthogonal pressure equation loop
                  → p correction → U correction → phi correction
```

`HbyA` 没有被强制保存为额外完整 Field：当前动量预测速度和 `rAU∇p` 足以形成数学等价的
`phiHbyA`，因此保留一个额外 `HbyA` 场只会增加内存带宽。压力循环使用具名的
`correctNonOrthogonal()` 语义，不让裸 `pass` 编号进入算法表达。

预测通量所需的 `grad(p)`、`rAU∇p`、面插值和 Rhie--Chow 修正在 Discretization 内部完成；
压力后的 \(U\) 与 \(\phi\) 修正由 `fvc::subtract` 表达。这样 `SimpleSolver` 只保留
动量、压力、速度修正、通量修正与连续性这五个算法步骤，不保存 `grad(p')`，也不含 cell、face
或 patch 循环。

## 共享边界与禁止依赖

两种 Solver 都可以使用：

```text
Mesh、Field、Boundary、Material、fvm、fvc、Equation、solve、RunTime::loop、diagnostics
```

两种 Solver 都不得依赖：

```text
MPI、rank、communicator、halo、ghost、LDU/CSR/Eigen、SparseAssembly、Field 原始存储
```

内部方程控制只用于从动量离散提取 `rAU` 和给无固定压力边界的压力方程设置参考值；它们位于
`src/internal`，不是二次开发 API。

## 与 OpenFOAM 的关系

BabelSim 借鉴 OpenFOAM 的 Field-centric、`fvm/fvc`、`UEqn/pEqn` 和
`correctNonOrthogonal()` 组织思想。OpenFOAM 的不可压 SIMPLE 压力步骤同样明确形成
`rAU`、`HbyA`、`phiHbyA` 并在非正交循环中求压力方程。BabelSim 保留这一数学边界，但不
复制 `tmp<>`、`fvMatrix`、对象注册表、复杂继承或运行时 Factory。

参考：

- [OpenFOAM simpleFoam pEqn.H 官方源码说明](https://cpp.openfoam.org/v7/incompressible_2simpleFoam_2pEqn_8H.html)
- [OpenFOAM 压力非正交修正源码示例](https://cpp.openfoam.org/v6/heatTransfer_2buoyantBoussinesqSimpleFoam_2pEqn_8H_source.html)
