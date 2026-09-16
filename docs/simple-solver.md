# SIMPLE 求解器

稳态 SIMPLE 和瞬态 SIMPLE 是两个独立的 Solver，分别位于
[src/physics/simple/main.cpp](../src/physics/simple/main.cpp) 和
[src/physics/transient_simple/main.cpp](../src/physics/transient_simple/main.cpp)。两者
不共享算法实现或模式开关；每个 main 都直接显示自己的数值流程。

## 稳态结构

1. 加载 U、p，创建面通量 phi；
2. 读取 physics 物性、solution 的迭代/松弛/容差和 problem.methods；
3. 读取按未知场的线性控制；
4. 加载可选 RANS 模型，创建 pPrime 和两个 Equation；
5. 检查 steady time method；
6. 进入一个 SIMPLE 外循环。

每次外循环都是：

~~~cpp
const VectorField U_previous = U;
momentumEquation.reset();
equ::div(momentumEquation, phi, rho);
equ::laplacian(momentumEquation, muEff, -1);
equ::source(momentumEquation, -math::grad(p));

const double rU = diagnostics::relativeResidual(momentumEquation, U);
// relax -> volumeScaledInverseDiagonal -> velocity solve

// 本 main 中显式组合 Rhie–Chow
const auto gradP = math::grad(p);
auto phiH = math::flux(U);
math::add(math::flux(math::interpolate(rAU * gradP)), phiH, math::FaceRegion::Interior);
math::subtract(math::flux(math::interpolate(rAU), p, gradP),
               phiH, math::FaceRegion::Interior);

pressureCorrectionEquation.reset();
equ::laplacian(pressureCorrectionEquation, rAU, -1);
equ::source(pressureCorrectionEquation, -math::div(phiH));
pressureCorrectionEquation.referenceIfUnanchored(0.0);
// pressure solve -> p/U/phi correction -> turbulence -> convergence
~~~

残差在组装后、松弛和求解前读取，避免再次组装动量方程。压力 Equation 可以有
non-orthogonal 修正子循环；它仍属于这个 main 的数值步骤。

## 瞬态结构

瞬态版本先创建 `time::TimeStepper` 和 U 的 `History`，再使用：

~~~cpp
while (time.value() < time.end()) {
    time.advance();
    UHistory.save(U, time.dt());

    for (int iter = 0; iter < maxIterations; ++iter) {
        // ddt + momentum assembly
        // residual before relaxation
        // velocity solve
        // this main's Rhie–Chow
        // pressure correction and optional non-orthogonal sub-loop
        // turbulence, diagnostics, convergence
        if (converged) break;
    }
    if (!converged) return SolverResult::notConverged();
}
~~~

瞬态每个物理步只保存一次历史，内层迭代不推进时间。瞬态特有的压力/速度耦合松弛保持
在自己的 main 中。

## RANS 交互

SIMPLE 不知道 k–epsilon、k–omega 或 SA 的闭合：

~~~cpp
auto turbulence = rans::load(problem, U, phi);
const auto& muEff = turbulence.effectiveViscosity();
equ::source(momentumEquation,
            math::div(turbulence.deviatoricStressRemainder(U)));
const auto transport = turbulence.solveTransport();
~~~

每个模型自己的输运 Equation、历史、边界和闭合常数位于 RANS 子目录。模型数值控制
位于 solution，physics 只保留物理常数和模型选择。Rhie–Chow 不能移入 RANS、math 或
equ 算子。
