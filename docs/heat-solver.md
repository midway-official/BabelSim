# Heat 与 Transport 求解器

Heat 和 Transport 是两个独立的单文件 Solver。它们不依赖 SIMPLE 或 RANS，也不把时间
推进隐藏在运行时循环中。

## Heat

仓库实现位于 [src/physics/heat/main.cpp](../src/physics/heat/main.cpp)。主流程对应：

$$
\rho c_p\partial_tT=\nabla\cdot(k\nabla T)+Q.
$$

~~~cpp
auto& T = problem.scalarField("T");
const auto& physics = problem.physics();
const double rho = physics.positive("density");
const double cp = physics.positive("heatCapacity");
const double k = physics.nonnegative("conductivity");
const double Q = physics.number("source");

const auto linear = readLinearControl(problem, T);
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
    if (!solved.converged()) return SolverResult{solved.status};
}
~~~

`ddt` 使用 Case 一次性安装的 methods.time；BDF2 首步由 Equation 使用 Euler 历史。
求解器只在需要的时间步调用 write，运行时不打印和判断物理收敛。

## Transport

Transport 在自己的 main 中加载 C、U，创建面通量 phi：

~~~cpp
auto& C = problem.scalarField("C");
auto& U = problem.vectorField("U");
auto& phi = problem.createFaceScalarField("phi");
phi = math::flux(U);

const double storage = problem.physics().positive("storage");
const double D = problem.physics().nonnegative("diffusivity");
const double S = problem.physics().number("source", 0.0);
~~~

每个时间步按 `ddt -> div(phi) -> laplacian -> source -> solve` 组装 C。phi 是程序创建
的面场，不会从初始目录读取同名文件。两种 Solver 的完整可运行代码都可以作为新标量
PDE 的模板，场加载、Equation、History、诊断和输出职责保持一致。
