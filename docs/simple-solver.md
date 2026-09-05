# 稳态与瞬态 SIMPLE 求解器

本文说明两个内置层流不可压缩求解器的代码边界、数学流程和维护方法。

- `solver simple`：稳态 SIMPLE；
- `solver transientSimple`：瞬态 SIMPLE，每个物理时间步内迭代压力速度耦合。

二者都是具体求解器，不属于 BabelSim 公共 SDK。`include/babelsim/` 中没有
`simple.h`、`simple_control.h` 或 `transient_simple.h`。普通用户通过 Case 选择求解器，
普通 Solver 作者通过公共 `Case + Field + eqn + math + solve + diagnostics` 组合自己的算法；
只有维护内置 SIMPLE 时才阅读 `src/physics/simple*` 中的私有文件。

## 1. 为什么 Heat/Transport 只有一个文件，而 SIMPLE 有多个文件

Heat 和 Transport 是方程驱动求解器：读取场和物性后，在时间循环中求一个方程即可，
所以各自只需要一个 `main.cpp`。SIMPLE 是算法驱动求解器，一个外迭代包含动量预测、
压力方程、速度修正、通量修正和连续性检查，并且多个步骤共享 `pPrime`、`rAU`、
`phiHbyA` 等算法量。将它们全部塞入 main 会降低可读性，因此只在求解器目录内做职责拆分。

这不是向 Framework 增加新层。私有算法对象只是把同一算法的跨步骤状态限制在正确目录；
它不会安装为公共头，也不能被其他 Solver 当作稳定接口依赖。

## 2. 当前文件组织

```text
src/physics/
├── heat/main.cpp
├── transport/main.cpp
├── simple_common.h                 # 两个 SIMPLE 共用的私有物性、控制和结果类型
├── simple/
│   ├── main.cpp                    # 稳态入口和注册
│   ├── algorithm.h                 # 模块私有步骤接口
│   ├── algorithm.cpp               # Case 适配、步骤顺序、收敛和日志
│   ├── momentum.cpp                # 稳态动量方程与动量插值
│   ├── pressure.cpp                # 压力方程及速度/通量修正
│   └── state.h                     # 只供本目录使用的算法状态
└── transient_simple/
    ├── main.cpp                    # 物理时间循环、校正循环和注册
    ├── algorithm.h                 # 模块私有步骤接口
    ├── algorithm.cpp               # 时间步边界、步骤顺序、收敛和日志
    ├── momentum.cpp                # 含 ddt 的动量方程与动量插值
    ├── pressure.cpp                # 压力方程及速度/通量修正
    └── state.h                     # 只供本目录使用的算法状态
```

没有 `simple_discretization.cpp` 放在通用离散层，也没有 SIMPLE 专用 Runtime、MPI、
Matrix 或 Case reader。`simple_common.h` 位于 Physics 内且不随公共头发布；它只是消除
两个 SIMPLE 求解器之间完全相同的四个小数据类型，不包含数值执行实现。

## 3. 稳态入口就是 SIMPLE 外迭代

`src/physics/simple/main.cpp` 的核心流程为：

```cpp
int runSimple(Case& problem) {
    SteadySimpleAlgorithm simple(problem);
    problem.start();
    while (simple.loop()) {
        simple.solveMomentum();
        simple.solvePressure();
        simple.correctVelocity();
        simple.correctFlux();
        simple.checkContinuity();
    }
    return simple.converged() ? 0 : 2;
}

const SolverRegistration simple("simple", runSimple);
```

`SteadySimpleAlgorithm` 是本目录私有名称。Case 用户不会构造它；启动器根据
`case.bs` 中的 `solver simple` 调用 `runSimple`。

稳态动量方程直接表达为：

```cpp
solveWithResponse(
    eqn::div(rho, phi, U) ==
        -math::grad(p) + eqn::laplacian(mu, U),
    rAU, relaxed(velocityRelaxation));
```

对应

\[
\nabla\cdot(\rho\,\phi\,\mathbf U)
=-\nabla p+\nabla\cdot(\mu\nabla\mathbf U).
\]

`solveWithResponse` 同时求速度并得到体积/动量对角响应 `rAU`，但不向算法暴露 LDU、
CSR 或矩阵数组。

## 4. 瞬态入口分清物理时间和校正迭代

`src/physics/transient_simple/main.cpp` 的核心流程为：

```cpp
int runTransientSimple(Case& problem) {
    TransientSimpleAlgorithm simple(problem);
    while (problem.loop()) {
        simple.beginTimeStep();
        while (simple.loop()) {
            simple.solveMomentum();
            simple.solvePressure();
            simple.correctVelocity();
            simple.correctFlux();
            simple.checkContinuity();
        }
        if (!simple.converged()) return 2;
    }
    return 0;
}

const SolverRegistration transient_simple("transientSimple", runTransientSimple);
```

外层 `problem.loop()` 推进物理时间和场历史，内层 `simple.loop()` 只做当前时间步的
压力速度校正。`beginTimeStep()` 不复制 Field、不操作 MPI，只重置该时间步的算法状态。

瞬态动量方程复用稳态写法，只增加物理时间导数：

```cpp
solveWithResponse(
    eqn::ddt(rho, U) + eqn::div(rho, phi, U) ==
        -math::grad(p) + eqn::laplacian(mu, U),
    rAU, relaxed(velocityRelaxation));
```

对应

\[
\rho\frac{\partial\mathbf U}{\partial t}
+\nabla\cdot(\rho\,\phi\,\mathbf U)
=-\nabla p+\nabla\cdot(\mu\nabla\mathbf U).
\]

Euler/BDF2 的场历史由已有 FVM 时间离散处理。同一物理时间步内重复求解不会推进旧时间层。
当前数值范围是“时间步内迭代收敛的瞬态 SIMPLE”，不等同于 PISO/PIMPLE，也没有声称
实现额外的专用时间通量修正。

## 5. 压力修正和动量插值属于 Physics 算法

两个 SIMPLE 使用相同的数值结构：

1. 求动量方程，得到预测速度和 `rAU`；
2. 根据速度、压力梯度和 `rAU` 构造同位网格预测面通量 `phiHbyA`；
3. 求压力修正方程；
4. 修正压力、速度和面通量；
5. 计算全局连续性和速度变化。

预测通量使用公开数学积木组合：

```cpp
math::evaluate(math::grad(p), gradP);
math::evaluate(math::flux(U), phiHbyA);
math::evaluate(math::interpolate(rAU), rAUFace);
math::add(..., phiHbyA, math::FaceRegion::Interior);
math::subtract(..., phiHbyA, math::FaceRegion::Interior);
```

压力修正方程为：

```cpp
solve(
    -eqn::laplacian(rAU, pPrime) == -eqn::source(divPhiHbyA),
    hasFixedPressure ? EquationControl{} : referenceValue(0.0));
```

这些步骤具有明确的 SIMPLE 数值语义，因此保留在 `physics/simple*`，没有错误提升为
通用 Operator，也没有放入 Discretization 或 Runtime。

## 6. 状态所有权

| 状态 | 所有者 |
| --- | --- |
| U、p、phi | Case；算法只借用 |
| pPrime、rAU、phiHbyA、前一校正速度 | 对应 SIMPLE 私有 State |
| 梯度、面插值和预测通量散度的可复用 Field | 对应 SIMPLE 私有 State |
| 时间历史 | FVM/RunTime |
| LDU、线性工作向量 | Algebra/FVM |
| owned/ghost、halo buffer、MPI communicator | Parallel/Runtime |

算法对象必须比它借用的 Case/Field 更早析构。所有工作 Field 在构造时分配一次，
校正循环只更新数值，不重复创建整场临时数组。

## 7. 收敛语义

- `healthy`：数值有限，线性过程没有数值失败；
- `linear_converged`：动量分量和全部压力修正线性求解均达到线性容差；
- `converged`：前两项成立，并且全局连续性和速度相对变化达到 SIMPLE 外迭代容差。

这些状态不能互相替代。所有停止依据通过 `diagnostics` 做全局归约，因此所有 rank
执行相同数量的校正，不允许由本地残差分别决定流程。

## 8. 用户如何运行

普通用户不包含任何 SIMPLE 头，只在 Case 中选择：

```text
solver simple
```

或：

```text
solver transientSimple
```

瞬态 Case 的 `methods.bs` 必须选择 `time euler` 或 `time bdf2`，`control.bs` 必须给出
有效的 `startTime/endTime/deltaT`。稳态 Case 使用 `time steady`。两个求解器都读取 U、p、
密度、动力黏度、SIMPLE 控制以及通用 scalar/vector 线性求解配置。

串行和 MPI 使用同一个入口：

```bash
build/babelsim-solve -case <case>
mpirun -np 4 build/babelsim-solve -case <case>
```

## 9. Solver 作者与框架维护者的边界

若开发新的算法驱动 Solver，可以参考本目录的“步骤组织方式”，但不能包含这里的
`algorithm.h`、`state.h` 或 `simple_common.h`。新 Solver 应在自己的目录使用公共
Case/Field/eqn/math/solve API 表达自己的算法。

只有以下情况才修改 Framework：确实缺少可被多个 PDE 复用的数学算子、边界类型、
离散格式或线性后端。增加新的压力速度耦合步骤、松弛策略或停止准则仍属于新 Solver 自己。

## 10. 与 OpenFOAM 思想的关系

BabelSim 学习的是“主程序体现算法、动量和压力方程有清晰边界、Case 提供物理和数值配置”
的思想。没有复制 `UEqn.H/pEqn.H` 文本包含、对象注册表、复杂模板和历史兼容体系。

因此，BabelSim 的 main 仍能直接读出 SIMPLE 步骤，但具体算法类不对外发布，公共概念面
保持为 Field、Boundary、eqn、math、solve、Case 和 diagnostics。
