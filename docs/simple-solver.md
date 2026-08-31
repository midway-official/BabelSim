# 层流不可压缩 SIMPLE：算法驱动入口

BabelSim 的选择名称为 `simple`。它解决稳态、常密度牛顿流体的不可压缩方程，
并未实现湍流模型或瞬态压力速度耦合。核心用户入口是 `src/physics/simple/main.cpp`：

```cpp
int runSimple(Case& problem) {
    SimpleSolver simple(problem);
    while (simple.loop()) {
        simple.solveMomentum();
        simple.solvePressure();
        simple.correctVelocity();
        simple.correctFlux();
        simple.checkContinuity();
    }
    return simple.converged() ? 0 : 2;
}
```

调用者不列出 pPrime、梯度、rAU、面系数、历史速度或工作区，也不操作 RunTime。
SIMPLE 算法本身仍是有状态的；复杂性没有被删掉，而是放在需要理解它的实现边界。

## 文件按阅读目的划分

| 文件 | 谁应该阅读、为什么 |
| --- | --- |
| simple/main.cpp | 使用/组织 SIMPLE 的作者；完整算法主循环 |
| simple/momentum.cpp | 修改动量物理的人；动量 PDE 与预测通量步骤 |
| simple/pressure.cpp | 修改压力耦合的人；压力方程、非正交循环、速度和通量修正 |
| simple/create_fields.cpp | 维护初始化的人；Case 配置、场引用和固定工作区构造 |
| simple/convergence.cpp | 维护停止语义的人；全局诊断、日志 |
| simple/simple_solver.cpp | 维护算法执行的人；步骤顺序检查和单次迭代入口 |
| include/babelsim/simple.h | 普通调用者；简短算法接口 |
| src/internal/simple_state.h | 框架/算法维护者；完整私有状态 |
| simple_discretization.cpp | 数值核维护者；动量对角、压力参考和插值执行 |

此前 simple_case.cpp 混在一起的读入、外循环、日志和输出已经拆清：通用 Case 读入/输出，
main 组织循环，convergence 输出数值诊断。没有为五个步骤新增五个 Manager 或接口基类。

这不是所有新算法必须采用的目录模板。一个短耦合算法完全可以只写一个函数；
见 [双场测试例子](../tests/examples/coupled_scalar.cpp)。

## 动量与压力的数学边界

动量代码仍直接表达：

```cpp
simple::solveMomentumEquation(
    fvm::div(rho, phi, U) ==
        -fvc::grad(p) + fvm::laplacian(mu, U),
    relaxation, rAU);
```

rho 是密度，mu 是动力黏度。rAU 是带欠松弛的动量对角对应的体积/对角系数，
供同位压力速度耦合使用；其提取不让 Solver 接触离散矩阵。

预测通量 phiHbyA 由动量预测速度、rAU∇p 与 Rhie–Chow 风格的压力法向重构形成。
压力修正满足

\[
\nabla\cdot(rAU\nabla p')=\nabla\cdot\phi_{H/A}.
\]

代码采用等价的正对角符号：

```cpp
simple::solvePressureCorrectionEquation(
    -fvm::laplacian(rAU, pPrime) == -fvm::source(divPhiHbyA),
    fix_reference);
```

随后欠松弛更新 p、用 fvc 修正 U 和 phi。压力参考、非正交显式项与通量修正保持原有数值路径；
本轮没有以换一种 SIMPLE 公式来换取较短代码。

pPrime/rAU/phiHbyA/UPrevious 是算法状态；梯度、面系数和散度是预分配数值工作场。
不额外保存可由已有速度与压力梯度形成的 HbyA，避免无必要的整场存储。

## 循环和收敛

调用顺序必须是动量 → 压力 → 速度 → 通量 → 连续性；错序或漏步会报错。
loop 不偷偷执行额外外迭代；达到全局收敛或最大次数便停止。

三个状态严格分离：

- healthy：数值有限，没有线性数值失败；
- linear_converged：三个动量分量与每次压力修正线性求解全部满足各自容差；
- converged：前两者成立，且全局连续性、速度变化同时满足外迭代条件。

压力修正相对幅值仅是诊断量，不是默认停止准则。出现线性未收敛不能报告外迭代成功。
所有范数和停止依据都是全局量，不在各 rank 上自行决定是否继续。

## 与 OpenFOAM 的对照

对照 [OpenFOAM-8 simpleFoam 主循环](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/applications/solvers/incompressible/simpleFoam/simpleFoam.C)
及 [pEqn.H](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/applications/solvers/incompressible/simpleFoam/pEqn.H)：

| OpenFOAM 思路 | BabelSim |
| --- | --- |
| createFields | Case 命名场与 create_fields.cpp |
| SIMPLE loop | main.cpp 中五个有数值语义的步骤 |
| UEqn | momentum.cpp 中 fvm/fvc 方程 |
| rAU/HbyA/phiHbyA | 私有算法场、同位动量插值 |
| 非正交压力循环 | pressure.cpp 中压力修正循环 |
| 压力/速度/通量更新 | 明确的 correctVelocity / correctFlux |
| 写出由框架提供 | Case::finish，独立后处理 |

OpenFOAM 的压力步骤并不只是一条普通 laplacian：它还需要约束、参考值、预测通量与修正。
BabelSim 因此保留必要的专用数值语义，但不复制 tmp<>、fvMatrix、MRF 或湍流模型体系。

## 维护与兼容

simple_control.h 和带 RunTime/IncompressibleFields 参数的构造入口保留给内存型测试。
普通 Case Solver 不使用它们；公开主头不包含 runtime.h 或 simple_control.h。
这避免让每个使用 SIMPLE 的人先学习所有内部对象，同时保持现有数值回归直接可用。
