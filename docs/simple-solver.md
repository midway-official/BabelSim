# 层流不可压缩 SIMPLE：算法驱动入口

BabelSim 的选择名称为 `simple`。它解决稳态、常密度牛顿流体的不可压缩方程，
并未实现湍流模型或瞬态压力速度耦合。核心用户入口是 `src/physics/simple/main.cpp`：

```cpp
int runSimple(Case& problem) {
    SimpleSolver simple(problem);
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

调用者不列出 pPrime、梯度、rAU、面系数、历史速度或工作区，也不操作 RunTime。
SIMPLE 算法本身仍是有状态的；复杂性没有被删掉，而是放在需要理解它的实现边界。
`main.cpp` 包含 `babelsim/application.h`，只注册自己的 `simple` 名称。
通用启动器读取 Case 后按该名称调用 `runSimple`，不需要修改集中名单。

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
| src/physics/simple/state.h | 算法维护者；完整私有状态，只供本算法包含 |

此前 simple_case.cpp 混在一起的读入、外循环、日志和输出已经拆清：通用 Case 读入/输出，
main 组织循环，convergence 输出数值诊断。没有为五个步骤新增五个 Manager 或接口基类。

这不是所有新算法必须采用的目录模板。一个短耦合算法完全可以只写一个函数；
见 [双场测试例子](../tests/examples/coupled_scalar.cpp)。

七个算法文件都只依赖公开头文件与本目录的 `state.h`。没有单独的 SIMPLE 离散文件
放在通用框架中，也没有通过 `internal/` 调用数值核的例外。文件边界依据职责而不是
代码行数：动量插值只有一段短数学组合，保留在 `momentum.cpp`，无需另建类或文件。

## 四种状态分别归谁

| 状态 | 内容 | 所有者与边界 |
| --- | --- | --- |
| 物理状态 | U、p、phi | Case 或嵌入式调用者拥有；算法借用，必须活得比算法久 |
| SIMPLE 算法状态 | pPrime、rAU、phiHbyA、上一轮 U、步骤/迭代/收敛状态 | 本模块 `state.h`；跨算法步骤保持，不向通用 FVM 暴露 |
| 数学中间量 | gradP、rAUGradP、面插值、预测通量散度 | 本模块预分配的 Field；由公开 math/Field 运算更新，不含存储/通信信息 |
| 数值执行工作区 | 梯度重构临时量、通量暂存、矩阵、求解器工作向量、同步缓冲 | 通用 FVM/代数/并行模块独占，SIMPLE 不包含或借用这些内部类型 |

`NumericalWorkspace` 只是本算法几个可复用数学 Field 的分组，不是跨层共享的执行工作区。
没有把整套 SIMPLE 状态搬到框架、再把短入口当成完整算法。维护者仍可从本目录读到全部
压力速度耦合逻辑；通用框架不知道 pPrime、Rhie–Chow 或 SIMPLE 的停止策略。

## 动量与压力的数学边界

动量代码仍直接表达：

```cpp
solveWithResponse(
    eqn::div(rho, phi, U) ==
        -math::grad(p) + eqn::laplacian(mu, U),
    rAU, relaxed(relaxation));
```

rho 是密度，mu 是动力黏度。rAU 是带欠松弛的动量对角对应的体积/对角系数，
供同位压力速度耦合使用；其提取不让 Solver 接触离散矩阵。

预测通量 phiHbyA 由动量预测速度、rAU∇p 与 Rhie–Chow 风格的压力法向重构形成：

\[
\phi^{pred}_f = \operatorname{flux}(U)_f
 + S_f\cdot\operatorname{interpolate}(rAU\nabla p)_f
 - (rAU)_f (S_f\cdot\nabla p)_{\mathrm{discrete}}.
\]

后两个修正项仅作用于内部面。`momentum.cpp` 通过通用数学操作组合它们：

```cpp
math::add(math::flux(work.rAU_grad_p_face), phiHbyA, math::FaceRegion::Interior);
math::subtract(math::flux(work.rAU_face, math::reconstruct(p, work.grad_p)),
              phiHbyA, math::FaceRegion::Interior);
```

`flux(faceVector)` 直接计算面积向量点积；`flux(k,reconstruct(p,gradP))` 复用已算梯度，
使用 p 的扩散格式计算面积积分法向梯度。rAU 仍先按自身插值配置求面值，避免重构时
无意改变方法选择。`Interior` 是几何内部面，包含分区交界，不包括物理边界。
面遍历、几何索引、输入/输出同步属于通用 FVM，而不是上述算法代码。

这一操作组合属于 SIMPLE，不将其提升成核心 `MomentumInterpolation` API。
未来其他耦合算法可以复用这些数学操作；只有确实共享领域算法时才考虑提取领域模块。

压力修正满足

\[
\nabla\cdot(rAU\nabla p')=\nabla\cdot\phi_{H/A}.
\]

代码采用等价的正对角符号：

```cpp
solve(
    -eqn::laplacian(rAU, pPrime) == -eqn::source(divPhiHbyA),
    fix_reference ? referenceValue(0.0) : EquationControl{});
```

随后欠松弛更新 p、用 math 修正 U 和 phi。压力参考、非正交显式项与通量修正保持原有数值路径；
本轮没有以换一种 SIMPLE 公式来换取较短代码。

动量响应和参考规范是通用公开求解能力，不再由 SIMPLE 专用转发函数提供。
矢量 SolveResult 在后端合并三个分量；算法检查一次方程级结果即可。
SimpleSolver 构造不关闭 Case 的声明阶段，组合算法可继续声明额外场；主程序最后调用 start()。

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
日志通过 `diagnostics::report(message)` 输出一次，算法不查询主进程。
`numericalMethods()` 只返回只读方法配置，供稳态检查和非正交修正次数选择，
不返回执行器、矩阵、通信对象或临时工作区。

## 与 OpenFOAM 的对照

对照 [OpenFOAM-8 simpleFoam 主循环](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/applications/solvers/incompressible/simpleFoam/simpleFoam.C)
及 [pEqn.H](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/applications/solvers/incompressible/simpleFoam/pEqn.H)：

| OpenFOAM 思路 | BabelSim |
| --- | --- |
| createFields | Case 命名场与 create_fields.cpp |
| SIMPLE loop | main.cpp 中五个有数值语义的步骤 |
| UEqn | momentum.cpp 中 eqn/math 方程 |
| rAU/HbyA/phiHbyA | 私有算法场、同位动量插值 |
| 非正交压力循环 | pressure.cpp 中压力修正循环 |
| 压力/速度/通量更新 | 明确的 correctVelocity / correctFlux |
| 写出由框架提供 | Case::finish，独立后处理 |

OpenFOAM 的压力步骤并不只是一条普通 laplacian：它还需要约束、参考值、预测通量与修正。
BabelSim 因此保留必要的专用数值语义，但不复制 tmp<>、fvMatrix、MRF 或湍流模型体系。

## 维护与兼容

simple_control.h 和带 IncompressibleFields 参数的构造入口保留给内存型测试。
普通 Case Solver 不使用它们；公开主头不包含 runtime.h 或 simple_control.h。
这避免让每个使用 SIMPLE 的人先学习所有内部对象，同时保持现有数值回归直接可用。

嵌入式构造从 `SimpleSolver(run_time, fields, fluid, control)` 改为
`SimpleSolver(fields, fluid, control)`，所有仓库调用者已迁移。调用前仍须建立活动运行域，
如同调用公开 `solve/math`；算法不接收这个对象，也不能独立管理它的生命周期。
传入不同网格的场、错误位置或没有活动运行域会报错。Case 构造入口不变。

SimpleControl 只保存算法控制；线性求解配置唯一保存在 RuntimeControl，旧的
velocity_solver/pressure_solver 和 simpleRunTimeControl 已移除。内存测试直接配置
RuntimeControl.scalar_solver/vector_solver，再创建 RunTime。
`src/internal/simple_discretization.h` 与 `src/discretization/simple_discretization.cpp`
已经删除，不保留私有桥接或另一套 Rhie–Chow 实现。

`make test-external` 把本目录的所有 C++ 文件及 `state.h` 复制到临时外部工程，
仅提供公开头文件、预编译库，以普通 g++ 编译整个算法，再链接独立可执行程序。
该程序运行 1/2/4 进程 Poiseuille 算例，并检查收敛、单份日志与 U/p 一致性。
这比只编译一个调用现成 SimpleSolver 的入口更严格，但仍不等于已验证独立 PIMPLE 实现。

## 本次边界迁移的验收记录

以 `98ba61d` 为重构前基线，同网格、同物性、同格式和停止准则比较，未增加外迭代：

| 测试 | 验收结果 |
| --- | --- |
| 二维腔流，144 单元，1/2/4 进程 | 均为 137 次外迭代；同进程数新旧 U/p 最大分量绝对差小于 3e-16 |
| 重构后二维腔流，1 对 4 进程 | U/p 最大分量绝对差 7.70e-11 / 3.52e-11 |
| 三维腔流，216 单元，2 进程 | 57 次外迭代；同进程数新旧 U/p 最大差小于 2e-16 |
| 三维非正交腔流，125 单元 | 保持 79 次外迭代；连续性残差约 1.86e-8，存在三维展向流动 |
| Poiseuille，4340 单元，2 进程 | 保持 865 次外迭代；新旧 U/p 最大差约 1.49e-14 / 2.25e-14 |
| 公开 math，三维倾斜网格 | 22 项操作 × 3 种扩散格式，1/2/4 进程主动污染待同步数据后仍与串行参考一致；最大误差小于 8e-14 |

上述 U/p 差异按全局单元编号匹配、逐分量取最大绝对值，不是迭代残差，也不是与解析解的误差。
常规单元/求解器测试、MPI 测试、外部 SDK 构建、Case/后处理工作流和架构负向测试全部通过。
新增错误网格、位置、别名、无活动运行域等拒绝测试，防止接口简化后丢失安全校验。

2 进程 Poiseuille 的交错五次墙钟抽查中位数为重构前 35.37 s、重构后 35.01 s，
包含启动及写出。单次波动明显，只说明该本机样本未见明显退化，不能据此声称提速或
证明大规模 MPI 扩展性。梯度和整场工作数组继续复用；公开调用各自负责同步，
其通信次数不等于此前依赖调用者同步的私有融合核。
