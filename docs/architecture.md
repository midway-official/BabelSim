# BabelSim 架构与维护边界

本文面向框架维护者。Solver 作者请读 [DSL 与运行时用户手册](dsl-runtime-manual.md)。
本页只说明当前源码中的层次、允许的依赖和验证方式。

## 1. 设计目标

Physics code should read as executable mathematics and numerical algorithms：

- Case/Field 负责当前网格上的命名数据和生命周期；
- math 负责已经计算的场数学；
- equ 负责绑定未知量的离散 Equation；
- Solver 自己决定循环、修正、收敛和报告内容；
- Runtime 只安装方法、提供时间/并行执行和后端；
- monitor 只打印观测值，不认识任何具体算法。

运行时和后端不替 Physics 猜 PDE，不打印求解进度，不判断物理收敛。应用层将
SolverResult 映射为命令行退出码。

## 2. 当前层次

| 层 | 当前入口 | 职责 |
| --- | --- | --- |
| Case/IO | include/babelsim/case.h、src/io/case.cpp | 读取案例、拥有 Mesh/命名 Field、参数、输出 |
| 参数 | include/babelsim/config.h、src/io/config.cpp | typed getter、默认值校验、消费跟踪 |
| 场 | include/babelsim/field.h | scalar/vector/tensor、cell/face、边界和值运算 |
| 几何 | include/babelsim/geometry.h、src/geometry | 将 Mesh 的体积、面积、中心和法向表示为普通 Field |
| 方法 | include/babelsim/methods.h、src/io/numerics_reader.cpp | 一次性读取空间/时间离散及按场覆盖 |
| math | include/babelsim/math.h、src/discretization/operators.cpp | grad/div/flux/interpolate 等显式场运算 |
| equ | include/babelsim/equ.h、src/discretization/procedural_equation.cpp | Equation 生命周期、逐项装配、面通量和求解；唯一的方程入口 |
| 时间 | include/babelsim/time.h、history.h | TimeStepper、History；不存储算法状态 |
| 诊断/监视 | solver.h、monitor.h | 只读 residual/change/flux；通用 metric 报告 |
| Runtime/FVM | src/runtime、src/internal/fvm_execution.h | 活动运行域、同步、离散工作区、后端调用 |
| 后端 | src/backend、src/algebra、src/parallel | MPI halo、稀疏装配、线性代数和预条件 |
| Physics | src/physics/heat、transport、simple、transient_simple、RANS | 方程、算法循环、修正、收敛 |
| Application | include/babelsim/application.h、src/runtime/application.cpp | 注册分派、MPI 生命周期、错误映射 |

Physics 不能包含 src/internal、MPI、Eigen、CSR、其它 Solver 的私有头或原始 Field 存储。
后端可以替换，但不改变 Physics 的方程和字段语义。

## 3. 所有权和生命周期

Case 拥有网格、从文件加载的命名场和程序创建的命名场；返回引用在 Case 生存期内稳定。
加载、创建、查找使用三个不同 API：

~~~cpp
auto& T = problem.scalarField("T");
auto& phi = problem.createFaceScalarField("phi");
auto& again = problem.existingFaceField("phi");
~~~

Case 构造期间读取方法一次并创建 Runtime。problem.methods() 只读返回生效配置。参数
声明完成后调用 problem.validate()；start/time::start/setTime 会锁定声明阶段，随后不允许
创建新命名场。面标量场使用 `createFaceScalarField`，面向量/面张量分别使用对应的显式类型名称。
Equation 和 History 的析构顺序由 Case/Runtime 保证，Solver 不管理后端。

## 4. DSL 契约

math 运算要求同一 Mesh 和位置，必要同步由运行时完成；它们返回已经计算的 Field。
equ 操作立即修改已经绑定的 Equation，不延迟解析表达式。Equation 不暴露矩阵布局。
方程 API 只有 equ:: 这一层：表达式式方程、它的解释器和方程级控制都已删除，装配、
时间历史和求解时机全部写在 Physics 源码里，不存在第二套生命周期。div 绑定的面通量
同时成为该方程唯一的边界通量上下文，换通量要 reset() 后重新装配。
equ::faceFlux 读取实际离散项，因此 Physics 不应复制通量公式。

equ::laplacian(eq, coefficient, multiplier) 的统一定义是
multiplier * div(coefficient * grad(unknown))；调用点允许直接使用 -1 表示常用扩散左端。
math::interpolate 只有插值语义。Rhie–Chow、SIMPLE pressure correction 和其它专用算法
只能在相应 Physics main.cpp 中组合。逐项语义、全部算子与配置键见
[DSL 与运行时用户手册](dsl-runtime-manual.md) 第 6、7 节；本文档不再重述签名。

## 5. Solver 独立性

Heat、Transport、稳态 SIMPLE、瞬态 SIMPLE 都是独立注册函数。稳态 SIMPLE 使用一个主要
外循环；瞬态 SIMPLE 使用时间循环加自己的内层循环。两者不共享算法实现或模式开关。
non-orthogonal correction 可以作为各自压力步骤内的显式子循环。

RANS 模型 SA、k–omega、k–epsilon 分别拥有输运场、历史和闭合公式。SIMPLE 只通过
rans::load、effectiveViscosity、deviatoricStressRemainder、solveTransport 交互，不知道
具体闭合。

## 6. 配置和报告边界

physics 只描述物理；methods 只描述离散；solution 只描述线性求解、松弛、外迭代和模型
数值控制；control 只描述物理时间；output 只描述结果筛选。Parameters 的 typed getter
统一拒绝缺失、错误类型、非法范围、重复和未消费键。

diagnostics 只观察，monitor::Reporter 只格式化/按周期/按进程输出；收敛布尔值始终由
Solver 计算。运行时无打印逻辑。Case::output 或 output.bs 选择字段，write() 执行实际
写出，不替 Solver 宣告成功。

## 7. 维护流程

### 改动 DSL

先更新公共头和实现，再同步所有四个 Physics Solver、外部示例和测试。保持旧接口迁移表
与当前实现一致；不要为了表面去重复引入 Manager、Provider、Registry 或深层继承。

### 改动离散或后端

保持 math/equ 的公开数学定义不变，使用 operators/procedural_equation/FVM 测试验证
串行和 1/2/4 rank。后端可维护 MPI、稀疏矩阵和缓存，但不得把这些类型泄漏到 Physics。

### 改动 Physics

保留“加载/读取 -> 物性 -> 清空 Equation -> 逐项装配 -> 求解 -> 修正 -> 诊断/报告”
顺序。稳态和瞬态 SIMPLE 分开修改，Rhie–Chow 不下沉到通用层；每个 RANS 模型独立验证。

## 8. 验收命令

~~~bash
make -j4 all          # lib + solve + post
make -j4 test         # 串行单元/集成测试（含 test-architecture）
make test-workflow    # 新 Solver 工作流、时间序列、ParaView 读取
make test-external    # 仓库外用公共 include 构建 Solver（含负向 API 检查）
make test-rans        # RANS 方程与常数
~~~

涉及并行和非正交时运行 `make test-mpi`、`make test-mpi-poiseuille` 及对应回归。
完整验证入口见 [验证与维护检查](validation.md)。架构测试检查 include 闭包和 Physics
越界依赖；外部 Solver 测试确认只用公共 include。验证报告中的历史结果只表示当时提交的
证据，不替代当前构建测试。
