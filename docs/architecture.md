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
| 方法 | include/babelsim/methods.h、src/io/numerics_reader.cpp | 一次性读取时间和具名方程/项/算子离散配置 |
| math | include/babelsim/math.h、src/discretization/operators.cpp | grad/div/flux/interpolate 等显式场运算 |
| equ | include/babelsim/equ.h、src/discretization/procedural_equation.cpp | Equation 生命周期、逐项装配、面通量和求解；唯一的方程入口 |
| 时间 | include/babelsim/time.h、history.h | TimeStepper、History；不存储算法状态 |
| 诊断/监视 | solver.h、monitor.h | 只读 residual/change/flux；通用 metric 报告 |
| Runtime/FVM | src/runtime、src/internal/fvm_execution.h | 活动运行域、同步、离散工作区、后端调用 |
| 后端 | src/backend、src/internal/petsc_*、src/parallel | PETSc 稀疏系统生命周期，场同步、MPI halo 与并行结果写出 |
| Physics | src/physics/heat、transport、simple、transient_simple、piso、RANS | 方程、算法循环、修正、收敛 |
| Application | include/babelsim/application.h、src/runtime/application.cpp | 注册分派、MPI 生命周期、错误映射 |

Physics 不能包含 src/internal、MPI/PETSc、Eigen、CSR、其它 Solver 的私有头或原始 Field 存储。
后端可以替换，但不改变 Physics 的方程和字段语义。

默认数值后端使用 PETSc 的 Mat/Vec/KSP/PC。每个具名方程拥有长期存在的后端线性系统；首次装配
建立稀疏模式，之后按需要更新矩阵数值与右端，并可仅更新 RHS 后重复求解。PETSc 类型和句柄留在
backend/internal 层，不进入 Physics 或公开 Equation 抽象。线性系统选择从 `solution.bs` 的
`equation.<name>.kspType` 和 `pcType` 读取；压力方程可选 Hypre BoomerAMG 或 PETSc GAMG。
Eigen 仍用于离散算子中的局部最小二乘，不再是稀疏线性求解后端。

## 3. 所有权和生命周期

Case 拥有网格、从文件加载的命名场和程序创建的命名场；返回引用在 Case 生存期内稳定。
加载、创建、查找使用三个不同 API：

~~~cpp
auto& T = problem.scalarField("T");
auto& phi = problem.createFaceScalarField("phi");
auto& again = problem.existingFaceField("phi");
~~~

Case 构造期间读取方法一次并创建 Runtime。problem.methods() 只读返回生效配置。参数
声明完成后调用 problem.validate()；validate 只检查。start/time::start/setTime 会锁定声明阶段，随后不允许
创建新命名场。面标量场使用 `createFaceScalarField`，面向量/面张量分别使用对应的显式类型名称。
Equation 和 History 应在 Case 有效期内创建、使用和销毁；Solver 不管理 PETSc 句柄。

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
留在对应 Physics 的算法实现中，通用层只提供数学入口。
逐项用法见 [DSL 与运行时](dsl-runtime-manual.md) 和 [数值配置](numerical-configuration.md)。

## 5. Solver 独立性

Heat、Transport、稳态 SIMPLE、瞬态 SIMPLE、PISO 都是独立注册函数。稳态 SIMPLE 使用一个主要
外循环；瞬态 SIMPLE 使用时间循环加自己的内层循环。两者不共享算法实现或模式开关。
non-orthogonal correction 可以作为各自压力步骤内的显式子循环。

PISO 使用时间循环、动量预测及压力校正，与瞬态 SIMPLE 的步内收敛循环分别实现。
RANS 模型 SA、k–omega、k–epsilon 分别拥有输运场、历史和闭合公式。流动求解器通过
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

先更新公共头和实现，再核对全部五个内置 Solver、RANS、外部示例和测试。保持当前配置规范
与实现一致；不要为了表面去重复引入 Manager、Provider、Registry 或深层继承。

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
make test-full        # 顺序执行单元、工作流、外部 Solver、RANS、时间阶和并行测试组
make test-workflow    # 新 Solver 工作流、时间序列、ParaView 读取
make test-external    # 仓库外用公共 include 构建 Solver（含负向 API 检查）
make test-rans        # RANS 方程与常数
make test-piso-coupling # PISO 校正、耦合缺陷与接受判据专用回归
~~~

标量时间阶与跨案例一致性分别由 `make test-scalar-time` 和 `make test-simple-parallel` 覆盖。
涉及并行和非正交时运行 `make test-mpi`、`make test-mpi-poiseuille` 及对应回归。
完整验证入口见 [验证与维护检查](validation.md)。架构测试检查 include 闭包和 Physics
越界依赖；外部 Solver 测试确认只用公共 include。验证报告中的历史结果只表示当时提交的
证据，不替代当前构建测试。

## 9. 通用优化的落点

后端根据方程身份保留稀疏模式及线性系统。矩阵系数有变化时更新当前 A；矩阵不变时只换
RHS。预条件器可以跨矩阵版本复用，但最终始终检查当前 A 的真实残差。复用预算、迭代增长
与失败重试见 [数值配置](numerical-configuration.md)，不能把旧 PC 误写为旧矩阵。

MPI 读取由 rank 0 完成全局解析、划分与分发，工作 rank 只构造 owned 加完整所需 halo 的
局部网格。释放全局与传输暂存降低生命周期峰值；保留全局 ID、面方向、patch 及至少三层
ghost。METIS 同时约束单元数和面关联工作量，并修复精确 owned 配额；无 METIS 时仍使用
图生长。初始化与划分耗时单独记录，不计作迭代速度收益。

COO 行列在模式建立时暂存，向量分量顺序共用工作 Vec，场重构复用必要的值快照；优化
在网格/场/离散执行/后端内部实施，不添加按算例名、网格族或特定 rank 数分支。
工程证据及抽象边界审查见 [优化报告](reports/flow-optimization-2026-10-07.md)。
