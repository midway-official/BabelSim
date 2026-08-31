# BabelSim 架构与维护边界

本文面向框架维护者。普通 Solver 作者先读 [开发指南](solver-development.md)：
作者组合方程和算法；维护者负责这些数学操作如何离散、存储和并行执行。

## 本轮解决的问题

- 公开 fvc 曾混有不负责通信的单面内核。现在公开入口统一为整场 evaluate/subtract，
  同步输入、执行离散并同步结果；局部内核仅供离散维护者使用。
- 矢量源 Field、压力参考和对角响应原先只有内部算法能使用，现在提供通用数学 API。
  SIMPLE 自己也使用这些入口，不保留另一套专用求解转发。
- Case 校验曾顺便关闭声明阶段，SIMPLE 构造会妨碍组合算法继续创建场。
  现在 validate 只校验；start/loop 才关闭声明。
- 应用曾直接包含 MPI 生命周期和格式转换，外部 Solver 需要修改内置选择逻辑。
  现在应用只提供显式 SolverEntry 表，公共 runApplication 负责启动。
- Field/Mesh 的存储数组曾可被普通调用者修改。现在原始访问、分区和构造维护入口
  位于 src/internal；公开接口只保留数学场操作和只读几何查询。
- 原 RunTime 混合时间推进、表达式解释、LDU 装配和数值工作区。
  现在 FVM 执行实现独立，不包含 RunTime 定义。
- SIMPLE 曾通过 internal/simple_discretization.h 调用框架中的专用逐面核。
  现在 Rhie–Chow 组合完全位于 SIMPLE，复用公开 fvc 面通量/扩散通量与内部面增量。
  初始化和日志也不再包含 runtime.h；算法全部源文件遵守相同公开边界。

## 实际层次和文件边界

不是每个层次都必须是一个类，也不是一条继承链。Mesh/Field 是共享的数学数据，
应用运行层组合 Case、时间与 FVM 执行；不存在通用算子反向依赖 SIMPLE 的关系。

| 层次 | 主要文件 | 职责 | 禁止依赖/承担 |
| --- | --- | --- | --- |
| 几何/拓扑 | mesh.h、core/mesh.cpp | 几何、拓扑、构造校验、固定布局和缓存 | MPI 调用、物理模型 |
| 场/边界 | field.h | scalar/vector/tensor、位置、数学赋值、边界 | 离散、历史、通信 |
| 数学描述 | fvm.h、fvc.h、fvm_expression.cpp | 轻量项与 lhs == rhs | 矩阵分配、通信、执行 |
| 方法 | methods.h | enum 默认格式和按场名覆盖 | 物理算法状态 |
| 通用离散核 | operators.cpp、internal/boundary_evaluation.h | 梯度/通量/扩散/时间项及边界离散 | Case、SIMPLE、RunTime |
| FVM 执行 | internal/fvm_execution.h、discretization/fvm_execution.cpp | 表达式解释、历史、装配、工作区、同步、数值诊断 | Case、应用、具体 Physics、RunTime 定义 |
| 离散/装配 | discrete_equation.h、assembly.cpp | LDU、owned 行及稀疏结构 | 温度/压力物理意义 |
| 代数 | algebra/ | 串行和分布式求解、全局点积 | Case、Physics |
| 并行 | parallel/ | 分解、halo、归约、MPI 校验 | PDE、算法停止策略 |
| 时间与运行 | runtime/runtime.cpp | 活动运行域、时间、后端生命周期 | LDU/装配公式和工作数组 |
| 数学 API 绑定 | runtime/solver_api.cpp | 将公开 solve/fvc/diagnostics 交给活动 FVM 后端；失败日志 | 离散实现、物理方程 |
| 应用启动 | application.h、runtime/application.cpp | 参数、MPI 初始化/销毁、显式分派与失败退出 | 内置物理分支 |
| Case/IO | case.h、io/、parallel_writer.cpp | 配置、命名场所有权、时间输出 | 求解物理方程 |
| 结果/后处理 | result.h、result_reader.cpp、postprocess.cpp | 文件验证、合并、VTU/PVD/Tecplot | MPI 运行依赖 |
| Physics/算法 | physics/heat、transport、simple | 方程、耦合、修正、收敛 | 存储、MPI、代数实现 |
| 应用文件 | apps/babelsim_solve.cpp、babelsim_post.cpp | 内置选择表/调用公共启动函数 | MPI API、格式细节 |

RunTime 仍选择当前唯一 FVM 后端，不是假称已具备任意后端插件能力。
但 FVM 执行者只接收网格、方法、线性配置、并行能力和步长，不知道应用时间循环。
solver_api.cpp 是已有公开函数的绑定实现，不是新增 Facade/Manager 类。

## Public Solver API

- case.h：命名场、物性与算法参数、声明阶段、时间循环和输出选择。
- field.h：数学边界、fill/assign/assignScaled/assignProduct/addProduct/evaluate。
- fvm.h：隐式项、已知源、轻量表达式；Scalar/VectorEquationDefinition。
- fvc.h：显式描述和整场 evaluate/add/subtract；单独包含即可使用。
- solver.h：solve、solveWithResponse、relaxed、referenceValue、只读 numericalMethods、diagnostics。
- simple.h：现成 SIMPLE 的步骤接口。
- application.h：SolverEntry、runApplication；无需 Registry 或注册宏。

外部 Solver 仅需 include/ 和预编译 libbabelsim.a，不需要 -Isrc、Eigen 头或 MPI 头。
普通编译器可编译 Solver；最终数值程序使用 MPI 链接器链接框架依赖。
result_reader 则连 MPI 库也不需要。include/ 中还保留供维护者使用的代数/并行头，
不能把“可包含”误认为“普通 Solver 应依赖”。

src/internal/field_access.h、mesh_access.h 是维护通道，不作为外部 SDK 分发。
它们内联访问原连续数组，不增加虚调用或元素复制。Mesh 的整体赋值、所有权/patch 修改
均受限；普通 API 无法用替换 Mesh 或改数组的方式破坏已经绑定的 Field。

## 生命周期与阶段

Case 拥有 Mesh、稳定地址的命名 Field、RunTime；销毁顺序为 FVM 后端/历史、Field、Mesh。
返回的引用以及表达式借用的场必须活得比表达式的求值更久。

1. 声明：读物性和控制，读取/创建全部场，定义边界，选择输出。
2. validate()：校验未消费参数，不关闭声明阶段，可重复执行。
3. start()：校验并关闭声明；算法驱动主入口显式调用。loop() 会隐式调用 start()。
4. 计算：可读取已有命名场，禁止再创建新场；loop 按物理时间推进和写出。
5. finish()：只用于成功结果；公共启动器在返回 0 后调用。

SimpleSolver 构造不再关闭 Case，组合算法可继续声明别的物理场。
同一个线程仍只允许一个活动 Case/RunTime；多区域 Case、子通信器和 restart 尚未提供。
外部 Solver 不应绕过 Case 创建短命未知量后让时间历史继续使用它。

## 数学与同步契约

公开 fvc 的所有参与进程必须按相同顺序调用；先同步输入，执行指定方法，再同步结果。
同位输入/输出别名会破坏直接求值顺序，因此 laplacian(T,T)、div(phi,U,U) 等被拒绝。
Field 的局部数学赋值不通信；下一个 fvm/fvc 操作负责所需同步。
Field::evaluate(function) 按位置定义初值/物性/源，函数必须与分区、调用次数无关。

source(a,F) 是已知体源，不是隐式线性化；支持 cell scalar/vector。
referenceValue(c) 是相容、常数零空间标量方程的定值规范，不是任意方程的点值约束。
它锚定固定全局参考单元；矢量方程使用该规范会被拒绝。
solveWithResponse 返回 V/aP 数学场，沿用 SIMPLE 当前缩放形式的欠松弛对角，
不暴露矩阵布局。线性收敛、数值健康、SIMPLE 外迭代收敛分别判断。

数学表达式只复制小描述符，不复制完整矩阵/Field。系数数组仍在离散方程建立时分配；
稀疏结构、线性工作区、几何缓存及算法工作场复用，不宣称“零分配”。
时间历史只在下一物理步推进一次，重复 solve 不推进时间层。

## SIMPLE 的私有边界

main 是算法流程；momentum/pressure 写数学方程与修正；create_fields 管理固定状态，
convergence 处理全局结果和日志，state.h 仅属于该算法。
Rhie–Chow 的数学组合位于 momentum.cpp；原专用逐面核及私有桥接头已删除。
FVM 只知道面矢量通量、给定重构梯度的扩散通量和几何面区域增量，不知道 SIMPLE。
动量响应和压力参考已改用通用 solveWithResponse/solve，不再是 SIMPLE 独占的内部入口。

保留项必须如实区分：simple.h 的 unique_ptr 是本算法私有所有权；state.h 的预分配
工作场是数学 Field，不是与 FVM 共享的内部存储。非正交 for 是数学迭代，不是存储循环。
包括初始化、诊断在内的全部 SIMPLE 源码都不访问 internal/、RunTime、MPI、矩阵或
原始 Field 数组。配置查询和日志分别经 numericalMethods/diagnostics::report 绑定到底层。
这两个函数不持有算法状态、不形成新的管理对象。

物理状态由 Case/调用者拥有；算法状态与数学中间量由 SIMPLE 拥有；数值执行工作区
由 FVM/代数/并行独占。禁止把专用算法原样搬入通用核心，也禁止把逐面存储内核原样搬入 Solver。

## 接口迁移

| 旧接口/文件 | 当前替代 |
| --- | --- |
| Field::data/mutableData/values/at/operator[] | 普通作者用数学场操作；维护者用 internal/field_access.h |
| Mesh 公开数组、setOwnership/setPatches/addPatchFace、整体赋值 | 只读几何 API；维护者用 internal/mesh_access.h |
| fvc::integratedNormalGradient 单面函数 | fvc::evaluate(fvc::normalGradient(p), result)；局部核留在 operators |
| SIMPLE 专用求解转发 | solveWithResponse(eq,rAU,relaxed(alpha))；solve(eq,referenceValue(value)) |
| detail::solve 与 ScalarEquationControl | 公开 EquationControl；矢量装配的可选响应仍是内部实现 |
| apps/solver_selection.* | application.h 的显式表；内置表在 babelsim_solve.cpp，外部表在自己的 main |
| RunTime 数值实现 | internal/fvm_execution.h 与 discretization/fvm_execution.cpp |
| result_reader.h 经 parallel_writer.h 获取结构 | 独立 result.h，无 MPI 包含链 |
| internal/simple_discretization.h、discretization/simple_discretization.cpp | SIMPLE 中用公开 fvc 组合；通用面区域增量保留在 FVM |
| SimpleSolver(run_time,fields,fluid,control) | SimpleSolver(fields,fluid,control)；活动运行域仍由调用者建立 |
| SIMPLE 中 RunTime::current().methods()/primary() | numericalMethods()/diagnostics::report() |

此前删除的 thermal.h/transport.h 专用库式入口、equation.h 兼容别名及重复线性控制不恢复。
存储维护 API 是有意的源码不兼容变化；仓库调用者和测试已迁移，不长期保留两套接口。

## 维护验收

- make test-architecture：含引号/尖括号项目头的依赖闭包、包含环、上下层边界；全部
  Physics 源码禁止跨入 src/internal、其他模块私有文件或运行/代数维护头。
  注入私有存储、运行头等违规样本，检查门禁本身确实拒绝，而不是只做正向扫描。
- make test-external：临时目录中只用公开 SDK 编译外部方程/耦合/矢量 Solver，
  1/2/4 进程解析解验证；完整 SIMPLE 源码仓库外重编译及 Poiseuille 回归；
  原始存储等十类负向编译；普通 g++ 链接结果读取器。
- make test / test-workflow：算子、生命周期、Heat、Transport、SIMPLE、自动时间输出。
- make test-mpi / test-mpi-poiseuille：halo、分布式代数、非正交 fvc、SIMPLE 和结果差异。

这是可执行的边界与数值证据，不是 ABI 或任意外部程序正确性的形式化证明。
当前未实现 FEM/FDM 后端、张量未知量求解、一般约束系统、多区域或动态插件。
