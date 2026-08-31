# BabelSim 架构与维护边界

本文面向框架维护者；普通 Solver 作者先读 [solver-development.md](solver-development.md)。
精简的标准是：同一职责只有一个实现来源，依赖方向清楚，数值与生命周期契约不被削弱，
而不是把必要的工作区删掉，或者把几百行代码藏进另一个 Manager。

## 1. 本轮审计结论与处理

| 审计发现 | 处理 | 保留的数值能力 |
| --- | --- | --- |
| Heat/Transport 的 Case 入口与库式入口重复方程、时间循环和结果类型 | 删除 thermal.h、transport.h 及两个 transient_* 实现；内存测试直接用 fvm/solve | 常数/Field 系数、时间历史、串并行求解 |
| equation.h 只是旧名称转发，数学 Equation 与 LDU 容易混淆 | 删除兼容别名；唯一存储定义为 discrete_equation.h 中的 DiscreteEquation | 原来的连续 LDU 数组、跨分区面系数 |
| SimpleControl 与 RuntimeControl 重复保存线性配置 | SimpleControl 只保留算法参数；线性配置仅归 RuntimeControl | 相同求解方法、容差和预条件器 |
| 通用 operators.cpp 含 Rhie–Chow 专用核 | 专用核归 simple_discretization.cpp | 原公式、面遍历顺序与非正交法向梯度 |
| 专用离散文件反向定义 SimpleSolver 私有工作区的方法 | 工作区与编排回到 physics/simple；专用数值层只接收数学方程/Field | 固定分配并复用工作场 |
| 启动选择逻辑打包进核心库，声明位于公开数值头目录 | 声明/实现归 apps，只有启动器链接选择表 | 一个 Solver 函数加一项显式选择 |
| 四种方法覆盖重复相同解析逻辑 | 一个文件内模板检查 Field 名重复，保留独立 enum | 原配置语法、类型检查、错误诊断 |
| 非正交次数有单用途计数类，多个函数只做一次转发 | 使用普通修正循环和已有 fvc/Field API | 原求解次数与步骤顺序 |

没有增加新的框架类、继承层、Factory、Registry 或执行引擎。SIMPLE 的 State 不是新增抽象，
只是把原有私有状态放回所属算法目录。

## 2. 实际层次与文件边界

这不是一条严格的继承栈。Mesh/Field 是共享对象；RunTime 将数学操作、离散、代数与并行组合起来。
应按以下职责判断依赖是否合理，而不是只看目录深度。

| 层次 | 主要文件 | 负责什么 | 不负责什么 |
| --- | --- | --- | --- |
| 几何/拓扑 | mesh.h、vector.h；core/mesh.cpp | cell/face/patch、几何缓存、整数索引及局部分布映射 | MPI 调用、物理模型 |
| 场/边界 | field.h | 标量/矢量/张量存储、位置、数学边界、通用场赋值 | 时间历史、离散方程、消息通信 |
| 方法与数学描述 | methods.h、fvm.h、fvc.h；fvm_expression.cpp | 格式选择、轻量运算描述、lhs == rhs | 创建矩阵、求解或通信 |
| 物理/算法 | physics/heat、transport、simple | PDE、物性、耦合顺序、收敛判断 | MPI、稀疏矩阵、逐 cell/face 存储处理 |
| 通用离散 | operators.h；discretization/operators.cpp | 梯度、散度、对流、扩散、时间项与 BC 的离散 | SIMPLE 状态或具体热模型 |
| 专用数值步骤 | internal/simple_discretization.h；simple_discretization.cpp | 动量响应、压力参考、Rhie–Chow 核；复用通用离散/执行接口 | 持有 SimpleSolver、Case 或算法工作区 |
| 离散表示与装配 | discrete_equation.h、assembly.h；assembly.cpp | 唯一 LDU 表示、owned 行装配、稀疏结构缓存 | 物理方程选择、停止策略 |
| 串行代数 | linear_solver.h；algebra/linear_solver.cpp | A/b/x、Krylov、预条件器 | Mesh、Field、装配、物理模型 |
| 分布式代数 | distributed_solver.h；distributed_solver.cpp | owned 行、跨分区系数、halo matvec、全局点积 | 温度/压力等物理含义 |
| 并行 | parallel.h、mpi_support.h；parallel/ | 分解、halo、归约、MPI 生命周期检查 | PDE、SIMPLE 停止规则 |
| FVM 执行协调 | runtime.h；runtime/runtime.cpp | 时间历史、解释表达式、调用离散/代数、同步和诊断 | 读取物性、选择 Heat/SIMPLE |
| Case/IO | case.h、各 IO 头；io/、parallel_writer.cpp | 输入校验、命名场所有权、按时间写结果 | 离散公式、算法外迭代 |
| 应用 | apps/babelsim_solve.cpp、solver_selection.*、babelsim_post.cpp | 参数、MPI 进退、显式分派、独立后处理 | 物理方程与求解器工作场 |

分布式代数知道 Mesh 的 owned/ghost 映射和离散面的代数耦合，是当前 FVM 并行实现的必要契约，
不等于知道 Navier–Stokes。串行线性求解器不需要这些内容，已移除它经 assembly.h 间接引入的几何依赖。

SIMPLE 专用数值步骤可以调用通用执行接口；通用执行和通用算子不能反过来包含 SIMPLE。
专用步骤也不能包含 physics/simple/state.h。物理算法负责准备工作场，数值核只消费所需数据。

## 3. Solver API 与维护接口

普通 Solver 入口仍是：

- case.h：场、物性、时间和输出的统一问题入口；
- solver.h：solve、显式求值和全局 diagnostics；
- simple.h：使用现成 SIMPLE 时的算法步骤；
- Field、Boundary、fvm/fvc：由上述头提供的数学对象与表达。

include/babelsim 是可包含的头文件集合，不代表其中每个文件都应由普通 Solver 使用。
runtime.h、discrete_equation.h、operators.h、assembly.h、linear_solver.h、parallel.h 是维护接口。
目前用明确命名、包含闭包和源码检查区分，而不是再造一套包装 Field。

所有 include/babelsim 头的项目依赖均留在 include/，不再从公开头包含 src/internal 文件。
src/internal 只放框架内部的方程控制和专用数值桥接声明。
SIMPLE 私有状态在 src/physics/simple/state.h；应用选择声明在 src/apps/solver_selection.h。
动量、压力、主步骤及私有状态头均不包含 RunTime 定义，只有初始化与日志实现直接使用它。

测试构造一个小网格后可显式创建 RunTime，直接调用同一个 fvm/solve；不再为测试保留另一套
固定物理方程库入口。Case 工作流测试实际运行生产 Heat/Transport 入口，防止单元测试与入口脱节。

## 4. 生命周期与配置的唯一来源

Case 拥有局部 Mesh、稳定地址的命名 Field 和 RunTime。析构顺序为：
执行对象/历史/后端先释放，Field 再释放，Mesh 最后释放。外部引用不得超过 Case 的生存期。

带初值创建的中间场与输入场使用同一所有权机制，不默认写出。
新场在进入计算前声明；Case::validate 会检查参数消费并关闭声明阶段，因此不是 const 操作。
目前 SIMPLE 构造会调用该检查；新增耦合场/物性仍须在构造 SIMPLE 前声明/读取。
这是当前工作流约束，不应被隐瞒成“任意顺序都可组合”。

SIMPLE 保留固定分配的状态：

| 状态 | 内容 |
| --- | --- |
| 物理场引用 | U、p、phi，由 Case 或显式调用者拥有 |
| 算法场 | pPrime、rAU、phiHbyA、UPrevious |
| 数值缓存 | 压力梯度、面插值、散度 |
| 执行顺序与结果 | 步骤枚举、迭代次数、线性/健康/外迭代结果 |

SimpleControl 只保存最大次数、非正交次数、欠松弛和算法停止容差。
RuntimeControl 保存 Methods、TimeControl、scalar_solver/vector_solver。
不得再从 SimpleControl 复制一套“可能与实际运行不一致”的线性配置。

## 5. 数学、数值和性能契约

lhs == rhs 形成描述符列表；加减数学项不复制 LDU 或整场。
真正 solve 时才选择方法并调用离散核；边界条件进入同一次离散。
ScalarEquationDefinition/VectorEquationDefinition 是数学描述；
ScalarDiscreteEquation/VectorDiscreteEquation 是明确命名的代数系数载体。

fvm::source 是已知体源，不是源 Field 的隐式线性化；source(a,C) 复用描述符系数。
数学符号 A/b/x、U/p/phi/rAU 保留原命名，不做没有语义收益的长名替换。

几何、稀疏结构、线性求解工作区和 SIMPLE 整场缓存继续复用。
当前表达式列表仍有小容器分配，Runtime 每次构造离散方程仍会分配其系数数组；
本轮没有宣称“零分配”，也不通过删除安全检查或缓存来缩短代码。
时间历史只在进入下一物理步时推进一次，内迭代重复 solve 不修改旧时间层。

MPI、残差归约和通信顺序未改变。健康状态、线性收敛与外迭代收敛依然分别检查。
非正交循环改用普通 for，但正交求一次、非正交追加修正的次数与顺序不变。

## 6. 文件与编码规范

- 一段生产 PDE 不同时维护 Case 版和库函数版；确有外部复用需求再设计共享入口。
- 状态定义随所属算法放置；底层只接收最小数学/代数输入，不包含上层私有类。
- 匿名命名空间收纳文件私有函数，避免把未声明的辅助符号导出到整个库。
- 类成员使用 m_ 前缀；数学符号和简单值结构体字段保持直观命名。
- 每个成员有明确初始化；借用的 Mesh/Field 必须覆盖所有使用者的生命周期。
- 头文件直接包含必要声明，不通过偶然的传递 include 获取类型。
- 注释解释数学、生命周期和性能理由；不用计数器类替代普通数值循环。
- 不因精简改动浮点容差、预条件器、边界公式或进程间停止条件。
- 无兼容需求的旧 API 必须迁移所有仓库调用者后删除，不留只做别名转发的长期层。

## 7. 删除 API 的迁移表

| 已删除/调整 | 现在使用 |
| --- | --- |
| thermal.h、solveTransientHeat、solveHeatStep、HeatResult | 普通入口写热方程；内存测试用 RunTime + solve，返回 SolveResult |
| transport.h、solveTransientScalarTransport、ScalarTransportResult | 同一套 ddt/div/laplacian/solve |
| equation.h、Equation<T>、ScalarEquation、VectorEquation | discrete_equation.h、DiscreteEquation<T>、ScalarDiscreteEquation、VectorDiscreteEquation |
| SimpleControl.velocity_solver / pressure_solver | RuntimeControl.vector_solver / scalar_solver |
| simpleRunTimeControl | 显式配置 RuntimeControl，不再复制两份控制 |
| solvers.h | apps/solver_selection.h，仅供启动器链接 |
| src/internal/simple_state.h | physics/simple/state.h，只有该算法实现包含 |
| 两个单独的 equation_control 头 | internal/equation_control.h，统一内部方程控制 |

仓库调用者和测试均随迁移调整；这些变更对直接使用旧维护接口的外部代码不保持源码兼容。
普通 Case、Field、fvm/fvc、solve 和 SIMPLE 五步调用不变。

## 8. 验证及审计范围

make test-architecture 检查全部 src/include 项目头依赖：
缺失 include、包含环、公开 Solver 头带入实现、底层依赖 Physics/应用、
通用算子依赖 SIMPLE、串行代数依赖几何，以及被删重复入口重新出现。
它同时由 make test 和 make test-workflow 运行；没有复制两份架构检查逻辑。

该检查是包含图与源码级约束，不是所有函数调用关系、数值正确性或 ABI 的形式化证明。
Field 仍公开维护者所需的原始数据接口，故不能宣称 C++ 类型系统已完全禁止误用。
RunTime 当前是 FVM 执行协调器，尚不是任意离散后端的纯运行时。
FDM/FEM、一般材料模型、restart、多区域和块耦合也不是本轮增加的功能。

每轮精简后运行 make test、make test-workflow、make test-mpi、make test-mpi-poiseuille。
源码检查通过之后仍要核对数值、并行停止次数、全局守恒和真实后处理读取结果。

### 8.1 对剩余实现细节的逐项判断

不能把“搜索到 C++ 语法”直接等同于抽象泄漏。复查 src/physics 与 SIMPLE 公开头后，
按职责作如下判断：

| 出现位置 | 内容 | 判断与边界 |
| --- | --- | --- |
| simple.h 私有成员 | unique_ptr<State> | 保留明确所有权和不完整类型隔离；Solver 作者无需创建/销毁 State，没有新增包装层 |
| simple/create_fields.cpp | RunTime、Field 构造及初始化 | 属于现成算法的维护实现；普通入口只构造 SimpleSolver(problem)，动量/压力文件不包含 RunTime |
| simple/convergence.cpp | RunTime::current().primary() | 仅限制日志输出；健康、线性与外迭代判据仍通过全局 diagnostics 计算，不以本地值控制进程 |
| simple/convergence.cpp | 三个速度分量的 for | 数学结果合并，不是 cell/face 存储遍历；保留普通循环 |
| simple/pressure.cpp | 非正交修正 for | 算法迭代本身，保留顺序与求解次数；不是底层内存操作 |
| simple/momentum.cpp、pressure.cpp | fvc::evaluate、assign、assignProduct | 高层场运算并复用预分配缓存；不读取 data 或访问 cell/face 索引 |
| simple/momentum.cpp | detail::applyMomentumInterpolation | 框架内专用数值接口，仅供算法维护者调用；拓扑遍历和公式实现归专用离散文件 |
| internal/equation_control.h | 可选 mobility 指针 | 同步方程求解的借用输出，不拥有 Field；由专用数值入口传入，普通 Solver API 不暴露此控制对象 |

当前 src/physics 没有 MPI API、通信缓冲、Eigen、稀疏矩阵操作、原始 Field data 访问，
也没有逐 cell/face 手写循环。上述保留项不能被宣传为“整个 physics 目录只剩 PDE”；
该目录也包含 SIMPLE 的算法维护实现。普通 Solver 主入口与这些内部实现的认知边界必须保持明确。
