# BabelSim DSL 审计报告

审计对象：include/babelsim 的公共 Physics API、src/physics 下 Heat、Transport、稳态 SIMPLE、
瞬态 SIMPLE、RANS，以及 Case/Runtime/IO/离散边界。审计结论以当前源码和回归测试为准。

## 1. 保留的优点

- Field 与 Mesh 解耦，cell/face 位置和边界条件由 Field 描述；
- math 返回已经计算的显式场，equ 立即向绑定 Equation 加离散项；
- 普通 C++ 控制流保留，SIMPLE 的算法步骤在 main 中可见；
- faceFlux 从实际离散 Equation 取得一致面通量；
- Physics 不暴露 CSR、MPI、Eigen、CUDA、halo 和网格遍历；
- RANS 模型可独立拥有输运场、历史和闭合；
- 通用 monitor 只接收名称和值，不包含 SIMPLE 或 RANS 收敛逻辑。

## 2. 已解决的语义混乱

| 原问题 | 当前规则 |
| --- | --- |
| math 混入复制/历史/生命周期 | 场值复制直接用 Field 值语义；历史使用 time::History |
| correction 名称不能说明边界含义 | field::homogeneousLike 明确表示同网格齐次边界场 |
| Equation 创建语义不清 | equ::createEquation(unknown) 明确绑定未知量 |
| solve 再次传入未知量容易不一致 | equ::solve(equation, control) 使用已绑定未知量；旧重载保留兼容 |
| response 需要理解 SIMPLE 才能猜 | 删除复合响应 API；由 geometry::cellVolumes / Equation::diagonal 直接写 V/aP |
| interpolate 特殊重载隐藏 Rhie–Chow | interpolate 永远只插值；Rhie–Chow 在 SIMPLE main |
| load/create 重载混淆 | scalarField 等只加载，create*Field 只创建，existing*Field 只查找 |
| 初始化值可能被静默忽略 | 重复 create 或 load/create 混用直接抛错 |
| methods 反复解析 | Case 构造时读取一次，problem.methods() 只读查询 |
| 线性控制签名与行为不一致 | readLinearControl(problem, field) 支持 scalarSolver.Field/vectorSolver.Field |
| Physics 解析 token | Parameters 提供 word、boolean、number、positive、nonnegative、fraction、integer |
| 输出字段必须改 C++ | output.bs 提供 writeFields/excludeFields，程序场默认关闭 |
| Physics 直接 return 2 或控制打印 | SolverResult 返回状态，monitor 通用报告，application 做退出映射 |

## 3. 目标职责

Case 拥有网格、命名场、参数和输出；Field 提供值和边界；math 计算显式场数学；equ
组装离散 Equation；time 推进时间和历史；diagnostics 只观察；monitor 只打印 metric；
Runtime/FVM/后端负责同步、离散执行和线性代数；Physics 负责 PDE、算法循环、修正和收敛。

配置分类为 physics（物性和模型）、methods（空间/时间格式）、solution（线性和算法数值
控制）、control（物理时间）、output（结果选择）。所有字典项按类型读取并追踪消费。

## 4. SIMPLE 和 RANS 的边界

稳态 SIMPLE 与瞬态 SIMPLE 是两个独立 Solver，分别只有自己的 main.cpp。稳态一个主要
外循环，瞬态是时间循环加自己的 SIMPLE 内循环；压力 non-orthogonal 修正可在各自 main
中作为显式子循环。Rhie–Chow 只由各自 main 用公开 math/equ/Field 操作实现。

SA、k–omega、k–epsilon 各自实现输运和闭合。SIMPLE 只看到 effectiveViscosity、应力余项
和 TransportResult，不知道具体模型公式。模型常数属于 physics，松弛和容差属于 solution。

## 5. 迁移对照

| 不再作为新代码入口 | 当前入口 |
| --- | --- |
| math::copy/history/saveOld | Field 值复制、time::History::save |
| math::correction/createHomogeneousField | field::homogeneousLike |
| equ::Matrix | equ::Equation、equ::createEquation |
| Equation::volumeScaledInverseDiagonal | 删除；使用 `geometry::cellVolumes(mesh) / equation.diagonal()` |
| 含初值的 scalarField/vectorField 重载 | createScalarField/createVectorField |
| faceField/faceFlux 生命周期混合入口 | createFace*Field、equ::faceFlux |
| 网格内部几何数组对 Physics 不可见 | geometry::cellVolumes/faceAreaVectors/faceAreas 等 Field 视图 |
| Solver 自己读取 methods.bs | problem.methods() |
| entry().tokens | Parameters typed getter |
| Physics return 数字退出码 | SolverResult |
| simple.solve() 黑盒 | main 中显式 SIMPLE 步骤 |
| `LinearSystem` / `assemble()` | `SparseAssembly::update()` + `matrix()` 与显式 `assembleSource()` |
| `assembleMatrix()` 临时矩阵包装 | 可复用的 `SparseAssembly`（一次建立结构，多次更新值） |
| `solve(Eigen::SparseMatrix, ...)` 一次性后端入口 | `PreparedLinearSolver::compute/factorize/solve` |
| `Case::loadMethods()`、`loadMethods(Case&)` | 构造时加载，使用 `Case::methods()` |
| `Case::createFaceField()` | 按值类型命名的 `createFaceScalarField()` |
| `Case::loop()` 隐藏时间循环 | `time::start(problem)` + 显式 `TimeStepper`/算法循环 |

## 6. 后端静态清理

本轮按“先找调用者，再删除入口”的原则检查了 `src/backend`、`src/algebra` 和实现头文件。
`SparseAssembly` 仍被 Eigen/MPI 后端用于缓存 CSR 结构；`assembleSource` 仍被串行和分布式
线性求解路径使用；`PreparedLinearSolver` 仍被 Eigen 后端用于矩阵准备和重复求解。因此这些
能力属于后端稳定契约，不能因为 Physics 不直接调用就删除。

确认没有生产调用后，删除了三类只提供临时包装或重复语义的接口：

1. `LinearSystem` 与 `assemble()/assembleMatrix()`：它们每次重新创建矩阵并把矩阵/右端捆成
   一个后端结构，绕过了可缓存的 `SparseAssembly`。测试已改为显式更新装配器、取得矩阵并
   单独装配源项，仍覆盖标量/向量系数等价性。
2. `PreparedLinearSolver` 之外的 Eigen 一次性 `solve(A,b,x,config)`：该函数内部立即构造、
   准备并销毁求解器，不能复用预条件器，也没有 Physics 或后端调用者。测试改用同一个
   `PreparedLinearSolver` 完成相同的 CG、BiCGSTAB、AMG 和无预条件器验证。
3. 向量按值返回的 `assembleSource(equation)`：只有带输出参数的版本被 Eigen/MPI 路径使用；
   删除按值重载避免不必要的临时分配。

同时删除了 Case 中未被任何 Solver、测试或应用使用的历史兼容入口：`createFaceField`、
`loadMethods`（成员与自由函数）、单参数 `setTime` 和 `Case::loop`。它们不改变数值语义，
却会重新引入类型不明确、方法重复解析或隐藏时间推进。保留的 `write(Case, TimeStepper)`
和显式 `Case::setTime/write` 仍是当前 Physics/示例实际使用的 I/O 入口。

后端公开实现头仍由 architecture test 标记为 implementation-only；Physics 的包含闭包没有
引入 Eigen、CSR、MPI、DistributedLinearSolver 或 DiscreteEquation。后续若要进一步缩小
安装包，可把这些已标记的实现头从安装清单迁入 `src`，但这属于打包边界调整，不在本轮数值
重构范围内。

## 7. 验收

architecture_test 检查层次和包含边界；case_io/lifecycle 检查配置消费和场生命周期；
procedural/operators/backend 检查数学合同；Heat、Transport、SIMPLE、时间历史和 RANS
脚本检查实际流程；workflow、external 和 MPI 检查应用接入与并行一致性。当前验证命令及
结果见 [validation.md](../validation.md)。
