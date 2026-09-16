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
| response 需要理解 SIMPLE 才能猜 | Equation::volumeScaledInverseDiagonal 明确是 V/aP |
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
| equ::response | Equation::volumeScaledInverseDiagonal |
| 含初值的 scalarField/vectorField 重载 | createScalarField/createVectorField |
| faceField/faceFlux 生命周期混合入口 | createFace*Field、equ::faceFlux |
| Solver 自己读取 methods.bs | problem.methods() |
| entry().tokens | Parameters typed getter |
| Physics return 数字退出码 | SolverResult |
| simple.solve() 黑盒 | main 中显式 SIMPLE 步骤 |

## 6. 验收

architecture_test 检查层次和包含边界；case_io/lifecycle 检查配置消费和场生命周期；
procedural/operators/backend 检查数学合同；Heat、Transport、SIMPLE、时间历史和 RANS
脚本检查实际流程；workflow、external 和 MPI 检查应用接入与并行一致性。当前验证命令及
结果见 [validation.md](../validation.md)。
