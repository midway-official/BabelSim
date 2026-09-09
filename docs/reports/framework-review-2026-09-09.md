# BabelSim 全代码审查与开发体验评估

审查日期：2026-09-09。基线：`e112e597172b7fca5956504c4905e19add13d252`（非结构六面体网格迁移后）。

本报告对应两轮检查：第一轮追踪数值算法和并行实现；第二轮重新覆盖全部核心模块，并补查网格、字段语义、运行生命周期、配置、应用入口、结果读写和后处理。结论来自当前代码、实际完成的测试、独立诊断程序和模型方程对照。审查不等同于形式化正确性证明，也没有修改生产代码来掩盖问题。

**总体判断：BabelSim 已经是一个分层较清楚、可独立开发求解器的紧凑型有限体积框架，具有研究、教学和中小规模 CFD 开发价值。层流稳态能力有实际收敛基准支撑；RANS 尚有明确正确性缺陷；瞬态完整物理验证不足；MPI 主路径有效，但特定算子组合不一致，AMG 和网格启动流程具有明显规模边界。当前不能作为已充分验证的通用工业多物理平台使用。**

## 1. 优先处理的审查发现

P1 表示会产生错误方程、错误结果或错误成功状态，应优先修复；P2 表示特定使用路径下的数据或边界行为缺陷。下面七项均有具体触发条件，不能用现有测试通过替代修复。

### F1 · P1：RANS 动量方程缺少变黏度应力的一部分

位置：[稳态 momentum.cpp:16](/home/midway/BabelSim/src/physics/simple/momentum.cpp:16)、[瞬态 momentum.cpp:20](/home/midway/BabelSim/src/physics/transient_simple/momentum.cpp:20)。

两条路径都只把层流扩散替换成 `laplacian(muEffective, U)`，实现的是 `div(muEffective grad(U))`。对于不可压缩 Boussinesq 涡黏性闭合，偏应力还包含转置梯度贡献；当湍流黏度空间变化时，`div(muEffective grad(U)^T)` 一般不为零。应力形式可对照 [TMR 的 Wilcox 模型定义](https://tmbwg.github.io/turbmodels/wilcox.html)。

一个直接反例是 `U=(y,0,0)`、`muEffective=1+x²`：速度散度为零，当前分量 Laplacian 为零，遗漏项却为 `(0,2x,0)`。该项具有非零旋度，不能普遍通过重新定义压力消去。这是数学上的方程不完整，不只是离散精度或收敛速度问题。

此发现**不把各向同性 k 项的压力吸收约定单独判错**；代码明确选择 Wilcox1988m，应尊重对应压力/生产项约定。真正的问题是变涡黏度的偏应力没有完整进入稳态和瞬态动量方程。

建议：提供通用张量散度或完整黏性应力离散能力，由私有物理模块组织闭合；以非均匀黏度、无散速度制造解检查各分量及边界通量。不要把某个湍流模型专用公式塞进框架公共接口。

### F2 · P1：派生黏度场没有继承正确的物理边界值

位置：[model.cpp:73](/home/midway/BabelSim/src/physics/RANS/model.cpp:73)、[SA 派生场计算](/home/midway/BabelSim/src/physics/RANS/spalart_allmaras.cpp:134)、[Field 值操作](/home/midway/BabelSim/include/babelsim/field.h:125)。

`fill/addScaled/evaluate/assign` 改变场值，保留目标场原有边界条件。这个通用语义本身可以成立，但 RANS 把内部涡黏度写入新场后，没有同时构造对应边界值；`muEffective` 仍使用默认 zeroGradient。

已复现：分子动力黏度 `0.001`，SA 壁面 `nuTilda=0`，内部 `nuTilda=0.01`。壁面有效黏度应回到 `0.001`，实际面插值得到 `0.00836425`。这会改变壁面黏性通量；SA 的派生扩散系数也有同类边界语义风险。此处不是要求任意字段自动复制 BC，而是模型必须显式提供派生物理量的正确边界行为。[标准 SA 壁面约束](https://tmbwg.github.io/turbmodels/spalart.html) 与此判断一致。

证据：[sa-boundary.cpp](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/sa-boundary.cpp)、[日志](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/sa-boundary.log)。

建议：定义模型派生场的边界计算规则，分清分子量、湍流量和总量；验证面上系数、壁面剪切、湍流变量扩散通量，而不只检查单元内部数值。

### F3 · P1：裁剪后的零变化会掩盖湍流方程未收敛

位置：[k_epsilon.cpp:58](/home/midway/BabelSim/src/physics/RANS/k_epsilon.cpp:58)，k-omega 和 SA 有相似的“求解、裁剪、计算变化量”结构。

当前线性求解状态描述裁剪前的线性系统；随后把场截到下限，再用截断后的场计算 `relativeChange`，最后返回此前的成功状态。因此，“线性求解成功 + dTurb=0”不能证明最终保存的场满足湍流方程。

已复现：静止流、`rho=1`、`mu=0.001`，`k=epsilon=kMin=epsilonMin=0.001`。连续五次更新均输出 `converged=1, dTurb=0`，而裁剪后 k 方程单位体积残差为 `-0.001`：梯度和生产项为零，耗散项没有平衡。使用较高但合法的下限是为了明确暴露判据缺陷；默认小下限也不消除这个逻辑问题。

证据：[rans.cpp](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/rans.cpp)、[日志](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/rans.log)。

建议：增加通用隐式线性源项能力，使破坏项采用合适的 sink 线性化；报告裁剪数量和幅度；对最终更新场计算非线性残差或明确的可接受性判据。残差、更新量和限制器活动应分别报告。

### F4 · P1：串行线性求解器会把相对残差 1 报告为收敛

位置：[linear_solver.cpp:96](/home/midway/BabelSim/src/algebra/linear_solver.cpp:96)，特别是 110–115 行的 roundoff 判据；同文件固定绝对 breakdown 阈值也是触发链的一部分。

已复现：一阶、条件数为 1 的系统 `A=[1e-20], b=[1e-20]`，初值零，禁用预条件器，`atol=1e-30, rtol=1e-12`。正确解为 1；CG 和 BiCGSTAB 都返回 `converged=1, iterations=0, x=0, relativeResidual=1`。

终结检查使用 `64*epsilon*max(1, ||b||+||A||*||x||)` 放宽目标，给很小的有量纲系统引入约 `1e-14` 的绝对下限，覆盖了用户显式要求。MPI 路径没有相同的放宽规则，串行/并行状态语义也不一致。

证据：[scale.cpp](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/scale.cpp)、[日志](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/scale.log)。

建议：使 breakdown 判断具有尺度一致性；成功必须满足公开容差合同。确实受舍入限制时，返回独立的停滞/精度受限状态，不应静默当作 Converged。增加同一线性系统跨数量级缩放和串行/MPI状态一致性测试。

### F5 · P1：Green–Gauss 复合算子需要的 halo 超过默认两层

位置：[Green–Gauss 修正](/home/midway/BabelSim/src/discretization/operators.cpp:240)、[Laplacian 内部梯度](/home/midway/BabelSim/src/discretization/operators.cpp:1106)、[默认 halo 宽度](/home/midway/BabelSim/include/babelsim/parallel.h:36)。

修正 Green–Gauss 对 owned 单元的梯度需要两层输入邻居。但复合面算子还会消费第一层 ghost 单元的修正梯度，等效依赖继续向外扩一层。内核直接构造的中间梯度没有经过一次独立同步，默认两层 halo 不足以保证结果完整。公开 `math::grad` 的结果同步与“其他算子内部临时梯度”是两条不同路径。

基于现有 parallel_math_test 构造非线性场和非仿射扭曲网格，确保串行和各 rank 物理边界完全一致，四进程结果如下：

| 梯度/halo | 串行与 MPI 最大绝对差 |
|---|---:|
| Green–Gauss，默认 2 层 | 2.46513e-4 |
| 同一 Green–Gauss，3 层 | 0 |
| LeastSquares，2 层 | 8.88178e-16 |

证据：[gg.cpp](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/gg.cpp)、[gg.log](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/gg.log)、[三层对照](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/gg3.log)、[LS 对照](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/ls.log)。诊断程序保留有限性检查，打印偏差；它以零退出不代表该路径通过正确性验收。

建议：在 FVM 执行层分解并同步中间量，或声明并满足算子组合真正需要的 halo 深度。保持离散内核不直接依赖 MPI。测试必须包含非线性场、扭曲几何和多种方法组合；线性场和仿射网格会掩盖此类缺陷。

### F6 · P2：inletOutlet 在入流时没有向扩散项提供入口 Dirichlet 条件

位置：[标量扩散边界装配](/home/midway/BabelSim/src/discretization/operators.cpp:1271)，显式扩散和向量路径需一起检查。

对流项按流向把 inletOutlet 入流侧当作给定入口值；扩散装配只处理 FixedValue/FixedGradient，把 InletOutlet 落到无贡献分支。因此，同一个边界在同一对流扩散方程中的解释不一致。

单个单位六面体反例：左面入流通量 -1、入口浓度 1，右面出流 +1、浓度固定 0，扩散率 1。左面设 fixedValue 时得到 `C=0.6`；仅改为入口值相同的 inletOutlet，得到 `C=0.333333`，两次均报告线性收敛。缺少的正是左边界扩散贡献。

证据：[inletoutlet.cpp](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/inletoutlet.cpp)、[日志](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/inletoutlet.log)。

建议：建立基于当前面通量的边界求值结果，让对流、扩散、梯度共享同一物理解释，并明确无通量上下文时的行为。添加同一 inletOutlet 面入流/出流切换的组合方程测试。

### F7 · P2：后处理会把旧结果静默映射到修改后的网格

位置：[postprocess.cpp:194](/home/midway/BabelSim/src/io/postprocess.cpp:194)、[result_reader.cpp](/home/midway/BabelSim/src/io/result_reader.cpp)。

后处理读取当前 case 网格，用全局单元数匹配历史结果。CSV 中虽带单元坐标，读取后没有用它验证当前几何，也没有对应网格指纹。

已复现：保留 heat 原结果，把网格所有顶点 x 坐标平移 +10，单元数量及编号不变；`babelsim-post ... -time original -format vtk` 仍成功生成输出。结果值被放在新位置，没有提示几何已不一致。

证据：[故意不匹配的输入与结果](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/post-case/case.bs)、[日志](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/post.log)。

建议：保存网格指纹或每次运行的网格快照，校验历史结果的几何来源；同时记录配置、版本、收敛状态和完整写入标记。现有全局 ID、分片完整性、时间及有限性检查有价值，应继续保留。

## 2. 软件性质和实际抽象层次

BabelSim 使用 C++17、MPI、Eigen3 和 GNU Make，采用 MIT 许可，提供静态库和 `babelsim-solve`、`babelsim-post` 两个应用入口。它已超出单个 SIMPLE 程序：外部项目可以依赖公开方程接口定义独立求解器，并沿用同一运行、离散和并行系统。

网格是显式非结构六面体：每个单元八顶点、六个四边形面。当前“二维算例”以薄的三维六面体层表示。任意编号和一般连接关系受到支持，但这不等于支持四面体、任意多面体、动态网格或任意维数。

| 层次 | 主要职责 | 当前边界评价 |
|---|---|---|
| Mesh / Field / Vec3 / Tensor3 | 几何、拓扑、场位置、边界和数值对象 | 所有权与普通物理代码隔离；字段值与边界状态的合同不够醒目 |
| eqn / math / Methods | 隐式方程项、显式算子和离散方法选择 | 物理表达清楚；支持范围集中在对流扩散类方程 |
| FvmExecution | 表达式解释、时间历史、同步调度、离散方程组织 | 合理的框架中枢；组合算子模板依赖没有充分显式化 |
| Operators / LDU assembly | 面通量、梯度、矩阵系数和源项 | 物理中立；大文件集中多种类型/BC分支，维护风险较高 |
| ComputeBackend / Eigen assembly / algebra | 稀疏装配、预条件、Krylov、场同步与归约 | 粗粒度替换接口有真实测试；不同后端的数值语义仍需统一 |
| ParallelContext / HaloExchange | 分区、owned/ghost、面关系、通信 | 普通 solver 无需接触 MPI；规模与模板合同有局限 |
| RunTime / Case / Application | 时间循环、字段注册、配置、结果、错误和 MPI 生命周期 | 减少样板代码；单活动上下文及注册表限制复杂耦合 |
| physics | heat、transport、稳态/瞬态 SIMPLE、私有 RANS | 算法归属正确；稳态/瞬态有重复公共计算步骤 |

关键调用关系是：应用选择物理入口，物理构造 `eqn/math` 表达式，运行层把表达式交给 FVM 执行层，离散层形成代数系统，后端求解并更新字段。Case 提供字段和配置，MPI 细节留在框架/后端内。

当前 physics 没有直接包含 MPI、Eigen、LDU/CSR 或内部字段存储；RANS 模型常数在 physics 字典中，模型工厂和模型细节保留在私有目录。这是值得保持的结构。`solveWithResponse` 提供通用对角响应，并不因注释出现 SIMPLE 或 aP 就构成物理算法侵入框架。

“后端可替换”目前主要是构建时替换 `COMPUTE_BACKEND_SOURCES`，不能直接理解为稳定的动态插件 ABI 或无修改 GPU 移植。局部最小二乘仍使用 Eigen，CPU 场计算也不是完全独立的设备执行抽象。公共安装目录同时包含普通 solver API 和维护接口；建议清楚区分两类支持承诺。

## 3. 数值算法逐项判断

### 稳态 SIMPLE

已追踪动量预测、Rhie–Chow 面通量构造、压力修正、压力松弛、速度/面通量修正、湍流更新及全局终止判据。压力采用物理压力，黏度是动力黏度，phi 是体积通量，密度进入动量对流和时间项，量纲基本一致。压力修正方程与通量修正的符号配对正确。

这里的 `rAU` 按框架定义取当前缩放行表示下的 `V/aP`；原始物理源已另外乘以松弛系数，因此它不直接等于对原始物理源的响应。不能机械套用另一套松弛矩阵的 `alpha*V/aP` 公式判错。公开接口和文档对该约定有说明。建议补充对速度松弛系数变化的独立收敛测试，验证该约定下预测、修正与停止判据的一致性。

全 Neumann 压力使用全局参考单元约束，装配方式在兼容的常数零空间问题上能固定参考并保留对称结构。没有发现常规层流入口/壁面/压力出口用例中系统性的符号或密度错误。

结论：**当前受测层流稳态 SIMPLE 基本成立**，有二维/三维、非正交、串行/MPI回归及实际收敛基准支持。该判断不覆盖 F1–F3 的 RANS、不覆盖 F5/F6 所有特殊组合，也不等于复杂流动均已验证。

### 瞬态 SIMPLE 与时间离散

时间历史在步边界推进，单步内反复求解不会错误推进 old/older；Euler 和固定步长 BDF2 系数、BDF2 首步 Euler 启动路径有测试。Euler 可处理缩短的末步；BDF2 对固定步长及时间区间有约束，不应声称支持任意变步长二阶格式。

新增完整 transient SIMPLE 的单元离散扩散 ODE 检查：三个 `dt=0.1` 时间步，对照已离散空间系统的 Euler/BDF2 递推，最大误差分别 `8.94e-12`、`6.76e-12`。它检查时间系数和求解状态机，**不等于连续物理问题的时间收敛阶验证**。

瞬态算法属于每步内重复 SIMPLE，没有独立 PISO/PIMPLE 校正结构，也没有独立的瞬态面通量历史修正。仅凭这点不能直接定性错误；仍需要 Taylor–Green、振荡通道等包含非平凡压力耦合的基准，进行 dt 缩小、内迭代收敛和面通量一致性检查。

结论：**时间离散和基本步进有依据，完整非定常流动精度尚未充分验收**。当前证据不足以宣称在所有 dt 下获得预期压力/速度时间阶。

### 湍流模型

| 模型 | 代码对照结果 | 可接受范围与不足 |
|---|---|---|
| SA | 正变量标准 SA 结构；检查了 fv1/fv2、ft2、cw1、cb2 梯度项等 | 不是 SA-neg；需要真实 wallDistance，近壁派生量受 F2 影响 |
| kOmega | 明确是 Wilcox1988m，常数与模型命名基本一致 | 不是 SST，也不是 Wilcox2006；不能混用其边界和验证期望 |
| kEpsilon | 标准高 Re k-epsilon 的常数、生产和耗散结构基本对应 | 当前无完整自动壁函数策略；显式破坏项依赖松弛并受 F3 影响 |

模型公式对照：[SA 定义](https://tmbwg.github.io/turbmodels/spalart.html)、[Wilcox 版本定义](https://tmbwg.github.io/turbmodels/wilcox.html)、[标准 k-epsilon 文档](https://doc.openfoam.com/2212/tools/processing/models/turbulence/ras/linear-evm/rtm/kEpsilon/)。

已执行稳态/Euler × none/SA/kOmega/kEpsilon × 1/2 进程共 16 个接线 smoke run，全部退出成功，说明字典选择、SIMPLE 调用、模型字段更新及输出链路确实存在。它们不构成湍流物理基准。当前仓库没有对应模型的成套壁面摩阻、速度剖面、网格收敛和参考解验收。

结论：**RANS 不是空壳，但不能判为正确完成**。F1–F3 和近壁工具缺口必须先解决，之后再用已收敛的参考算例验收。残差下降或 smoke run 成功均不能替代此过程。

### 算子与装配

内部面采用等大反向贡献，扩散正交项与非正交修正分开；固定值/固定梯度普通边界的主要符号检查通过。LeastSquares 使用局部加权 3×3 系统；Green–Gauss 有面位置修正；可选 Orthogonal、Corrected、LimitedCorrected 扩散，Upwind、LinearUpwind、Central 对流。

这些是实际实现的离散，非占位接口。但 LinearUpwind 的延迟修正和 Central 不自动保证有界性，LimitedCorrected 限制非正交扩散修正也不等于完整的对流 TVD。高 Peclet、强畸变、间断系数的稳定性和阶数应单独验证。当前标量系数面插值不能直接推广为已验证的间断材料界面调和平均或各向异性张量扩散。

时间项带系数时要分清 `c*d(phi)/dt` 与 `d(c*phi)/dt`；不能据当前恒密度/热容用法就声称支持一般可压缩守恒时间项。

结论：**基础离散结构合理且覆盖较广；边界组合 F6、并行组合 F5 是明确缺陷，高阶/强畸变/间断介质还缺独立验收。**

### 线性代数、预条件器和 AMG

CG/PCG、BiCGSTAB 有真实的稀疏矩阵迭代、预条件应用和最终真残差检查；不能因为某个结果收敛就忽略 F4 的成功判据缺陷。CG 的应用范围仍是满足要求的对称正定系统，不能用于一般对流矩阵。

串行 IC 和 ILUT 依赖 Eigen 实现，ILUT 当前使用固定 drop/fill 参数。MPI 下 IC/ILUT 是对本 rank 内部块的局部因子分解，相当于块式预条件，不是全局 IC/ILUT；这属于合法实现选择，但通常随分区增加而变弱。

| AMG 路径 | 实际算法 | 评价 |
|---|---|---|
| 串行 | 非平滑聚合；分片常数 P，R=Pᵀ，Ac=PᵀAP；加权 Jacobi 平滑；递归 V-cycle；粗层 SparseLU | 是真正的轻量代数多重网格预条件器；对系数跳变、各向异性、非对称系统的鲁棒性证据不足 |
| MPI | 图连接驱动聚合，可累计多次粗化映射；细层平滑加复制式粗解校正；各 rank 建立并求解相同粗系统 | 更准确地说是分布式两层聚合预条件器；聚合阶段数不等于完整分布式多层 V-cycle |

MPI 聚合主要依赖拓扑连接，未采用与串行相同的系数强度信息。粗系统限制到 2048 行，存在粗矩阵稠密缓冲和复制求解；全局图信息也有每 rank 复制成本。它适合受控规模，但没有证明网格无关收敛、大核数弱扩展或复杂系数鲁棒性。AMG 在当前设计中是 Krylov 的预条件器，不是独立主求解算法。

结论：**预条件器和 AMG 确有实现，基本数学结构成立；“AMG 已实现”不能扩展成“成熟可扩展 AMG 已实现”。** 应报告网格尺寸、rank、迭代数、达到同一容差的总时间和峰值内存；固定迭代吞吐不能冒充收敛求解性能。

### MPI

owned/ghost 分区、全局单元编号、处理器面数据、归约、分布式 SpMV 和结果合并均有真实实现和多进程测试。常规层流、热传导、输运结果一致性较好。应用层能管理 MPI 初始化/释放，并对 rank 异常采取整体终止，降低局部失败后其他 rank 永久等待的风险。

但有三项必须分开看：

1. **正确性**：F5 说明不能宣称所有离散方法组合都具有串行/MPI等价性。
2. **内存规模**：启动时广播完整网格再分区，AMG也复制全局图/粗系统；启动和建层不满足完整的每 rank O(N/P) 存储目标。
3. **性能实现**：`KrylovHalo::begin()` 实际调用阻塞的 `exchangeFirstLayer/Alltoallv`，`finish()` 主要改变状态。接口名称和内外矩阵拆分并不代表已经实现通信与计算重叠。Alltoallv 还具有随 rank 数增长的元数据成本。

结论：**MPI 主链路有效，有中小规模证据；正确性存在特定缺口，大规模扩展能力还不能背书。**

## 4. 第二轮全代码检查得到的工程评价

第二轮不仅重读 physics/algebra，也追踪了 mesh_reader、mesh 几何构造、配置解析、字段生命周期、应用错误路径、写出与后处理；审查文件清单和哈希保存在证据目录。它发现了 F6/F7，并核实下列边界。

| 范围 | 当前优点 | 需要改进 |
|---|---|---|
| 网格 | 显式拓扑、几何缓存、非法输入检查、乱序非结构回归 | 六面体限制需醒目；缺分布式读入和分区缓存；更复杂畸变需质量验收 |
| 字段 | 持有稳定引用、区分 Cell/Face、隐藏内部存储 | 值操作不含BC、halo有效性等合同应更明确；缺模型局部命名空间 |
| 方程接口 | 表达式接近方程，外部 solver 已证明可用 | 缺通用隐式源项、张量算子、块耦合和多字段组合能力 |
| 方法/求解配置 | 离散格式可按字段覆盖，普通字典拒绝未知项 | 线性求解主要按 scalar/vector 配置；字段覆盖名缺使用核验 |
| 时间/运行 | old/older管理和错误传播较完整 | 单活动 RunTime 限制多区域耦合；没有通用事务式步失败回退 |
| 输入/初值 | case 文件划分直观，模型参数集中 | 内部初值主要 uniform，patch值常量；缺一般非均匀场文件输入 |
| 输出 | 分 rank 写出、全局 ID、有限性和分片完整性检查 | F7；缺 restart、网格/配置/版本来源、完整落盘标志 |
| 构建/测试 | 架构、外部API、负例和多进程测试都实际存在 | 未见仓库内CI配置；需数值参考构建、编译器矩阵和持续物理验收 |

尤其是 SA：它需要空间分布的真实 wallDistance，但当前内置字段输入仅支持 uniform internal。开发者可以用 C++ 构造非均匀场，普通用户却不能仅靠一般 case 文件导入真实复杂几何距离。再加上没有内置壁距生成、壁函数与 y+ 工作流，“有模型选择开关”距离“用户能可靠设置湍流算例”仍有明显差距。

当前结果文件保存单元场，不包含完整的面通量、old/older等求解历史，因此应叫结果输出，不能当作可恢复计算的 checkpoint。默认本机优化、LTO、fast-math 配置适合性能实验；建议同时保留参考数值构建，持续比较误差与成功状态。直接使用含 Eigen 类型的维护接口还要匹配其编译/对齐配置，普通 solver API 应继续屏蔽这些要求。

## 5. 三种角色的开发和使用体验

### 框架维护者：可以维护，但正确性合同落后于模块划分

代码体量紧凑，主依赖方向清楚，有架构测试和外部构建测试。维护者能追踪从方程到矩阵再到结果的完整链路，不需要先理解巨型对象体系。这一点很好。

主要成本来自：大型 operators.cpp 中重复的 scalar/vector/BC 路径；串行/MPI求解实现的判据漂移；字段值、边界、halo、临时量之间的隐含约定；输入输出结果来源与完成状态不完整。F2/F4/F5/F6 都体现了“模块各自合理，组合合同不完整”。

**评价：日常维护难度中等，数值核心维护仍依赖作者经验。** 优先把共同合同变成少量明确接口和跨路径验收，而不是继续增加抽象类。可按通量、梯度、边界求值、隐式装配拆分内部实现，保留稳定公共 API。

### 物理 solver 开发者：支持的方程家族内比较容易，超出后会迅速遇到能力边界

写热传导、被动标量、多个弱耦合对流扩散方程时，Case 管字段、eqn 表方程、RunTime 管历史、后端自动并行，这是一条真正可用的开发路径。外部求解器无需包含 MPI/Eigen/internal，是实际测试结论。

新增类似现有 RANS 的模型，可以实现私有 Model、读取 physics 参数、创建变量并注册工厂；SIMPLE 不需要知道每个模型常数或输运细节。这种边界应保留。稳态和瞬态有各自状态机也合理，但应提取共用的私有动量/压力通量计算内核，避免修一次漏另一份。

限制也很具体：刚性反应缺 `Sp/SuSp`；完整应力/各向异性扩散缺通用张量算子；压力和湍流变量共享 scalar 求解配置；多个模型使用全局字段名容易冲突；多字段代数往往需要临时场；多区域和强耦合系统没有现成执行合同。

**评价：标量输运类新物理可以敏捷开发；新的数学结构还不能仅靠写一个 physics 目录完成。** 这是可明确改进的能力边界，不应通过让 solver 作者访问私有矩阵或 MPI 来绕过。

### 求解器用户：层流教程可用，独立搭建可信复杂算例还不够容易

case 组织、字典未知键拒绝、命令行求解和 ParaView 输出降低了基本使用门槛。已有 cavity/poiseuille 等实例，使用户能从可验证问题开始。

复杂问题的准备成本仍高：网格格式和单元类型有限；一般非均匀初场和 wallDistance 难输入；没有 restart；近壁模型策略需用户自行补齐；线性容差、外迭代变化量和真正方程残差不够容易区分；F4/F7 会直接影响对“成功”和历史结果的信任。

**评价：具备 CFD 背景的研究用户能使用受验证案例；面向普通工程使用者还缺完整算例准备、校验和结果追溯流程。** 用户体验的优先改进不是增加更多模型名字，而是让现有模型有可靠默认配置、明确适用范围和可重复验收结果。

## 6. 敏捷扩展能力：按实际任务区分

| 新任务 | 当前需要的工作 | 敏捷程度 |
|---|---|---|
| 新热方程、被动标量、非刚性显式源项 | physics 内字段、表达式、参数与循环 | 较高 |
| 多个松耦合标量 | 多字段注册与循环、明确历史和收敛 | 中等至较高 |
| 新涡黏性输运模型 | 私有模型及工厂；还需修复共同应力/BC/残差问题 | 接口较容易，可信验收成本高 |
| 刚性化学反应/强源项 | 先补通用隐式源项及非线性收敛机制 | 当前偏低 |
| 变黏度完整应力、各向异性材料 | 先补张量离散与边界合同 | 当前偏低 |
| 可压缩、固体力学、强块耦合、多区域 | 守恒变量、通量、方程类型和执行能力扩展 | 当前不属于现成能力 |

不能把“通用”理解为任何 PDE 都无需扩展框架；合理目标是让通用数学能力一次加入后，可被多个 physics 复用。每次能力增加都应同时给出串行/MPI一致性、边界行为和外部 solver 示例。

## 7. 改进顺序与验收标准

**第一阶段：修复会影响正确性的共同路径。**

- 修 F4：跨尺度线性系统准确返回状态，串行/MPI采用一致容差语义。
- 修 F5：声明或调度组合算子模板；非线性扭曲网格在 1/2/4 ranks 与串行对照。
- 修 F1/F2/F3：完整应力、派生量壁面规则、后更新残差与裁剪诊断；再验收各模型。
- 修 F6/F7：入流/出流组合边界和历史网格来源；不得静默生成错误成功结果。

**第二阶段：补齐能被多个物理模块复用的数学与配置能力。**

- 通用线性隐式源项、张量通量/散度和变量系数制造解。
- 按方程/字段选择线性 solver、preconditioner、容差及必要参数；检查覆盖名是否真正使用。
- 统一边界求值和派生字段边界构造；模型局部字段命名及多字段计算工具。
- 提取稳态/瞬态共享私有数值内核；保留清楚的独立状态机。

**第三阶段：建立物理可信度和用户完整工作流。**

- 瞬态速度/压力 dt 阶验证，内迭代收敛与连续性检查。
- SA/kOmega/kEpsilon 分别提供对应版本、近壁策略、壁距和收敛的参考算例；报告摩阻/剖面及网格敏感性。
- 支持非均匀场输入、壁距生成、restart/checkpoint、运行来源与完整结果标志。
- 自动化执行架构、API负例、MPI方法组合和有明确误差阈值的物理验证；保存参考构建与性能构建结果。

**第四阶段：在前述正确性稳定后推进并行规模。**

- 分区后分发/直接读入本地网格，减少每 rank 全量几何和图复制。
- 真正非阻塞 halo 与内部 SpMV 重叠；用测量区分通信、归约、装配、预条件和写出成本。
- 使粗层分布与层数对应真实多层算法，或接入具有清晰合同的成熟分布式代数后端；保持 solver API 不变。
- 用相同收敛目标测总求解时间、迭代增长和每 rank 内存，进行强/弱扩展评估。

## 8. 本次实际完成的验证与证据边界

| 执行 | 结果 |
|---|---|
| `make -j4 test test-workflow test-external` | 退出 0；串行、架构、工作流、外部 solver、负例通过 |
| `make -j4 test-mpi` | 退出 0；1/2/4 ranks 数学/求解链路通过，含乱序非结构网格 |
| `make -j2 validate-cavity validate-poiseuille` | 退出 0；两项均先达到求解器终止条件，再计算误差 |
| 16 个模型接线 smoke run | 退出 0；只作为模型选择/更新/输出接线证据 |
| 独立尺度、湍流、halo、边界、后处理诊断 | 暴露 F2–F7；F1 为解析反例 |
| 完整 transient SIMPLE 离散 ODE | Euler/BDF2 与离散递推一致；不替代物理时间阶验证 |

主要量化结果：

- Re=100，64×64 cavity，MPI 4：2908 次外迭代达到终止条件，质量指标 `1.350446e-14`，`dU=3.150412e-7`，`dP=9.970449e-7`；Ghia 对比 `u_Linf=0.0041522014`、`v_Linf=0.0085114402`，均满足该基准 0.01 阈值。
- Poiseuille：973 次外迭代收敛；`u_max=1.49884024`，出口剖面 `Linf=0.00085364265`，`L2=0.0004627372`。
- 原有 parallel_math：22 个操作 × 3 种扩散方法，1/2/4 ranks 最大差分别为 0、0、`2.84217e-14`；这没有覆盖 F5 新构造的 Green–Gauss 组合。
- heat 串行/两进程输出最大绝对差 `4.1218644e-12`。

这里的层流误差都来自已达到本算例终止条件的计算，没有把固定次数中途结果当作最终精度。日志中相对线性残差大于相对阈值不必然失败，因为绝对阈值也可以满足；F4 的反例则同时违反了用户规定的容差，二者应分清。

证据入口：[README 与复现说明](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/README.md)、[测试日志](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/serial.log)、[MPI 日志](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/mpi.log)、[物理基准日志](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/validation.log)、[复现脚本](/home/midway/BabelSim/docs/reports/framework-review-2026-09-09-evidence/reproduce.py)。

本次产物仅为报告和审查证据。列出的缺陷尚未修复；“已完成 review”不表示“框架已全部正确”。
