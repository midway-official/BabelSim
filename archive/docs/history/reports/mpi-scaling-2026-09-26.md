# BabelSim 当前 CPU/PETSc 后端的 MPI 扩展性研究

> 历史开发记录：保留原报告的接口、命令、数据与结论，只用于追溯对应时期的工作。
> 不作为当前实现或性能的使用说明；现行文档见 [文档中心](../../../../docs/README.md)，归档范围见 [历史索引](../README.md)。

测试日期：2026-09-26。基线提交：`fa7d7c3d31f578f452182aa946344a06d249239c`。本报告评估当前 CPU 版 PETSc 求解器在单机 MPI 下的强扩展，并加入 SIMPLE、PISO 流动算例。原始逐次/逐 rank 性能记录、输入字典、网格生成参数和汇总表保存在 [证据目录](mpi-scaling-2026-09-26-evidence)；当前机器和既有 AMD ROCm 环境评估见 [ROCm 可行性报告](rocm-feasibility-2026-09-26.md)。

## 结论

当前框架在 10 万至 100 万单元上能从 MPI 获得实际求解加速，但效率取决于工作量和算法：约 10 万单元的热方程在 4 rank 最快，继续加到 8 rank 后变慢；50 万至 100 万单元在 4 rank 分别达到约 2.4× 和 3.0× 加速。SIMPLE 和 PISO 在 4 rank 达到约 2.1× 和 2.5×，但其 Krylov 迭代数会因分区改变而上升，且 SIMPLE 的 8 次外迭代没有收敛。

503 万单元网格已生成，但在当前 WSL 约 11 GiB 可见内存下，求解器尚未进入有计时的求解阶段就占用约 10.8 GiB RSS；为避免 VM 耗尽内存，测试被手动停止。因此，**本研究没有 500 万单元的求解时间或 MPI 效率**。这说明当前分布式网格读入/广播/分区路径是进入超大规模求解前的现实门槛。

## 测试对象与计时方法

测试运行于 WSL2（Linux 6.6.87.2）、AMD Ryzen AI MAX+ 392 / Radeon 8060S 主机可见的 CPU 环境：12 个物理核、24 个逻辑线程，Open MPI 4.1.2，CPU PETSc 3.25.5 + Hypre。所有 MPI 测试固定 `OMP_NUM_THREADS=1`、`OPENBLAS_NUM_THREADS=1`、`MKL_NUM_THREADS=1`，rank 绑定物理核，使用同一个 BabelSim 二进制（SHA-256 见证据 JSON）。没有使用 ROCm/GPU；这些数据是后续 GPU 对照的 CPU 基线。

生成网格均为面显式存储的三角柱网格：交错对角线的二维三角网格沿第三维挤出，主要规模为 101,306、500,094、1,024,000 和 5,030,912 个单元。非规则版本对内部网格点施加最大为平面网格间距 0.15 倍的平滑扰动；101,306 单元另有未扰动对照。它们是多面/面邻接形式存储的非结构网格，但都属于同一三角柱拓扑族，并非四面体、混合单元或随机连通性的全面覆盖。热方程使用生成的一步线性求解；SIMPLE/PISO 在同一几何族上使用入口、出口、壁面和对称面边界。网格生成器和流动算例生成器分别见 [generate_mpi_scaling_case.py](../../../../tools/generate_mpi_scaling_case.py) 和 [generate_mpi_flow_case.py](../../../../tools/generate_mpi_flow_case.py)。

**求解时间定义：**每次运行取所有 rank 中最大的 `solverComputeSeconds` 为 `T_p`；它等于 `max(0, solverSeconds - local.outputSeconds)`。计时从求解器入口开始，因此排除了 Case 构造期间的网格/字段读取、MPI 初始化以及写出时间；求解器内部的方程装配、PETSc/Hypre 预条件器建立和线性求解仍计入。该定义对应代码中的计时边界和写出扣除逻辑（[application.cpp](../../../../src/runtime/application.cpp#L234)、[application.cpp](../../../../src/runtime/application.cpp#L91)）。强扩展加速比 `S_p = T_1/T_p`，并行效率 `E_p = S_p/p`。SIMPLE/PISO 虽固定了外迭代数或时间步数，但 Krylov 迭代数仍随分区变化，因此它们的效率同时反映了并行开销和数值迭代工作量变化，不能解读为纯通信效率。

每个 rank 数值先做一次预热（1m 热方程例外），再重复 1–3 次；下表 `n` 是计入的重复数。50 万和 100 万单元的大多数结果只有一次正式测量，适合看趋势，不足以刻画运行波动。100k warped 热方程、101k SIMPLE/PISO 有 2–3 次重复，结论相对可靠。具体逐次数据见 `run-attempts.csv`，逐 rank 计数器见 `rank-metrics.csv`，每组未经汇总的 harness 元数据见 `raw-metadata/`。

## 求解时间与 MPI 效率

### 热方程强扩展

| 网格 / 单元数 | MPI ranks | 求解时间 `T_p`（秒） | 加速比 | MPI 效率 | 重复数 |
|---|---:|---:|---:|---:|---:|
| 扭曲三角柱 / 101,306 | 1 | 0.9612 | 1.00× | 100% | 3 |
|  | 2 | 0.6245 | 1.54× | 77.0% | 3 |
|  | 4 | 0.4448 | 2.16× | 54.0% | 3 |
|  | 8 | 0.4735 | 2.03× | 25.4% | 3 |
| 未扭曲对照 / 101,306 | 1 | 0.9450 | 1.00× | 100% | 2 |
|  | 4 | 0.4125 | 2.29× | 57.3% | 2 |
| 扭曲三角柱 / 500,094 | 1 | 6.9574 | 1.00× | 100% | 1 |
|  | 2 | 4.4774 | 1.55× | 77.7% | 1 |
|  | 4 | 2.8703 | 2.42× | 60.6% | 1 |
| 扭曲三角柱 / 1,024,000 | 1 | 19.2025 | 1.00× | 100% | 1 |
|  | 2 | 10.4775 | 1.83× | 91.6% | 1 |
|  | 4 | 6.5063 | 2.95× | 73.8% | 1 |

### SIMPLE 与 PISO

流动边界条件为入口定速、出口定压、上下壁面无滑移、前后面 symmetry，密度为 1、动力黏度为 0.01。SIMPLE 固定最多 8 次外迭代；PISO 执行一个 `dt=0.01` 的瞬态时间步和 2 次 corrector。两者的每 rank 设置和生成输入见 `case-inputs/`。

| 算法 / 网格 | MPI ranks | 求解时间 `T_p`（秒） | 加速比 | MPI 效率 | 重复数 | 求解状态 |
|---|---:|---:|---:|---:|---:|---|
| SIMPLE，101,306 单元 | 1 | 6.4776 | 1.00× | 100% | 2 | 达到 8 次上限，未收敛 |
|  | 2 | 4.7901 | 1.35× | 67.6% | 2 | 达到 8 次上限，未收敛 |
|  | 4 | 3.0509 | 2.12× | 53.1% | 2 | 达到 8 次上限，未收敛 |
| PISO，101,306 单元，1 步 | 1 | 1.4946 | 1.00× | 100% | 2 | converged |
|  | 2 | 0.9836 | 1.52× | 76.0% | 2 | converged |
|  | 4 | 0.6052 | 2.47× | 61.7% | 2 | converged |
| PISO，500,094 单元，1 步 | 1 | 11.6395 | 1.00× | 100% | 1 | converged |
|  | 2 | 7.8003 | 1.49× | 74.6% | 1 | converged |
|  | 4 | 4.5695 | 2.55× | 63.7% | 1 | converged |

SIMPLE 的时间是固定 8 次外迭代的吞吐测试，不应解读为收敛解的耗时或不同 MPI 数下达到同一物理解的成本。PISO 的 `converged` 表示求解器在指定单步内满足其停止条件，不代表这是独立验证过的基准流动物理问题。

## 主要时间开支

计数器将 `preconditionerSeconds` 记为 `KSPSetUp` 和 BJacobi 子求解器 setup 的耗时；`linearSolveSeconds` 是 `KSPSolve` 的耗时。它们在代码中先 setup、后求解，计时阶段顺序不重叠（[petsc_linear_system.cpp](../../../../src/backend/petsc_linear_system.cpp#L226)、[petsc_linear_system.cpp](../../../../src/backend/petsc_linear_system.cpp#L304)）。`KSPSolve` 内含矩阵向量乘、预条件器应用、Krylov 更新与 PETSc/Hypre 内部的归约通信。JSON 的 `preconditionerApplySeconds` 为 0 是因为后端没有单独为 PCApply 计时，不代表求解时没有应用预条件器。

| 算例与 ranks | 总求解时间 | KSPSetUp | KSPSolve | 其他求解工作 |
|---|---:|---:|---:|---:|
| 热方程 101k，1 rank | 0.961 s | 0.505 s（53%） | 0.428 s（45%） | 约 0.028 s（3%） |
| 热方程 101k，4 ranks | 0.445 s | 约 0.267 s（60%） | 约 0.158 s（35%） | 约 0.020 s（5%） |
| 热方程 101k，8 ranks | 0.473 s | 约 0.337 s（71%） | 约 0.117 s（25%） | 约 0.020 s（4%） |
| 热方程 1.024m，1 rank | 19.203 s | 10.205 s（53%） | 8.604 s（45%） | 约 0.394 s（2%） |
| SIMPLE 101k，1 rank | 6.478 s | 约 3.376 s（52%） | 约 2.619 s（40%） | 约 0.483 s（8%） |
| PISO 101k，1 rank | 1.495 s | 约 0.487 s（33%） | 约 0.755 s（51%） | 约 0.252 s（17%） |

因此，这批求解工作主要花在预条件器 setup 和 Krylov 解算；热方程约 95–98% 在这两段。PISO 的场更新、动量/压力校正外层计算占比相对大一些。这里的“其他”是总求解计时扣除两个阶段后的余量，包含字段/方程运算、矩阵更新、RHS 和真残差检查等；PETSc 后端 `assemblySeconds` 只计装配计划生成矩阵/RHS 值的部分，不覆盖所有 FVM 离散运算，不能把它当作完整方程装配成本。

## 性能差异的原因分析

**小热算例从 4 rank 加到 8 rank 不再划算。**101k 扭曲网格在 4 rank 的 0.4448 秒已是本组最低；8 rank 用时升至 0.4735 秒，慢约 6.4%。热算例在各 rank 下均为 1 个线性系统、4 次 Krylov 迭代，所以回退不是迭代次数增多造成的。计数器显示，8 rank 的每 rank 线性求解阶段从 4 rank 的约 0.158 秒降至约 0.117 秒，但预条件器阶段由约 0.267 秒升至约 0.337 秒；总计时仍上升。rank 变多还增加分区邻接（rank 0 邻居数从 3 增至 5），每个 rank 只处理约 1.3 万个单元，局部计算更难摊薄进程调度、预条件器和 MPI/PETSc/Hypre 协调成本。这支持“小问题已开始被固定开销主导”的解释。报告中的 `globalReductions` 是 BabelSim 自身记录的规约次数，不包含 PETSc/Hypre 内部计数；当前计时器也不能把 PETSc/Hypre 内部通信单独分离出来，因此这是基于分阶段计时与分区数据的原因判断，不是直接测出的通信占比。

**问题变大后，更容易摊薄并行固定成本。**1.024m 热算例的 2 rank 效率达到 91.6%，4 rank 仍有 73.8%；500k 热算例 4 rank 效率约 60.6%。一个线性系统在这些热测试中均保持 4 次 Krylov 迭代，说明主要收益来自较大的局部矩阵/向量工作量被分到更多 CPU 核，而算法迭代负担没有变重。测试仅覆盖到 4 rank 的中大型网格，不能推断 8/12 rank 下的最佳点。

**流动算法的迭代工作量会随分区变化。**101k SIMPLE 在 1/2/4 rank 下记录的 Krylov 迭代数从 122 增至 165、180；PISO 从 20 增至 22、23，500k PISO 从 22 增至 25、25。SIMPLE 各 rank 数量均记录 32 次线性求解调用，PISO 单步均为 5 次；外迭代/时间步工作数固定，变化的是线性迭代工作量。应用自身的规约计数在 SIMPLE 从 1 rank 的 104 次增至多 rank 的 402 次，在 PISO 从 25 次增至 104 次，也体现并行路径增加了全局协调操作。两者线性求解器配置含 BJacobi 和 Hypre；域分解后局部问题与预条件器作用改变，可能导致 Krylov 迭代数上升，进而削弱强扩展效率。SIMPLE 4 rank 的 53.1% 尤其应结合迭代增量和未收敛状态解读。规约计数只包括 BabelSim 自身调用，不是 PETSc/Hypre 内部通信的完整统计。

**读入/分区时间不在表内，但会限制实际可用规模。**按用户要求，`caseSetupSeconds` 未加进上述求解时间。它在 101k 热算例约为 1 rank 0.70 秒、2 rank 约 1.5 秒、4 rank 约 1.3 秒、8 rank 约 1.7 秒；1.024m 算例则约为 8.6、18.2、16.9 秒。多 rank 时的 Case 构造明显更慢，尽管这些时间不影响表中 MPI 效率。代码路径由 rank 0 读取完整网格、广播全局网格，再由各 rank 分区（[parallel_context.cpp](../../../../src/parallel/parallel_context.cpp#L513)）。数据提示网格读入与分发本身也需要优化，才能让端到端任务获得与纯求解计时相同的收益。

**网格扭曲影响小于算法/分区因素，但样本有限。**相同 101,306 单元下，未扭曲与扭曲几何的 1 rank 时间分别是 0.9450 与 0.9612 秒，4 rank 为 0.4125 与 0.4448 秒。差异不足以说明严重网格几何惩罚；但对照只有 2 次重复，且仍是同一网格拓扑，因此不能据此概括其他单元质量或连通模式的影响。

**500 万单元失败的直接观察是容量边界，不是数值失败。**5,030,912 单元输入网格文件为约 780 MB；单 rank 构造阶段运行约 107 秒时，进程树峰值 RSS 约 10.78 GiB，接近此 WSL 环境可见的约 11 GiB RAM。没有输出 solver performance JSON，不能计算加速比；该进程被手动停止以免耗尽 VM 内存。因为这是单 rank 尝试，代码走 `readDistributedMesh` 的单进程分支并直接读取网格（[parallel_context.cpp](../../../../src/parallel/parallel_context.cpp#L520)），没有执行 MPI 网格广播；所以这次峰值反映单进程网格读取/对象构建的容量压力，不能归因于多 rank 广播或重复分区。运行没越过初始化阶段，也不能断言各项内存开销的精确来源。详细记录见 [5m 尝试说明](mpi-scaling-2026-09-26-evidence/five-million-attempt.txt)。

## 优化顺序与结论边界

建议按“先看清 KSP 内部，再调 PC/KSP，最后优化网格分发”的顺序做下一轮：

| 优先级 | 优化动作 | 针对的观测与验证要求 |
|---|---|---|
| 1. 补齐诊断 | 在受控诊断构建中启用 PETSc event log，按 rank 拆出 `KSPSetUp`、`KSPSolve`、MatMult、PCApply；另用 MPI profile/trace 和 Hypre 诊断补充 collective 与 Hypre 内部阶段。 | 当前记录说明 setup 与 KSP 总耗时最大，但 PCApply 与 PETSc/Hypre 通信仍混在 `KSPSolve` 中。应用使用 `PetscInitializeNoArguments()` 并显式设置 KSP/PC，没有调用 `KSPSetFromOptions()`；不能假定直接传 `-log_view` 或 `-ksp_type` 会生效。应通过受控代码开关在 PETSc 初始化前启用 profiling，避免测试设置绕过 BabelSim 的容差配置。PETSc 的 [PetscLogDefaultBegin](https://petsc.org/release/manualpages/Log/PetscLogDefaultBegin/) / [PetscLogView](https://petsc.org/release/manualpages/Log/PetscLogView/) 支持按事件采集 profile。 |
| 2. 降低 KSP 同步成本 | 热方程用 BCGS；流动的动量方程用 BCGS、压力校正用 CG。若 profile 证实归约延迟显著，测试 PETSc `pipebcgs`、`pipecg`，或 CG 的 single-reduction 选项。 | `pipebcgs` 可重叠非阻塞规约与矩阵向量乘/预条件器应用，但只允许右预条件，MPI 还需要有异步进展；CG 只适用于矩阵与 PC 满足 SPD 条件。PETSc 文档说明了 [PIPEBCGS 的同步/预条件约束](https://petsc.org/release/manualpages/KSP/KSPPIPEBCGS/) 和 [CG single-reduction / pipelined 变体](https://petsc.org/release/manualpages/KSP/KSPCG/)。本项目的配置目前限定为 CG、BCGS、GMRES、FGMRES，且不解析 PETSc 命令行选项，需增加显式、可验证的配置支持；同容差检查真残差、物理守恒和总 `T_p`。热方程每次只有 4 次 Krylov 迭代，pipeline 收益预计有限；优先在迭代数较多的 SIMPLE 上 profile 后试验。 |
| 3. 调低预条件器成本 | 对 Hypre BoomerAMG 评估 coarsening、插值、平滑和层次参数；对动量方程比较当前 BJacobi/局部 ILU 与带 overlap 的 ASM 或 AMG/GAMG。 | 热方程 KSPSetUp 占总求解约一半；SIMPLE 8 外迭代记录 16 次 PC setup。Hypre 文档指出 aggressive coarsening 可降低层次复杂度和内存，非 Galerkin 粗网格稀疏可减少通信，但属实验性能力；strong threshold 依问题而异，扩散问题默认值表现较好不代表压力/动量矩阵最佳。建议逐项调参，记录 setup、PC apply、Krylov 次数、AMG 层数/算子复杂度和总时间，参考 [BoomerAMG 调优说明](https://hypre.readthedocs.io/en/stable/solvers-boomeramg.html)。PETSc 会在矩阵值改变时自动重建 PC；[KSPSetReusePreconditioner](https://petsc.org/main/manualpages/KSP/KSPSetReusePreconditioner/) 可以强制复用旧 PC，即使矩阵值变化，但迭代数可能显著增加。因此只把它作为 SIMPLE 相邻外迭代之间的受控实验：按 `setup + KSPSolve`、外迭代收敛和真残差判断收益，不要无条件滞后更新。当前精确相同矩阵的安全复用路径应保留。 |
| 4. 匹配 rank 数与分区 | 先对 100k 热方程选 4 rank；对 500k–1m 比较 2、4、8、12 rank 的绝对时间和效率。再评估以切割面/halo 和 PC 成本加权的 k-way 图分区。 | 100k 热算例 8 rank 已比 4 rank 慢；1m 的 2 rank 效率高、4 rank 绝对耗时更低。当前分区能平衡单元数，但不保证面数、ghost 或 AMG 工作均衡；[MPI 审查](mpi-audit-2026-09-26.md)也记录了这项边界。ParMETIS 的 [PartKway 算法说明](https://karypis.github.io/glaros/files/sw/parmetis/manual.pdf)以 edge-cut 为目标；对本求解器这只是候选指标，还需测实际每 rank 方程装配、PC setup 和 solve 时间，验证切边减少是否转化为求解加速。不要因缩短三层 halo 就牺牲目前复合离散所需的正确性。所有大网格点至少做 3 次正式重复。 |
| 5. 单独改进端到端容量 | 改为只在 root 读全局网格并分区一次，再将局部网格发给各 rank；评估并行读入/分区文件。若目标是 5m 单 rank，还需降低全量 mesh 对象的构造峰值内存。 | 1m 的 case setup 已达约 8.6–18.2 秒，且不计入本报告求解时间；5m 单 rank 在初始化阶段就触及 WSL 内存上限。此项可能改善启动时间与可计算规模，却不会直接改变当前口径下的 `T_p`，应将端到端时间单独报告。 |

对当前单机 CPU 跑法，约 100k 单元热方程先用 4 rank；8 rank 没有性能收益。500k–1m 单元的热方程和 100k–500k 的 PISO 在已测范围内用 4 rank 有明显收益。SIMPLE 并行收益偏低，后续应增加外迭代至收敛，并比较达到相同残差/守恒误差的总成本。大于 1m 单元先改进分布式网格读入与分区路径，再扩展到 8/12 rank 和重复测试。

结论只适用于当前 CPU PETSc/Hypre、当前单机 WSL 环境、当前单一三角柱非结构网格族和已测试 rank 数。较完整的统计还应增加四面体/多面体或不同连通性的网格、多个网格质量水平、500k/1m 每点至少 3 次正式重复，以及达到同一数值收敛标准的 SIMPLE/PISO。此处 CPU 数据不能预测 AMD ROCm/HIP 的速度、GPU-aware MPI 收益或 GPU 内存容量；应把它作为之前 [ROCm 迁移可行性研究](rocm-feasibility-2026-09-26.md)所建议后续 GPU 对照的 CPU 基线。
