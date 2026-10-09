# Re=3900 圆柱绕流：复用 Re=100 几何与网格方法并局部加密

该算例沿用 `../cylinder_re100` 的圆柱、计算域、边界分区、环带加外部
Gmsh 三角网格的生成方法，以及单层三角柱挤出形式。原算例与原网格不改动；
新网格独立生成，原网格 SHA256 记录在 `case_parameters.json`。

这里计算二维 URANS：瞬态平均流动加 Wilcox1988m k–omega 湍流闭合。
展向只有一层，前后边界为 symmetry。它不能解析三维湍流尾流的小尺度、
展向失稳或转捩过程。Re=3900 的真实圆柱尾流具有三维结构；本算例不作为
LES/DNS 或三维实验结果的定量复现。

## 物理与模型

| 参数 | 值 |
| --- | --- |
| 圆柱直径 D、来流 U、密度 rho | 1、1、1 |
| Re_D | 3900 |
| 运动 / 动力黏度 | 1/3900 = 2.5641025641e-4 |
| 计算域 | x=[-10,25]D，y=[-10,10]D |
| 挤出厚度 | 0.1D，一层 |
| 入口湍流强度、长度尺度 | I=1%，l=0.07D |
| k_inlet=1.5(I U)^2 | 1.5e-4 |
| omega_inlet=sqrt(k)/(betaStar^0.25 l) | 0.3194382825 |
| 入口 mut/mu | 1.83134 |

选择现有 `turbulenceModel kOmega`，因为它支持 k=0 的近壁积分配置，
且无需本仓库尚未提供的 k–epsilon 壁面函数。代码实际模型是 Wilcox1988m，
常数 betaStar=0.09、beta=0.075、gamma=5/9、sigmaK=sigmaOmega=0.5。
本算例未添加 SST 混合函数、交叉扩散或剪切应力限制器。
模型依据为 [Turbulence Modeling Resource 的 Wilcox 定义](https://tmbwg.github.io/turbmodels/wilcox.html)。
Wilcox 模型对自由来流 omega 有敏感性；这里给出的 I、l 是算例假设，后续
定量比较需检查入口参数敏感性。

圆柱 U=0、k=0，p 零梯度。omega 采用 TMR 说明的有限距离壁面近似
omega_wall=60 nu/(beta d1^2)，使用实际壁面法向首单元中心距离
d1=0.001333288712D，得到 omega_wall=115392.3387。
该壁面近似见 [TMR 的壁面条件说明](https://tmbwg.github.io/turbmodels/sst.html)；
使用这条边界条件不会将输运模型改成 SST。
入口 U=(1,0,0)，出口 p=0 且 U/k/omega 零梯度；上下远场和前后边界为 symmetry。
初始 U/p 沿用 Re=100 的无散流函数启动方法和一次性局部扰动；
初始 k 向壁面衰减，omega 包含近壁渐近分布。所有初场按全局单元 ID 生成。

Re=3900 尾流及展向分辨率的参考为
[Kravchenko & Moin, 2000](https://doi.org/10.1063/1.870318) 和
[2022 年展向长度与网格研究](https://doi.org/10.1371/journal.pone.0266065)。
这些文献用于说明流动与维度要求，不表示本算例已经与其吻合。

## 网格与数值设置

| 项目 | Re=100 原网格 | 本算例 |
| --- | --- | --- |
| 单元数 | 20,486 | 49,350 |
| 圆周分段 | 160 | 384 |
| 环带层数 | 12 | 28 |
| 第一层厚度 | 0.015D | 0.004D |
| 层增长率 | 1.10 | 1.07 |
| 尾流目标尺寸 | 0.11D | 0.08D |

生成端质量检查全部通过：最小体积 1.63617e-6，最小三角角 23.626 度，
最大内部非正交角 38.800 度，最大归一化偏斜 0.497982；共形环带接口完整。
原生读取器独立确认 49,350 个三角柱、0 个六面体及正体积，最大面积闭合误差
约 9.99e-17。质量数据见 `mesh/mesh_quality.json` 和 `mesh/reader_quality.json`。

原正式计算 `formal_dt0p01_tight` 已结束：完成 902 个物理时间步至 t=9.02，
在 t=9.03 的压力线性求解中退出码为 2；最后完整快照为 t=9.00。
显式真实残差 `1.000002882754994e-14` 略高于绝对门槛 `1e-14`。
这一线性验收失败本身不能证明流场发散。原输入、结果和日志保留；
本次通用软件及求解设置优化已通过选定配置的数值验收，见
[本次优化报告](../../docs/reports/cylinder-optimization-2026-10-09.md)。
新正式运行标签为 `formal_dt0p01_optimized_v3`，从 t=0 开始，8 rank、BDF2、
dt=0.01、终点 t=120。运行状态以新标签的状态文件为准。

原正式计算从 t=0 的独立初场开始，
dt=0.01，目标 tU/D=120，共 12,000 个物理时间步。
时间推进为隐式 BDF2，第一步使用 Euler 启动。
动量对流使用 linearUpwind；k/omega 对流使用 upwind，因此不能宣称全系统
空间二阶精度。梯度为 leastSquares；动量及湍流扩散 limitedCorrected，压力
扩散 corrected。为抑制大时间步下的分裂及系数滞后误差，使用现有 PISO
的 Picard 外层迭代：每一外层重新组装动量方程，校正压力速度，然后更新
k、omega 及有效黏度；在同一物理时间层反复耦合，达到条件后提前结束。
`maxIterations=200` 是上限，不是每步固定迭代次数。

原正式输入的压力至少校正 3 次，按 `couplingTolerance=1e-7` 最多增加至 50 次；
每次压力校正执行 3 次额外非正交修正，共 4 次线性求解。
质量容差 1e-8；外层要求速度相对变化不超过 1e-7、动量相对残差不超过
1e-6、压力校正中的最大相对变化不超过 1e-4、湍流变化与残差均不超过
1e-6。压力线性求解使用 CG/Hypre，relativeTolerance=1e-11、
absoluteTolerance=1e-14。动量线性求解使用 BiCGStab/Block Jacobi，
relativeTolerance=1e-12、absoluteTolerance=1e-15，并使用 warm start。
速度和湍流松弛系数均为 1。

`turbulenceCoupling=iterated`，`turbulenceMaxIterations=200` 保留；
原输入在多外层模式下每一外层只求解一次湍流输运，通过整个外层循环
满足湍流容差。这里的 200 不表示每个外层另外执行 200 次湍流内迭代。
状态脚本分别记录完成的物理时间步和日志报告的外层迭代，不把中间未收敛
的外层记录当作时间步失败。

本次优化后的当前 `numerics/solution.bs` 显式设置：

```text
nCorrectors 3
maxCorrectors 50
couplingTolerance 1e-7
nonOrthogonalCorrections 0
initialNonOrthogonalCorrections 3
turbulenceCorrectors 8
```

一次性潜势初始化仍执行 4 次压力求解；步内每次 PISO 压力校正执行一次 corrected
压力求解。湍流输运内校正上限为 8，达到原湍流判据后提前退出。全部外层与真实
线性残差门槛保留。8 rank、100 步至 t=1 的窗口为 191.380 秒，五场/独立涡量
最大相对差 5.49e-6，力系数估计最大绝对差 5.66e-7，均通过验收。
1/2/4 rank 的 100 步候选较慢，且 k/mut 的部分相对 L2 指标超过本次 1e-5 预算，
所以当前正式配置保留 8 rank；这组对照不用于宣称各进程数的严格等价。

## 生成、运行与输出

从仓库根目录运行：

```bash
OPENBLAS_NUM_THREADS=1 python3 cases/cylinder_re3900_urans/generate_mesh.py
build-petsc/check_cylinder_re3900_mesh cases/cylinder_re3900_urans/mesh/cylinder_re3900.mesh > cases/cylinder_re3900_urans/mesh/reader_quality.json
python3 cases/cylinder_re3900_urans/run.py --label new_formal_label --ranks 8
```

生成脚本依赖 numpy、matplotlib 和 gmsh，并生成网格、U/p/k/omega 初场、
模型字典和派生参数。`check_cylinder_re3900_mesh` 用本算例的 `check_mesh.cpp`
和当前 `libbabelsim.a` 编译，与原算例使用相同的独立读取审计方法。
运行脚本拒绝覆盖结果；复现时需用尚未使用的 label。
原正式输入副本为 `results/inputs/formal_dt0p01_tight`，记录 Git HEAD、
求解器与网格 SHA256；原 worker 已结束。每个 MPI rank 使用一个线程，
默认 8 ranks。可以用 `--label` 创建独立输出，或 `--foreground` 在前台执行。
`--binary` 可指定已验证的冻结求解器；运行清单记录其路径与 SHA256，worker 启动前
校验二进制未发生变化。清单同时保留全局及 `equation.*` 命名方程设置和离散方法。
`--end-time`、`--dt`、`--write-interval` 只修改本次输入副本，便于复核短算。
求解器在输入副本内使用相对 `results` 目录，`results/<label>` 是指向该次
实际输出目录的链接。当前正式结果的链接为 `results/formal_dt0p01_optimized_v3`。

- `results/formal_dt0p01_optimized_v3_status.json`：状态、PID、已收敛物理时间、外层迭代与实测范围。
- `results/formal_dt0p01_optimized_v3_solver.log`：完整求解日志，包括实际 Courant 数、质量与耦合残差。
- `results/formal_dt0p01_optimized_v3/<time>/rank-*`：每 25 步，即每 0.25 个对流时间保存 U/p/k/omega/mut。
- `results/formal_dt0p01_optimized_v3_performance`：各 rank 的执行报告。
- `validation/formal_dt0p01_optimized_v3/analysis`：成功结束后自动生成涡量、速度、k、mut/mu 图，
  近壁 y+、尾流探针、升阻力估计、有限时窗频谱及动画。

短算实测每个快照约 23.84 MB；t=120 共 480 个快照，预计 CSV 约 11.4 GB，
加上输入、日志与图像。结束时间估计随实际速度更新，不是完成承诺。
守护程序检查完整全局 ID、有限值、正值和明显发散范围；这些检查不构成物理验证。

运行中查看最新完整快照：

```bash
OPENBLAS_NUM_THREADS=1 python3 cases/cylinder_re3900_urans/postprocess.py --latest-only
```

后处理通过全局 ID 合并任意 rank 数，独立最小二乘梯度的仿射精确性已经检查。
涡量为 omega_z=dUy/dx-dUx/dy，使用对称色标。升阻力和 y+ 使用一层切向
速度梯度估计，未做壁面牵引积分的网格收敛验证。
只有足够的后期采样才给出频谱；启动瞬态不报告为已验证的涡脱落。

## 本次验证

四进程 `validation/smoke_dt0p001` 已完成 200 步、t=0.2：退出 0，200 条
PISO 记录均 accepted=true / linear=ok / turbulenceSettled=true / converged=true；
最大 Co=0.402744，最大相对质量不平衡 3.2185e-13。4 个保存快照包含完整有限
U/p/k/omega/mut，k 与 omega 均为正，未丢失全局 ID。
首单元 y+ 一层估计最大 1.11382、最后快照 P95=1.08724，接近 1。
这是启动短算证据，长时间周期、三维湍流结构、模型适用性及网格/时间步独立性
仍需分别评估。当前正式运行状态以 `results/formal_dt0p01_optimized_v3_status.json` 为准。

八进程 `validation/pilot_mpi8` 已完成 1,000 步、t=1：退出 0，20 个完整有限
快照，全部 PISO 步接受、线性求解通过、湍流内迭代收敛。最大 Co=0.456421，
最大相对质量不平衡 4.08965e-13，最大动量耦合残差 4.43096e-7；
近壁 y+ 估计最大 1.30731，t=1 的 P95=1.24187。
该试算采用此前 dt=0.001 配置的 1 次最小压力校正和 1 次非正交修正。
在 t=0.2，与上述 4 进程/较多校正的配置比较，U/p/k/omega/mut 的相对 L2
差异最大为 4.21303e-4，详见 `validation/pilot_mpi8/comparison_at_t0p2.json`。
该对照同时改变了 rank 数和校正数，不能作为单独的 MPI 一致性或时间/空间
收敛检验。

后台启动与自动后处理通过 `results/launcher_smoke_status.json` 的完整流程
复核：10 步到 t=0.01，求解与后处理均退出 0，生成图像、CSV 和动画。
原来的 dt=0.001 长算及 dt=0.005 / 0.01 试算已经按用户要求停止，
全部输入、已保存结果和日志保留。停止时物理推进时间分别为 9.217、
3.505、5；停止与 PID 身份核验见 `validation/formal_dt0p01/predecessors_stopped.json`。
旧试算在共同时间 t=2.5 的近尾流压力面积加权差异为 17.0% 和 29.3%，
不能只据 `converged=true` 判断大时间步精度。对照见 `validation/timestep_study`。

正式启动时，首先发现单独增加外层迭代会被压力变化判据卡住：
压力求解仍用原相对容差 1e-8 时，变化量在反复迭代中停于 1e-5 到 1e-4
量级。失败日志与输入保留在 `results/formal_dt0p01*` 的对应历史目录；
判据与后续求解设置调整记录于 `validation/formal_dt0p01`。
收紧压力线性求解容差后，停滞仍存在；同时收紧动量求解容差后，
首步在 72 次外层迭代后通过，压力变化量降为 4.21365e-13、动量残差
1.78022e-16、湍流残差 4.76707e-7。最终设置使用上述 1e-4 的压力变化
判据。是否能持续到 t=120 由状态文件记录，首步收敛不构成完整长算验证。

本次新增算例和生成文件为本地文件，未提交或推送；原 Re=100 算例及其他已有工作保留。
