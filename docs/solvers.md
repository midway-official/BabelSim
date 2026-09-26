# BabelSim 内置求解器手册

本手册逐个说明当前仓库里的内置求解器与湍流模块：**控制方程与算法**、**需要的场与边界条件**、
**全部配置键（含默认值）**、**收敛与失败语义**、**怎么跑**，以及**验证到什么程度、证据在哪**。

求解器作者想自己写一个新求解器，读 [DSL 与运行时用户手册](dsl-runtime-manual.md)（尤其是第 1、7、9
节）；本文只讲已经装好的这些求解器。验证状态用
[验证与维护检查](validation.md) 的五级词汇描述：源码存在 → 调用链接入 → 回归通过 →
短跑稳定 → 物理验证。

## 1. 状态总表

| 求解器                  | 覆盖的问题         | 最高验证级别           | 主要证据                                | 复现命令                                              |
| -------------------- | ------------- | ---------------- | ----------------------------------- | ------------------------------------------------- |
| `heat`               | 瞬态/稳态导热（标量）   | 物理验证（解析解）        | `T=2t` 逐点精确、1/2/4 rank 一致           | `make test-workflow`                              |
| `transport`          | 标量对流-扩散（速度给定） | 回归通过 + 短跑稳定      | 单胞解析值、1/2/4 rank 逐时刻一致              | `make test-workflow`、`make test-mpi`              |
| `simple`             | 稳态不可压层流       | 物理验证（基准数据 + 解析解） | Ghia Re=100/400/1000、Poiseuille 抛物线 | `make validate-cavity`、`make validate-poiseuille` |
| `transientSimple`    | 瞬态不可压（含 RANS） | 短跑稳定             | 1/2/4 rank 的 Euler/BDF2 逐步一致        | `make test-simple-parallel`、`make test-mpi`       |
| `piso`               | 瞬态不可压（算子分裂）   | 短跑稳定（单算例 6000 步） | Re=1000 层流平面射流，100,000 正交六面体；质量误差 ≤ `1.48e-14` | [`cases/planar_jet`](../cases/planar_jet/README.md)、`make test-architecture` |
| RANS（k-ω / k-ε / SA） | 湍流黏性输运        | 回归通过（方程核对 + 收敛阶） | TMR/OpenFOAM 逐项核对、Euler≈2/BDF2≈4    | `make test-rans`                                  |

RANS 是模块而不是注册名：由 `simple` / `transientSimple` / `piso` 按 `physics` 字典里的
`turbulenceModel` 启用（第 7 节）。

怎么把上表任意一个求解器跑起来，见 [第 2 节](#2-怎么跑从构建到看图)（构建、命令行、结果
目录与后处理）；每个求解器的具体命令、预期输出与改法写在各自章节的 **运行** 小节。

## 2. 怎么跑：从构建到看图

这一节是所有求解器共用的操作流程；各求解器自己的命令、预期输出与改法在对应章节。

### 2.1 构建

```bash
make -j4          # 产出 build-petsc/babelsim-solve 与 build-petsc/babelsim-post
```

依赖 C++17 编译器、PETSc MPI 构建、Eigen 3、MPI-3 实现和 GNU Make。Eigen 用于局部最小二乘梯度，
PETSc 用于稀疏线性系统。默认优化含 `-march=native`，
换机器（尤其异构集群）需按 [根 README](../README.md#构建与运行) 重新编译。

### 2.2 命令行

```text
build-petsc/babelsim-solve -case <算例目录> [-time <运行名>] [-performance <目录>]
```

| 参数 | 必填 | 含义 |
|---|---|---|
| `-case <目录>` | 是 | 算例目录，相对当前目录或绝对路径都行；`case.bs` 里的路径都相对该目录 |
| `-time <运行名>` | 否 | 给本次运行命名，决定结果写到哪个子目录（见 2.4）。名字是任意字符串（`run1`、`mpi4`），与物理时间、进程数无关 |
| `-performance <目录>` | 否 | 每个 rank 一份 `rank-XXXX.json`：阶段耗时、Krylov 迭代数、halo 字节数、分区信息（字段含义见 [performance/README.md](performance/README.md)） |

**`-time` 不是 MPI 必需的**：`mpirun -np 4 build-petsc/babelsim-solve -case cases/heat` 也能跑。
它解决的是“两次运行的结果不要混在一起”：不加 `-time` 时所有运行都往
`results/<物理时间>/` 写，不同进程数的 rank 文件落在同一目录，后处理会以
`rank directory count does not match metadata` 之类的一致性问题拒绝读取，也分不清
哪份结果来自哪次运行。串并行对比就靠它：

```bash
build-petsc/babelsim-solve -case cases/poiseuille -time serial
mpirun -np 4 build-petsc/babelsim-solve -case cases/poiseuille -time mpi4
python3 tools/compare_parallel_results.py \
  cases/poiseuille/results/serial cases/poiseuille/results/mpi4 \
  --atol 5e-6 --rtol 5e-6
```

### 2.3 算例目录

五个求解器都从 `case.bs` 读同一套结构（逐键说明见用户手册第 3 节）：

```text
cases/<名字>/
├── case.bs                  # solver <注册名> 与其余文件的相对路径
├── mesh/*.mesh              # BABELSIM_MESH 3 polyhedral 网格，patch 带角色（wall/inlet/outlet/symmetry…）
├── fields/initial/*.field   # 每个场的初值与边界条件
├── physics/*.bs             # 物性/模型常数（各求解器读哪些键见对应章节）
├── numerics/methods.bs      # 空间/时间离散格式
├── numerics/solution.bs     # 线性求解与算法数值控制
├── control.bs               # startTime / endTime / deltaT
├── output.bs                # directory / timeName / writeInterval / 字段筛选
└── results/                 # 运行生成，不纳入 Git
```

跑新问题时复制最接近的内置算例再改文件即可，不需要重新编译：

```bash
cp -r cases/heat /tmp/my-heat
# 编辑 /tmp/my-heat 下的 case.bs、control.bs、physics/*.bs …
build-petsc/babelsim-solve -case /tmp/my-heat
```

初始内部场默认可写 `internal uniform (...)`。也可用 `internal file <路径>` 读取按全局单元编号排列的场值：普通文本每行是 `globalCellId` 后接该单元的标量、三个速度分量或九个张量分量。还可直接给求解器的结果 CSV 文件，或给某个结果时间目录（自动合并其中的 `rank-*/<场名>.csv` 分片）：

```text
internal file ../../results/run/5.0
```

路径相对当前 `.field` 文件；文件中每个全局单元必须恰好出现一次，因此同一网格可从不同 MPI rank 数的结果场重启。将 `control.bs` 的 `startTime` 设为该快照时间，并给 `output.bs` 设一个新的 `timeName`，避免与先前结果混写。该功能恢复场值，不保存时间离散历史；BDF2 重启后的第一步自动用 Euler 启动，下一步起恢复 BDF2。

三条硬约束（启动时校验，违反直接报错，不静默回退）：

- `solver` 必须是注册名之一（`heat`/`transport`/`simple`/`transientSimple`/`piso`），
  否则报 `unknown BabelSim solver`；
- 每个实际使用的方程必须在 `equation.<name>.*` 中提供完整线性配置（`make test-workflow` 覆盖）；
- `physics/` 与 `solution.bs` 里的键必须全部被本次运行的求解器消费，拼错或留了别的
  求解器才用的键（例如把 `simple` 算例原样改成 `piso`，留下 `pressureRelaxation`）
  会报 `unused or unknown entry`。

### 2.4 结果写到哪

结果根目录由 `output.bs` 的 `directory` 决定（相对算例目录，约定 `results`）；
`-time` 决定根目录下的子目录（最终目录名用 `-time` 的值，省略时用 `output.bs`
的 `timeName`）：

| 命令 | 时间序列 | 最终结果 |
|---|---|---|
| `build-petsc/babelsim-solve -case cases/heat` | `cases/heat/results/0.01/`、`0.02/` … | `cases/heat/results/final/` |
| `build-petsc/babelsim-solve -case cases/heat -time mpi4` | `cases/heat/results/mpi4/0.01/` … | `cases/heat/results/mpi4/` |

每个时刻目录内按 rank 分目录，每个 rank 只写自己拥有的单元：

```text
results/mpi4/0.01/
├── rank-0000/{T.csv, metadata.bs, mesh.geometry}
├── rank-0001/…
```

`*.csv` 列为 `global_id,x,y,z,value0[,value1,value2]`；`metadata.bs` 记物理时间、
rank 数与全局单元数，是后处理的一致性依据。写出规则：

- 稳态 `simple` 收敛后只写一次：最终目录 + 序列目录下一个以 `startTime`（通常 `0`）命名的
  快照（如 `cases/poiseuille/results/serial/0/`）；瞬态按 `writeInterval`（默认每步）写，
  最终时刻总是写；
- 只有收敛的时间步才落盘，失败的步不会产生结果文件；
- 默认写出“从文件加载的 cell 场”，程序自建场（`phi`、`mut` 等）要在 `output.bs` 的
  `writeFields` 白名单里列出（`excludeFields` 是黑名单，见用户手册第 3.2 节）。

### 2.5 看图：babelsim-post

```text
build-petsc/babelsim-post -case <算例目录> [-time <选择>] -format <vtk|tecplot> [...]
```

`-time` 选择要处理的结果：省略时取 `output.bs` 的 `timeName`（未命名运行的最终结果）；
`-time mpi4` 处理命名运行的最终结果；`-time mpi4/latest` 只处理最后一个物理时刻；
`-time mpi4/all` 处理整条序列并额外生成 ParaView 时间序列 `series.pvd`；
未命名运行用 `latest` / `all`。产物写在 `<case>/post/`：

```bash
build-petsc/babelsim-post -case cases/poiseuille -format vtk tecplot     # post/final.vtu、post/final.dat
build-petsc/babelsim-post -case cases/heat -time mpi4/all -format vtk    # post/mpi4/*.vtu + post/mpi4/series.pvd
```

`.vtu` 用 ParaView 打开；`.dat` 是八顶点六面体兼容的 Tecplot FEBRICK。任意面数的 polyhedral 网格使用 VTK
输出。后处理会核对每个 rank 的 metadata
与全局单元完整性并检查结果与网格一致，缺 rank、混入别的运行或网格不匹配都会报错。

### 2.6 退出码与失败语义

`0` 求解完成（所有 rank 的线性/算法判据都满足）、`2` 未收敛或数值失败、`1` 配置或运行
错误。MPI 下取所有 rank 的最坏状态，单个 rank 失败不会被其它 rank 掩盖。

## 3. `heat`：瞬态/稳态导热

```text
ρc_p ∂T/∂t = ∇·(k∇T) + Q
```

每个时间步组装一次并线性求解：`ddt(ρc_p, T)`、`laplacian(k, T, -1)`、`source(Q)`。
`methods.time = steady` 时时间项不加入，退化为“每步重解一次稳态导热”。密度与热容以乘积
`ρc_p` 进入时间项，`T` 初值来自 `fields/initial/T.field`。

**场与边界**

| 场   | 类型     | 位置   | 说明                                                   |
| --- | ------ | ---- | ---------------------------------------------------- |
| `T` | scalar | cell | 未知温度；边界可用 `fixedValue` / `zeroGradient` / `symmetry` |

**配置键**

| 文件         | 键                               | 类型          | 说明            |
| ---------- | ------------------------------- | ----------- | ------------- |
| `physics`  | `density`                       | positive    | ρ             |
| `physics`  | `heatCapacity`                  | positive    | c\_p          |
| `physics`  | `conductivity`                  | nonnegative | k，可为 0（纯瞬态蓄热） |
| `physics`  | `source`                        | number      | Q，可为负         |
| `solution` | `equation.temperature.*` | 完整配置 | 温度方程的线性求解设置 |

**收敛语义**：每步只有一次线性求解，`SolveResult` 即状态；未收敛返回退出码 2。

**运行**：内置算例是 32 单元的一维棒（两端 `fixedValue`，其余 4 个 patch 为 `symmetry`），
`control.bs` 给出 0→0.05、dt=0.01，共 5 步。

```bash
make -j4                                               # 首次运行前先构建
build-petsc/babelsim-solve -case cases/heat                  # 串行
mpirun -np 2 build-petsc/babelsim-solve -case cases/heat -time mpi2
mpirun -np 4 build-petsc/babelsim-solve -case cases/heat -time mpi4
```

控制台每个时间步一行 `heat time=<t> residual=<r>`，本算例残差在 `1e-16` 量级，
退出码 0 表示 5 步全部收敛：

```text
$ build-petsc/babelsim-solve -case cases/heat -time demo
heat time=0.01 residual=6.26041e-17
heat time=0.02 residual=1.96501e-17
…
heat time=0.05 residual=6.1994e-17
```

结果写到 `cases/heat/results/<物理时间>/rank-XXXX/`；带 `-time mpi4` 时整批写进
`cases/heat/results/mpi4/`（见 2.4），两种进程数的结果互不干扰。看图：

```bash
build-petsc/babelsim-post -case cases/heat -time demo/all -format vtk    # post/demo/series.pvd 等
```

改这个算例：`control.bs` 的 `endTime`/`deltaT` 管时长与步长；`physics/thermal.bs` 的
`conductivity`/`source` 管物性（把 `source` 改成 2、所有边界改为绝热/对称，就是验证用
的 `T=2t` 场景）；`numerics/methods.bs` 的 `time` 改成 `steady` 后时间项不再加入，
每步重解一次稳态导热。改完直接重跑，不需要重新编译。

**验证**：`make test-workflow` 把算例改成“绝热 + 均匀源 2”，此时解析解为 `T = 2t`，测试要求
1/2/4 rank 的所有单元在 0.02/0.04/0.05s 与解析解相差小于 `1e-11`；`make test` 里的
`heat_solver_test` 另用单胞 FVM 手算结果（`T=10/22`、变系数 `T=20/32`）作精确回归。

## 4. `transport`：标量对流-扩散

```text
∂(storage·C)/∂t + ∇·(phi·C) = ∇·(D∇C) + S,    phi = U·Sf（面积积分通量）
```

速度场 `U` 是**给定**的：启动时算一次面通量 `phi = math::flux(U)`，此后保持不变；只有 `C` 是未知量。

**场与边界**

| 场   | 类型     | 位置   | 说明                                                              |
| --- | ------ | ---- | --------------------------------------------------------------- |
| `C` | scalar | cell | 未知浓度；入口 `fixedValue`、出口 `zeroGradient`、侧壁 `symmetry`            |
| `U` | vector | cell | 预设速度场；入口 `fixedValue`、出口 `zeroGradient`、壁面 `fixedValue (0 0 0)` |

求解器自己创建并登记 `phi`（face 标量）作为输出场。

**配置键**

| 文件         | 键                               | 类型                              | 说明                                |
| ---------- | ------------------------------- | ------------------------------- | --------------------------------- |
| `physics`  | `storage`                       | positive                        | 时间项系数（蓄积系数）                       |
| `physics`  | `diffusivity`                   | nonnegative                     | D，可为 0                            |
| `physics`  | `source`                        | number                          | 体源 S                              |
| `solution` | `equation.transport.*`          | 完整配置                        | 输运方程的线性求解设置                    |
| `methods`  | `equation.transport.convection` | upwind / linearUpwind / central | 或使用 `equation.transport.term.convection.convection` 覆盖项 |

**收敛语义**：与 `heat` 相同——每步一次线性求解。

**运行**：内置算例是 32 单元的通道（入口 `fixedValue` 浓度、出口 `zeroGradient`），
`control.bs` 给出 0→0.05、dt=0.01，共 5 步。

```bash
build-petsc/babelsim-solve -case cases/transport                 # 串行
mpirun -np 4 build-petsc/babelsim-solve -case cases/transport -time mpi4
```

每个时间步一行 `transport time=<t> residual=<r>`（与 `heat` 同格式，本算例残差在
`1e-17` 量级）。结果目录结构同 2.4；看图：

```bash
build-petsc/babelsim-post -case cases/transport -format vtk tecplot   # 最终结果
```

改这个算例：`physics/transport.bs` 的 `storage`/`diffusivity`/`source` 是物性；
`numerics/methods.bs` 可给 `transport` 方程单独换对流格式（
`equation.transport.convection linearUpwind`），
格式只影响对流项的离散精度，求解流程不变。速度场 `U` 来自 `fields/initial/U.field`，
不参与求解——换速度只需改这个文件，不需要动求解器。

**验证**：`make test-workflow` 检查 1/2/4 rank 在同一物理时刻的 `C` 一致（阈值 `1e-8`）与时间序列
目录；`make test-mpi` 的 `parallel_transport_test` 用零通量单胞算例比对解析值 `C=0.3`（2 rank，
误差 < `1e-12`）。仓库内没有对标准对流的基准数据验证。

## 5. `simple`：稳态不可压 SIMPLE

每个外迭代做五个显式步骤（源码 `src/physics/simple/main.cpp`，无隐藏的算法封装）：

1. 组装动量方程 `div(phi,U,ρ) - div(μ_eff∇U) = -∇p`（湍流时另加显式应力余项）；
2. 测残差 → 欠松弛 `equ::relax(αU)` → `equ::scale(αU)` 保持行归一化 → 求解，得 `V/aP = rAU`；
3. Rhie–Chow：在内部面上把 `phi` 关联到压力梯度；
4. 压力修正方程 `-div(rAU ∇p') = -div(phi_H)`，先 `referenceIfUnanchored(0)`，解完后
   `p += αP·p'`、`U -= rAU·∇p'`、`phi = phi_H + flux(p')`；
5. 湍流输运（若启用）→ 用 `rU / dU / dP / mass` 判收敛。

`methods.time` 必须是 `steady`，否则直接报错。

**场与边界**

| 场   | 类型     | 位置   | 说明   |
| --- | ------ | ---- | ---- |
| `U` | vector | cell | 未知速度 |
| `p` | scalar | cell | 未知压力 |

求解器创建 `phi`（face）、`pPrime`（修正场，与 `p` 同边界但齐次）、`muEffective`
（湍流时另创建 `mut`）。推荐边界组合（`cases/cavity`、`cases/poiseuille` 即此约定）：

| Patch 角色       | `U`                 | `p`            |
| -------------- | ------------------- | -------------- |
| 入口             | `fixedValue`        | `zeroGradient` |
| 出口             | `zeroGradient`      | `fixedValue`   |
| 壁面             | `fixedValue`（静止则 0） | `zeroGradient` |
| 对称/二维退化        | `symmetry`          | `symmetry`     |
| 通用 / processor | `zeroGradient`      | `zeroGradient` |

封闭腔体（无压力出口）依赖压力修正的零空间参考值，不需要额外设置。

**配置键**

| 文件         | 键                               | 默认值  | 说明                                |
| ---------- | ------------------------------- | ---- | --------------------------------- |
| `physics`  | `density`                       | 必填   | ρ                                 |
| `physics`  | `dynamicViscosity`              | 必填   | μ；**层流也必须给**，它是 `muEffective` 的初值 |
| `physics`  | `turbulenceModel` 等             | 可选   | 见第 7 节                            |
| `solution` | `maxIterations`                 | 1000 | 外迭代上限                             |
| `solution` | `velocityRelaxation`            | 0.7  | αU                                |
| `solution` | `pressureRelaxation`            | 0.3  | αP                                |
| `solution` | `nonOrthogonalCorrections`      | 1    | 压力修正中非正交项的显式遍数                    |
| `solution` | `continuityTolerance`           | 1e-8 | 连续性（面通量守恒）相对不平衡                   |
| `solution` | `velocityTolerance`             | 1e-7 | 相邻外迭代 `‖ΔU‖`                      |
| `solution` | `momentumTolerance`             | 1e-6 | 动量方程相对残差 `rU`                     |
| `solution` | `pressureCorrectionTolerance`   | 1e-6 | `p'` 相对 `p` 的幅值                   |
| `solution` | `equation.pressureCorrection.*` / `equation.momentum.*` | 完整配置 | 压力修正与动量的独立线性配置 |

**收敛语义**：`converged = 线性全部收敛 && rU ≤ momentumTolerance && dU ≤ velocityTolerance
&& dP ≤ pressureCorrectionTolerance && mass ≤ continuityTolerance`；启用湍流时还要求
`dTurb ≤ turbulenceTolerance && rTurb ≤ turbulenceTolerance`。达到收敛才写出最终结果；
迭代耗尽返回 `notConverged`（退出码 2），不写出“看起来完成”的结果。

**运行**：`simple` 没有时间循环，一次运行给出的就是收敛后的稳态场；`control.bs` 的
`endTime` 不参与推进（`cases/cavity` 用 `endTime 0`），`methods.time` 必须是 `steady`。

```bash
build-petsc/babelsim-solve -case cases/cavity                        # 64² Re=100 方腔，串行
mpirun -np 4 build-petsc/babelsim-solve -case cases/poiseuille        # 通道流，4 rank
mpirun -np 4 build-petsc/babelsim-solve -case cases/cavity -time mpi4 # 已有串行结果时用 -time 分开存
```

控制台在首迭代、之后每 100 次、以及收敛时各报一行迭代指标：

```text
$ build-petsc/babelsim-solve -case cases/cavity
SIMPLE 1 mass=1.34188e-09 dU=1 rU=1 dP=3.33333 linP=9.9512e-09 … converged=false
SIMPLE 100 mass=1.3434e-13 dU=0.00332187 rU=0.00586256 dP=0.00877011 … converged=false
…
SIMPLE 2908 mass=1.28267e-14 dU=3.1504e-07 rU=7.39902e-07 dP=9.97043e-07 … converged=true
```

`converged=true` 后才写出 `cases/cavity/results/final/rank-XXXX/{U.csv,p.csv}`（带
`-time mpi4` 时是 `results/mpi4/`），退出码 0；迭代耗尽未收敛则退出码 2、不写结果。看图：

```bash
build-petsc/babelsim-post -case cases/cavity -format vtk tecplot     # post/final.vtu、post/final.dat
```

改这个算例：Re 由 `physics/simple.bs` 的 `density`/`dynamicViscosity` 决定（腔体边长 1、
盖速 1 时 Re = ρ·U·L/μ，`cases/cavity` 是 1/0.01 = 100）；松弛因子与收敛容差在
`numerics/solution.bs`，放宽 `momentumTolerance` 等可以更快结束；`numerics/methods.bs` 的
`equation.momentum.convection`/`equation.momentum.gradient` 影响精度与鲁棒性（一阶迎风稳但耗散大，验证数据见下表）。

**验证（本仓库验证程度最高的求解器）**

| 内容                                | 结论                                                                                 | 证据                                                                                     |
| --------------------------------- | ---------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| Ghia 方腔 Re=100/400/1000，128² 壁面加密 | 中心线两个速度分量的 Ghia 采样点 RMS < 顶盖速度 1%                                                  | [Ghia 方腔验证报告](reports/cavity-ghia-validation.md)                                       |
| 网格收敛                              | 32²→64²→128² 整线 RMS 差异逐级下降 4.1–6.2 倍                                               | 同上，`cases/cavity/validation/plot_cavity_convergence.py`                                |
| 对流格式                              | Re=1000、64²：一阶迎风 `u/v` RMS `6.20e-2/7.31e-2` → 二阶 `linearUpwind` `2.64e-3/6.51e-3` | 同上                                                                                     |
| 串并行一致                             | 1 vs 2/4 rank 最大绝对差 `3.63e-9`–`1.98e-8`（U）、`1.46e-9`–`1.37e-8`（p）                  | 同上；`make test-simple-parallel`                                                         |
| Poiseuille 解析解                    | 出口剖面与抛物线比较，L2 门限 `5e-3`                                                            | `make validate-poiseuille`                                                             |
| 快速回归                              | 2D 64²、3D 8×8×6、2D/3D 扭曲网格（非正交修正）                                                  | `make test`、[simple-parallel 一致性报告](reports/simple-parallel-consistency-2026-09-10.md) |
| Re=10000                          | `256²` 准稳态快照，只检查离散稳定性，**不作为精度验收**                                                  | [Re=10000 报告](reports/cavity-re10000.md)                                               |

报告明确限定：验证支持 **Re ≤ 1000、稳态层流、二维退化三维**；不证明更高 Re、非稳态转捩或任意网格上的精度。

## 6. `transientSimple` 与 `piso`：瞬态不可压

两者共用与 `simple` 相同的动量预测 / Rhie–Chow / 压力修正部件，区别在时间步内的组织方式：

| <br />          | `transientSimple`                   | `piso`                                                           |
| --------------- | ----------------------------------- | ---------------------------------------------------------------- |
| 时间步内结构          | 反复做带欠松弛的动量/压力迭代，直到步内收敛              | 一次动量预测 + `nCorrectors` 次压力修正（修正步不欠松弛）                            |
| `maxIterations` | 步内迭代上限（默认 1000）                     | 预测–修正流程的额外外层遍数，**默认 1 即标准 PISO**                                 |
| 时间格式            | `euler` / `bdf2`（BDF2 首步自动降为 Euler） | 同左；必须瞬态                                                          |
| 时间步接受判据         | 步内迭代达到与 `simple` 相同的收敛组合            | 要求线性求解成功、守恒及湍流输运收敛；`maxIterations > 1` 时额外要求外层收敛 |
| 动量预测欠松弛         | `velocityRelaxation`（默认 0.7）        | 同左，但修正步始终施加完整修正                                                  |
| 压力欠松弛           | `pressureRelaxation`（默认 0.3）        | 不适用（修正不做欠松弛）                                                     |

`piso` 使用的键：`maxIterations`(1)、`nCorrectors`(2)、`nonOrthogonalCorrections`(1)、
`velocityRelaxation`(0.7)、`continuityTolerance`(1e-8)、`velocityTolerance`(1e-7)、
`momentumTolerance`(1e-6)、`pressureCorrectionTolerance`(1e-6)。
启用湍流且 `maxIterations=1` 时，`turbulenceMaxIterations`（默认 1000）限制
同一时间层的模型输运内迭代；保持 `turbulenceRelaxation` 和模型接口不变，历史只在
物理步开始时推进一次。模型变化量与初始残差均须达到 `turbulenceTolerance`，
且所有湍流线性方程必须收敛。超过内迭代上限或线性迭代上限时返回
`notConverged`，不推进下一时间步。`maxIterations>1` 时仍由原有外层迭代收敛模型，
不使用此内迭代设置。守恒判据不满足时返回
`notConverged`（退出码 2），不会带着质量不平衡继续推进。

场与边界、`physics` 键与 `simple` 相同（`density`、`dynamicViscosity`，湍流时可加模型键）。
瞬态算例的 `control.bs` 必须给出正的 `deltaT`；`output.bs` 的 `writeInterval` 控制写出间隔。

**运行**：`transientSimple` 的内置算例是 `cases/naca0012`（112k 单元、transientSimple + k-ω，
0→8.0、dt=0.001，共 8000 步，完整跑一遍很慢；第一次跑先把 `endTime` 改小）：

```bash
cp -r cases/naca0012 /tmp/naca-short
# 编辑 /tmp/naca-short/control.bs：endTime 8.0 → 0.01（先跑 10 步看流程）
mpirun -np 4 build-petsc/babelsim-solve -case /tmp/naca-short -time smoke
```

PISO 的层流示例是 [`cases/planar_jet`](../cases/planar_jet/README.md)：二维平面射流，展向只有一层六面体，
前后面使用对称边界；没有启用 RANS/LES 湍流模型。喷口宽度 `D=1`、出口方向域长 `20D`，
网格为 `500×200×1=100000` 个正交六面体，`Re_D=1000`。算例采用 `deltaT=0.01 D/Uj`、
`endTime=60 D/Uj`、BDF2（首步由求解器以 Euler 启动）、动量对流 `linearUpwind`、每步三次压力修正，
并每 100 步保存一次 `U` 与 `p`。在仓库根目录运行：

```bash
python3 cases/planar_jet/generate_mesh.py
mpirun -np 2 build-petsc/babelsim-solve -case cases/planar_jet
python3 cases/planar_jet/plot_flow.py
python3 cases/planar_jet/make_vorticity_gif.py
```

生成器会重建网格、几何质量报告和初始速度场；网格与算例参数、运行记录、后处理脚本及图像均保存在该目录。
已记录的运行从 `t=0` 到 `t=60` 共 6000 步，2 个 MPI rank 正常退出；6000 条 PISO 记录均为 `linear=ok`
且求解器报告 `converged=true`，最大相对质量误差为 `1.47991e-14`（容差 `1e-8`），写出了 60 个时刻的结果。
`converged=true` 对 `maxIterations=1` 的 PISO 时间步主要表示守恒判据通过，不是物理精度判据。可查看
[运行摘要](../cases/planar_jet/validation/run_summary.json)、[网格报告](../cases/planar_jet/mesh/mesh_quality.json)
和[涡量 GIF](../cases/planar_jet/validation/planar_jet_vorticity.gif)。

同样改法对 `transientSimple` 也适用（它保留 `solution.bs` 的 `pressureRelaxation`）。
`transientSimple` 在步内迭代收敛时报
`Transient SIMPLE <iter> mass=… dU=… rU=… dP=… linear=ok converged=true`；`piso` 每步报
`PISO <遍数> mass=… dU=… rU=… dP=… linU=… linP=… converged=…`，`maxIterations 1` 时
时间步要求线性求解成功、守恒量 `mass` 达标，以及启用时的湍流输运收敛；
速度的物理步间变化无需趋零（第 6 节开头）。两者都按 `writeInterval`
写时间序列（见 2.4），看图用 `build-petsc/babelsim-post -case <算例> -time mpi4/all -format vtk`。

**复制算例改 `solver` 时的坑**：`numerics/solution.bs` 与 `physics/*.bs` 里属于原求解器的键
必须删干净（例如 `pressureRelaxation` 之于 `piso`），否则启动即报 `unused or unknown entry`；
反过来启用 RANS 时，`turbulenceModel` 要求的场（`k`/`omega` 等）必须在 `fields/initial/`
里有对应文件，缺失或名字写错会在启动时报读取 `fields/initial/<名字>.field` 失败。

**验证状态**

- `transientSimple`：**短跑稳定**。[1/2/4 rank 一致性报告](reports/simple-parallel-consistency-2026-09-10.md)
  在 2D、3D、扭曲方腔上用 Euler 与 BDF2 各跑 5 个时间步，逐时刻比较全部 U 分量与 p
  （阈值 `5e-6 + 5e-6·max|·|`），并复核每个时间步都收敛；`make test-mpi` 与
  `make test-simple-parallel` 可在当前源码上复跑。`cases/naca0012` 是可运行示例，
  **没有**对应的定量验证报告。
- `piso`：**单算例长时程数值运行完成；按本手册五级词汇最高仍记为“短跑稳定”**。
  [`cases/planar_jet`](../cases/planar_jet/README.md) 的 100,000 单元算例完成 6000 步，质量判据通过，
  但尚无针对射流基准的定量比较、网格/时间步收敛研究或专门的 PISO 回归用例。探针后处理显示波动在
  本次运行后段减小；每 1 个 `D/Uj` 保存一次不足以确认持续周期性涡脱落。因此这组结果说明该配置完成了
  数值推进，不构成物理验证，也不能据此断言存在稳定涡街。结构独立性仍由 `make test-architecture`
  覆盖；若要形成物理结论，还需按 [validation.md](validation.md) 第 5 节补充基准比较、网格与时间步研究、
  更密的探针采样和失败路径验证。

## 7. RANS 湍流模块

在 `physics/*.bs` 写 `turbulenceModel <名称>` 即启用（`simple` / `transientSimple` / `piso` 自动接入）：

| 名称（忽略大小写与下划线）                   | 模型                                       |
| ------------------------------- | ---------------------------------------- |
| `none` / `laminar` / `off`（或不写） | 层流                                       |
| `kOmega` / `wilcox1988`         | Wilcox 1988m 高 Re k-ω，无 SST 交叉扩散         |
| `kEpsilon` / `standardKepsilon` | 标准高 Re k-ε                               |
| `SA` / `spalartAllmaras`        | 标准 Spalart–Allmaras（正变量、无 trip 源、保留 ft2） |

写错模型名会报 `unsupported turbulenceModel; expected none, SA, kOmega or kEpsilon`。

**运行**：RANS 不是独立命令，它跟随所在的动量求解器——第 5、6 节的命令照旧，区别只在
`physics/*.bs` 多一行 `turbulenceModel`、`fields/initial/` 里多几个输运场。仓库里唯一内置的
RANS 算例是 `cases/naca0012`（transientSimple + k-ω，短跑做法见第 6 节）。

**需要的场**（放在 `fields/initial/`，初值与壁面渐近值都要给）

| 模型  | 输运场                      | 说明                                                 |
| --- | ------------------------ | -------------------------------------------------- |
| k-ω | `k`、`omega`              | 壁面 `omega` 需按高 Re 渐近式给大值（`cases/naca0012` 用 `1e4`） |
| k-ε | `k`、`epsilon`            | 壁面同样按渐近式给                                          |
| SA  | `nuTilda`、`wallDistance` | `wallDistance` 是几何到最近壁面的距离场，由用户提供                  |

模块自己创建 `mut`（涡黏）与 `muEffective`（`μ + μ_t`，层流分支等于 μ）两个程序自有场：启用模型时
`mut` 自动登记为输出场，`muEffective` 默认不写出（需要时在 `output.bs` 的 `writeFields` 里列出）。

**模型常数（`physics/*.bs`，全部有默认值）**

| 模型  | 键                                                                                                                                                                            |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 公共  | `density`、`dynamicViscosity`（进入 `muEffective`）                                                                                                                               |
| k-ω | `kOmegaBetaStar` 0.09、`kOmegaBeta` 0.075、`kOmegaGamma` 5/9、`kOmegaSigmaK` 0.5、`kOmegaSigmaOmega` 0.5、`kMin` 1e-12、`omegaMin` 1e-12                                           |
| k-ε | `kEpsilonCmu` 0.09、`kEpsilonC1` 1.44、`kEpsilonC2` 1.92、`kEpsilonSigmaK` 1.0、`kEpsilonSigmaEpsilon` 1.3、`kMin` 1e-12、`epsilonMin` 1e-12                                       |
| SA  | `saCb1` 0.1355、`saCb2` 0.622、`saSigma` 2/3、`saKappa` 0.41、`saCw1`、`saCw2` 0.3、`saCw3` 2.0、`saCv1` 7.1、`saCt3` 1.2、`saCt4` 0.5、`saNuTildaMin` 1e-14、`saWallDistanceMin` 1e-12 |

**耦合与收敛**：动量方程用 `muEffective` 做隐式扩散，涡黏应力的转置/无迹余项显式加入右端
（不要重复加整项）。每条输运方程在松弛前记录归一化残差、用 `turbulenceRelaxation`（默认 0.7）
欠松弛，并由 `turbulenceTolerance`（默认 1e-6）与外层收敛判据一起决定时间步/外迭代是否接受。
`k`、`omega`、`epsilon`、`nuTilda` 有下限保护，避免除零和非物理负值。

**验证**：[RANS 方程核对报告](reports/rans-equation-verification.md) 逐项对照
TMR 的 SA/Wilcox1988 定义与 OpenFOAM v7 的标准 k-ε 实现，并给出时间格式的收敛阶检验
（Euler 相邻加密误差比 ≈ 2，BDF2 ≈ 4）；`make test-rans` 在当前源码上复跑三个模型的
稳态/瞬态、1/2/4 rank 比较与下限裁剪拒绝。报告同时限定：**这不等于自动壁函数或复杂工程
壁流已经过验证**；k-ε 的原始论文未能逐式核读，采用的是标准实现交叉核对。

## 8. 怎么选

| 问题                                    | 选择                                                                  |
| ------------------------------------- | ------------------------------------------------------------------- |
| 只算温度/标量、速度已知或不需要流场                    | `heat`（纯导热）或 `transport`（对流-扩散）                                     |
| 稳态不可压层流（腔体、通道、外流）                     | `simple`                                                            |
| 瞬态不可压、时间步内需要迭代收敛                      | `transientSimple`                                                   |
| 瞬态不可压、希望每个时间步便宜且守恒（大 Co 需减小 `deltaT`） | `piso`                                                              |
| 湍流                                    | 在上面的动量求解器上加 `turbulenceModel`；高 Re 内外流优先 k-ω，工程边界层可用 SA，自由剪切层可用 k-ε |

## 9. 复现全部验证

```bash
make -j4
make test                    # 基础回归（含 heat/transport/simple 与三维/非正交腔体）
make test-architecture       # 分层门禁
make test-workflow           # heat/transport/coupled 的 1/2/4 rank 与时间序列
make test-external           # 仓库外 Solver 构建
make test-mpi                # MPI 网格、halo、算子、线性求解、SIMPLE 与标量输运
make test-simple-parallel    # 稳态/瞬态 SIMPLE 的 1/2/4 rank 一致性
make test-rans               # 三个 RANS 模型
make validate-cavity         # Re=100 Ghia 快速回归
make validate-poiseuille     # Poiseuille 解析解
```

## 10. 验证边界

- 上述“通过”只覆盖报告写明的网格、Re 数、格式与进程数；换问题类型、格式或网格后需要重新验证。
- 未收敛或迭代耗尽的结果不会被写出为成功结果；报告与结论只应引用收敛状态。
- 本框架使用 `-ffast-math`（见 [根 README](../README.md#构建与运行)），不保证逐位可复现；
  跨进程数的比较阈值按绝对/相对容差给出，不是 bit 级。
- 新求解器或新物理的最低验证线见 [validation.md](validation.md) 第 5 节。
