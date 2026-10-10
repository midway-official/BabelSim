# BabelSim

BabelSim 是用于有限体积 PDE 求解与不可压缩流动计算的 C++17 软件。它提供正式求解程序、几何后处理、文本算例配置和公开开发接口，支持串行与 MPI；用户可以运行现有模型，也可以编写新的方程和耦合算法。

本项目采用源码交付。所需 PETSc、Eigen 与 Hypre 源码统一放在 [dependencies](dependencies/README.md)，安装脚本检查系统环境并从交付目录构建依赖与应用。完整操作、配置查询及二次开发从 [文档首页](docs/README.md)进入。

## 安装与一键构建

首次准备机器，在 BabelSim 根目录执行：

```bash
./dependencies/install.sh
```

Ubuntu/Debian 缺少系统依赖时，脚本使用 apt-get 安装，系统安装阶段按需使用 sudo。系统环境已准备好时，一条命令完成构建：

```bash
make
```

产物为 `build/babelsim-solve`、`build/babelsim-post` 与 `build/libbabelsim.a`，PETSc/Hypre 安装到 `build/petsc`。Eigen 直接从本地头文件使用。必要环境包括 GCC C++17、GNU Make、MPI C/C++/Fortran、Python ≥ 3.9、pkg-config 和 BLAS/LAPACK；依赖版本、精简范围及许可证见 [依赖说明](dependencies/README.md)。

```bash
./dependencies/install.sh --check
./dependencies/install.sh --clean --jobs 4
```

第二条清理默认构建目录后重新编译。可移植构建、调试构建、可选 METIS 及已有 PETSc 用法见 [安装手册](docs/user-guide/installation.md)。默认优化针对当前 CPU，交付预编译程序时需验证目标环境。

## 第一次计算

创建完整小型热传导算例，调用求解器，再导出可在 ParaView 打开的时间序列：

```bash
python3 docs/examples/create_case.py --template heat --output work/heat
build/babelsim-solve -case work/heat -time serial
build/babelsim-post -case work/heat -time serial/all -format vtk
```

结果保存于 `work/heat/results/serial`，时间序列为 `work/heat/post/serial/series.pvd`。`-case` 选择输入目录，case 中的 solver 选择算法；求解 `-time` 是运行名。逐文件解释、成功条件和 MPI 对比见 [首次运行教程](docs/user-guide/quick-start.md)。

```bash
mpirun -np 4 build/babelsim-solve -case work/heat -time mpi4
```

## 提供的模型与工具

| 能力 | 入口 |
| --- | --- |
| 热传导、给定速度的标量对流扩散 | `heat`、`transport` |
| 不可压缩稳态与瞬态流动 | `simple`、`transientSimple`、`piso` |
| 层流及 RANS 闭合 | `none`、`SA`、`kOmega`（Wilcox1988m）、`kEpsilon` |
| 多面体网格、单元/面场和 MPI | 原生 `BABELSIM_MESH 3` 与公开 Field 接口 |
| 空间与时间离散 | 逐方程/逐项配置，Euler/BDF2 等 |
| PETSc 线性系统 | CG、BiCGSTAB、GMRES、FGMRES，Hypre/GAMG 等预条件器 |
| 几何后处理 | 多面体 VTK、时间 PVD，适用六面体的 Tecplot |
| 扩展开发 | Case、Field、math、equ、time 与求解器注册 |

当前交付依赖集针对 CPU/MPI、实数双精度。工程使用范围、边界条件与模型限制见 [产品说明](docs/user-guide/overview.md)。软件不自动解析单位或认证工程物理精度，实际应用需相应的网格、时间步和参考验证。

## 二次开发

无需修改框架后端即可编译一个只使用公开头文件的独立求解器：

```bash
make -f docs/examples/Makefile -j4
python3 docs/examples/create_case.py --template custom-heat --output work/customHeat
build/examples/custom-heat -case work/customHeat -time external
```

[开发教程](docs/developer-guide/extending-solvers.md)详细解释方程装配、时间历史、错误状态、公开 SDK 链接、内置注册和双场耦合。完整源码与自动验收入口位于 [docs/examples](docs/examples/README.md)。

## 文档与交付

| 阅读目的 | 入口 |
| --- | --- |
| 从安装到完成计算 | [用户手册：12 个连续章节](docs/user-guide/README.md) |
| 查参数、语法、默认值及格式 | [配置参考](docs/reference/README.md) |
| 新增求解器或维护框架 | [开发指南](docs/developer-guide/README.md) |
| 运行可复制代码 | [示例目录](docs/examples/README.md) |
| 源码交付、重建和验收 | [交付说明](docs/user-guide/delivery.md) |

`docs/` 只保留现行交付文档。历史开发、审计与版本化实验证据单独存放在 `archive/docs/`，不作为当前参数或产品性能承诺。

```bash
python3 tools/check_docs.py --json build/docs-check.json
python3 docs/examples/verify_examples.py --output build/tutorial-validation
make test-architecture
```

需要完整框架验收时，按 [测试矩阵](docs/developer-guide/verification.md)执行 `make test-full` 与额外 `make test-piso-coupling`；物理基准另有方腔和通道入口。构建、回归、数值推进完成与物理验证分别记录。

BabelSim 采用 [MIT License](LICENSE)，依赖遵循各自附带许可。源码来源与分发文件见 [dependencies](dependencies/README.md)。
