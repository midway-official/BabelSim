# BabelSim

BabelSim 是面向 CFD、PDE 和多物理场计算的 C++17 有限体积框架。Physics 用公开的
Case、Field、`math`、`equ` 和时间历史接口写方程及算法循环；网格分区、MPI 通信和
PETSc 稀疏系统由内部执行层负责。

## 当前实现

- 显式 face-based 三维 polyhedral 网格，可变顶点数的面、可变面数的单元和物理 patch；
  几何量在构造时缓存，原生 reader 使用 `BABELSIM_MESH 3`。
- cell/face 标量、向量和张量场，显式边界条件与整场数学。
- 梯度、插值、守恒面通量、对流、扩散和 Euler/BDF2 时间项；支持非正交与偏斜修正。
- 过程式方程装配：具名方程绑定未知场和数值配置，逐项组装后显式求解。
- 导热 `heat`、标量输运 `transport`、稳态 `simple`、瞬态 `transientSimple` 和 `piso`；
  SA、k–ω、k–ε 通过 RANS 模块接入流动方程。
- 长期存在的 PETSc Mat/Vec/KSP/PC：分离模式、矩阵更新、RHS 更新和重复求解。
  KSP 支持 CG/BiCGSTAB/GMRES/FGMRES，PC 支持 Hypre、GAMG、Block Jacobi 等。
- MPI owned/ghost 网格与 halo；rank 0 读取及划分全局网格，工作进程保留局部子网格。
  可选 METIS 约束单元数与面工作量，无 METIS 时使用确定性图生长。
- 每个 rank 写 owned 结果，后处理按 global ID 合并为 polyhedron VTK；
  八顶点六面体结果另可写出 Tecplot FEBRICK。

Eigen 用于局部最小二乘。Physics 不依赖 MPI/PETSc 头、稀疏布局或 Field 原始存储；
分层门禁由 `make test-architecture` 检查。

## 构建与运行

依赖 GNU Make、GCC 对应的 MPI C++ 编译器与 `gcc-ar`、Eigen 3、启用 MPI 的 PETSc。
压力使用 `hypre` 时，PETSc 需包含 Hypre。构建前设置实际安装前缀：

```bash
export PETSC_DIR=/path/to/petsc-install
make -j4
build-petsc/babelsim-solve -case cases/heat -time demo
build-petsc/babelsim-post -case cases/heat -time demo/all -format vtk
mpirun -np 4 build-petsc/babelsim-solve -case cases/poiseuille -time mpi4
build-petsc/babelsim-post -case cases/poiseuille -time mpi4 -format vtk tecplot
```

默认 `make` 产出 `build-petsc/libbabelsim.a`、`babelsim-solve`、`babelsim-post`；测试显式
执行。编译可用 `-j`，求解任务和性能测点一次执行一个。默认 `-march=native` 和快速
数学适合本机测量，不保证跨 CPU 可执行或逐位可复现。
依赖、可选 METIS、可移植构建、ASan/UBSan、外部链接见 [构建说明](docs/building.md)。

`-time demo` 命名整次运行：时间序列在 `results/demo/<物理时间>/`，最近写出的快照在
`results/demo/`。未命名运行的物理时间目录直接位于结果根目录，最近快照由
`output.bs` 的 `timeName` 命名。不同进程数或输入使用不同运行名。
文件、初场和重启见 [网格、场与结果文件](docs/mesh-and-results.md)。

## 测试与性能证据

```bash
make test-architecture
make test-full
make test-piso-coupling
```

命令逐项执行。`test-full` 顺序运行单元、工作流、外部 Solver、RANS、标量时间阶、
SIMPLE 一致性及 MPI 组；PISO 耦合专用组另行执行。物理基准入口为
`make validate-cavity` 和 `make validate-poiseuille`，范围见 [验证说明](docs/validation.md)。

最近工程优化覆盖 1,024,000 单元 PISO 三棱柱、512,000 单元 PISO Voronoi 和
500,094 单元 SIMPLE，比较 1/4/8/12 MPI 进程。时间排除初始化、首次完整迭代和输出；
强扩展使用各版本自己的串行中位数，内存使用全生命周期进程树峰值。完整测点和限制见
[2026-10-07—08 优化报告](docs/reports/flow-optimization-2026-10-07.md)。
划分先经 Python 合法网格图模型比较，再接入通用 MPI 层，见 [划分模型](docs/performance/partition-model.md)。

## 文档与扩展

从 [文档中心](docs/README.md) 进入现行手册：

- [DSL 与运行时](docs/dsl-runtime-manual.md)：Case/Field/方程/时间生命周期与示例。
- [数值配置](docs/numerical-configuration.md)：方程 → 项 → 显式算子的层级和 PETSc 配置。
- [内置求解器](docs/solvers.md)：方程、控制参数、接受判据和验证范围。
- [架构](docs/architecture.md)：维护边界及后端生命周期。
- [性能测量](docs/performance/README.md)：持续计算、强扩展、计数器和内存口径。
- [历史开发记录](docs/history/README.md)：过时计划、旧审查、旧接口与旧性能报告，仅用于追溯。

外部 Solver 使用普通函数和 `SolverRegistration`，main 调用 `runApplication`。
内置扩展放入 `src/physics/<name>/main.cpp`，Makefile 自动收集。
完整示例见 [coupled_scalar.cpp](tests/examples/coupled_scalar.cpp)。
后端替换属于内部维护，不需要改变 Physics 数学抽象。

MIT License，见 [LICENSE](LICENSE)。
