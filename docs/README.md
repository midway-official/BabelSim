# BabelSim 文档中心

本目录说明当前框架的使用和维护方式。2026-10-08 按公开头、Case/IO、内置 Physics、
PETSc/MPI 实现和 Makefile 核对。实验报告只证明其指定版本、输入与测量范围的结果。

## 阅读入口

| 任务 | 文档 |
| --- | --- |
| 配置依赖、构建、链接外部 Solver | [构建与运行](building.md) |
| Case、Field、`math`、`equ` 与时间历史 | [DSL 与运行时](dsl-runtime-manual.md) |
| 离散、KSP/PC、容差、预条件器复用 | [数值配置](numerical-configuration.md) |
| Heat、Transport、SIMPLE、PISO、RANS | [内置求解器](solvers.md) |
| 网格、初场、重启、并行结果与后处理 | [网格、场与结果文件](mesh-and-results.md) |
| 框架分层和后端维护 | [架构与维护边界](architecture.md) |
| 回归、解析解、MPI 和失败路径 | [验证与维护检查](validation.md) |
| 持续计算、强扩展、内存测量 | [性能测量](performance/README.md) |
| 划分与负载均衡候选比较 | [网格图模型](performance/partition-model.md) |
| 最近优化的原始证据 | [版本化实验报告](reports/README.md) |
| 旧接口、旧基准、开发过程 | [历史开发记录](history/README.md) |

## 快速运行

准备 MPI PETSc、Eigen 3、GNU Make 和 GCC 对应的 `mpic++`/`gcc-ar`。
PETSc 必须与应用使用同一 MPI 实现；`hypre` 方程还需其 Hypre 支持。

```bash
export PETSC_DIR=/path/to/petsc-install
make -j4
build-petsc/babelsim-solve -case cases/heat -time demo
build-petsc/babelsim-post -case cases/heat -time demo/all -format vtk
```

替换安装前缀；Makefile 默认值是开发机路径。默认 `make` 只构建库与应用。编译可使用
多个作业，测试组、求解任务和性能测点逐项执行，详见 [构建页](building.md)。

## 当前能力

Mesh 是三维显式 polyhedral 网格，内部面连接 owner/neighbour，face 顶点数和 cell
面数可变。薄域仍用三维网格。Field 支持 cell/face 的标量、向量和张量；`math` 立即
计算整场，`equ` 立即装配已绑定未知量的方程。时间循环、耦合和收敛判断由 Physics 表达。

生产稀疏后端是 PETSc Mat/Vec/KSP/PC，支持 CG、BiCGSTAB、GMRES、FGMRES，以及
Block Jacobi、Hypre、GAMG 等预条件器。Eigen 用于局部最小二乘。MPI 层提供局部网格、
halo 和归约；可选 METIS 提供多约束划分，无 METIS 时使用确定性图生长。

| 注册名 | 内容 | 参考 |
| --- | --- | --- |
| `heat` | 导热 | [cases/heat](../cases/heat/case.bs) |
| `transport` | 给定速度的标量对流扩散 | [cases/transport](../cases/transport/case.bs) |
| `simple` | 稳态压力速度耦合 | [方腔](../cases/cavity/case.bs)、[通道](../cases/poiseuille/case.bs) |
| `transientSimple` | 每步内迭代收敛的 SIMPLE | [源码](../src/physics/transient_simple/main.cpp) |
| `piso` | 瞬态预测与压力校正 | [平面射流](../cases/planar_jet/README.md)、[NACA0012](../cases/naca0012/README.md) |

SA、k–ω、k–ε 通过 `turbulenceModel` 接入，不是独立注册名。数值推进完成、串并行一致
和工程物理精度分别需要对应证据。

## 文档维护

API 以 [公开头](../include/babelsim/) 和实现为准。接口、配置、Make 目标变化时，更新
相应手册和示例，执行 [文档检查](../tools/check_docs.py)。收益须注明版本、相同输入、
工作量、绑定、退出状态、计时口径和重复次数。

`history/` 仅保存历史开发文件，旧配置与命令可能无法运行。保留原数字与证据，可修正
归档标记和文件链接。`reports/` 也有明确日期；冻结源码、原日志、哈希和清单保持原样。
快照中的旧文档不是自动更新的现行手册。

项目简介见 [根 README](../README.md)，许可证见 [MIT License](../LICENSE)。
