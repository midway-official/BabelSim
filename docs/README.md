# BabelSim 文档中心

> 面向 CFD、PDE 与多物理场的紧凑 C++17 有限体积框架：显式组装方程，整场数学，构建期可替换计算后端。

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](../LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](../README.md#构建与运行)
[![MPI](https://img.shields.io/badge/MPI-3.x-orange.svg)](../README.md#构建与运行)
[![Eigen](https://img.shields.io/badge/Eigen-3.x-8a2be2.svg)](../README.md#构建与运行)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](../README.md#构建与运行)
[![tests](https://img.shields.io/badge/tests-make%20test-brightgreen.svg)](validation.md)

BabelSim 用显式连接的三维非结构六面体网格求解不可压流动、对流扩散与多物理场问题。框架把
网格、场、离散算子、方程、运行时和计算后端切成可独立替换的层：Physics 只用公开的
Case / Field / `math` / `equ` 接口写方程和算法循环，不接触 MPI、halo、CSR/LDU、Eigen 或
Field 底层存储；这一边界由 `make test-architecture` 与 `make test-external` 自动守门。

本目录是文档总入口；项目简介、构建细节与完整命令见 [根 README](../README.md)。

**目录**：[特性](#特性) · [快速开始](#快速开始) · [内置求解器](#内置求解器) ·
[文档地图](#文档地图) · [目录结构](#目录结构) · [测试与验证](#测试与验证) ·
[参与贡献](#参与贡献) · [文档约定](#文档约定) · [许可证](#许可证)

## 特性

- **网格与几何**：三维非结构六面体网格、边界 patch 与 cell/face/vertex 拓扑；体积、逆体积、
  中心、面积向量、单位法向、正交系数、非正交修正向量、偏斜量与插值权重在建网格时预计算。
- **场与边界**：连续存储的 scalar/vector/tensor Field（cell/face），通用边界条件
  （`fixedValue`、`fixedGradient`、`zeroGradient`、`inletOutlet`、`symmetry`、计算边界迹）。
- **离散算子**：gradient、interpolation、flux、divergence、convection、diffusion、laplacian、
  time derivative；支持 Least-Squares 梯度、修正 Green–Gauss、修正面插值、非正交扩散与偏斜重构。
- **方程层**：`equ::` 显式组装——离散项立即写入绑定未知量的方程，`equ::solve` 只求解已组装
  系统；时间历史由 `time::History` 显式推进，没有延迟表达式或第二套生命周期。
- **线性代数**：LDU 稀疏装配，串行与分布式 CG / BiCGSTAB，IC / ILUT / AMG 预条件器；
  MPI Krylov 按 global cell ID 做稀疏 halo matvec 与融合归约，AMG 粗层支持跨 rank 图聚合。
- **并行**：`readDistributedMesh` 由 rank 0 解析原生 `.mesh`，并行层按单元邻接图构造每个 rank
  的 owned+ghost 局部网格与 processor patch；每个 rank 只写出 owned cell 的结果。
- **求解器**：稳态/瞬态不可压 SIMPLE、瞬态 PISO，Rhie–Chow 与压力修正留在各自求解器内部；
  RANS 湍流（Wilcox 1988 k-ω、标准 k-ε、Spalart–Allmaras）经模型接口与动量方程耦合。
- **输入输出**：原生 case / mesh / field 文件，通用并行结果写出，独立后处理生成
  ParaView VTK XML（含 `series.pvd` 时间序列）与 Tecplot FEBRICK。
- **可替换后端**：整组替换 `makeComputeBackend()` 工厂即可换掉默认 Eigen/MPI 装配与求解实现，
  Physics 与离散源码不变；不承诺动态插件或稳定后端 ABI。

## 快速开始

依赖：C++17 编译器、Eigen 3、MPI-3 实现与 GNU Make（默认面向 GCC 工具链：`mpic++` + `gcc-ar`）。

```bash
git clone git@github.com:midway-official/BabelSim.git
cd BabelSim
make -j4                                   # 构建 build/libbabelsim.a、babelsim-solve、babelsim-post

build/babelsim-solve -case cases/cavity                      # 串行稳态方腔
mpirun -np 4 build/babelsim-solve -case cases/poiseuille     # MPI 通道流
mpirun -np 2 build/babelsim-solve -case cases/heat           # 瞬态热传导（同一份 Solver 源码）
build/babelsim-post -case cases/heat -time all -format vtk   # 后处理：每个时刻的 VTK + series.pvd
```

默认 `make` 只构建、不跑测试；测试与验证必须显式调用 `make test*` / `make validate*`。
`-time <名称>` 为每次运行命名，结果落在 `cases/<name>/results/<名称>/`，可直接用
[`tools/compare_parallel_results.py`](../tools/compare_parallel_results.py) 比较不同进程数。

## 内置求解器

| 注册名 | 组织方式 | 说明 |
|---|---|---|
| `heat` | Equation-driven | 瞬态热传导；系数可以是常数或 cell 场 |
| `transport` | Equation-driven | 标量对流-扩散 |
| `simple` | Algorithm-driven | 稳态不可压 SIMPLE（`cases/cavity`、`cases/poiseuille`） |
| `transientSimple` | Algorithm-driven | 瞬态 SIMPLE（`cases/naca0012`，配合 k-ω RANS） |
| `piso` | Algorithm-driven | 瞬态 PISO |

每个求解器的控制方程、全部配置键（含默认值）、场与边界要求、收敛与失败语义，以及验证到
什么程度、证据在哪，见 [内置求解器手册](solvers.md)。RANS 不做成独立求解器：由动量方程
求解器读取 `physics` 字典里的 `turbulenceModel` 启用，线性配置统一为 `solution.bs` 中的
`scalarSolver` / `vectorSolver`。

## 文档地图

| 想做什么 | 读哪份 |
|---|---|
| 写一个新求解器 / 查 DSL 算子、配置键、文件格式 | [DSL 与运行时用户手册](dsl-runtime-manual.md) |
| 用内置求解器算题 / 查它的配置键与验证状态 | [内置求解器手册](solvers.md) |
| 改框架分层、公开接口、离散或后端（维护者） | [架构与维护边界](architecture.md) |
| 改动后跑哪些验证、验收标准是什么 | [验证与维护检查](validation.md) |
| 性能测量工具、构建开关与 JSON 字段 | [性能工具说明](performance/README.md) |
| 历史验证/性能证据（按时间归档，不代表当前接口） | [reports/](reports/) |

- **[dsl-runtime-manual.md](dsl-runtime-manual.md)** — 唯一的使用者手册。含最小可运行 Solver、
  Case 与全部配置键（`case.bs`、`physics`/`methods`/`solution`/`control`/`output`）、场文件与
  网格文件格式、Field/geometry/math/equ 全部算子与语义、时间与历史、线性求解契约、诊断与监视、
  结果与后处理、并行边界、内置求解器与 RANS、开发检查清单、旧接口迁移表。
- **[solvers.md](solvers.md)** — 内置求解器手册：`heat`/`transport`/`simple`/`transientSimple`/
  `piso` 与 RANS 模块的方程、算法、配置键、场与边界、从构建到看图（`babelsim-post`）的完整
  运行流程、验证证据与验证边界。
- **[architecture.md](architecture.md)** — 面向框架维护者：层次与依赖禁令、所有权与生命周期、
  DSL 契约、Solver 独立性、配置/报告边界、维护流程与验收命令。
- **[validation.md](validation.md)** — 验证入口与最低验收线：串行/并行测试目标、工作流与外部
  Solver 验证、MPI 一致性、架构门禁、新 Solver 的验证最低线。
- **[performance/README.md](performance/README.md)** — `-performance` 输出、构建开关
  （`ASYNC_HALO`、`CSR_SPMV` 等）、benchmark 驱动与 JSON 字段。
- **[performance/backend-audit-2026-09.md](performance/backend-audit-2026-09.md)** — 计算后端
  审计记录（历史，含 A/B 数据）。

历史证据存放在 `reports/`：案例验证（Ghia 方腔、Poiseuille）、框架审查（F1–F7）与实现映射、
RANS 方程核对、MPI 一致性、后端性能优化，以及已移除的 GMRES 后端基准。它们保留原始数字与
结论用于追溯；接口、配置键和性能数字**不保证**与当前版本一致，实现细节以用户手册和源码为准。

## 目录结构

```text
include/babelsim/   公开头：mesh / field / geometry / math / equ / case / runtime / solver
src/core            Mesh 存储、拓扑与几何缓存
src/io              case / mesh / field 读取，结果写出与后处理，监视器
src/discretization  离散算子、FVM 执行层、显式组装方程、Field 数学
src/backend         Eigen/MPI 后端、稀疏装配、AMG
src/algebra         线性求解器与分布式 Krylov
src/parallel        MPI 上下文、halo 计划、并行结果写出
src/runtime         RunTime、Solver API 桥接、Application 分派
src/physics         内置求解器：heat / transport / simple / transient_simple / piso / RANS
cases/              自包含算例：cavity、poiseuille、heat、transport、naca0012
tests/              串行、并行、工作流、外部与架构门禁
tools/              网格生成、并行结果比较、性能汇总、解析解验证脚本
docs/               本目录：手册、架构、验证与归档报告
```

## 测试与验证

```bash
make test                 # 几何、算子、case/field IO、Heat/标量输运/SIMPLE、通用输出
make test-architecture    # 头依赖与分层门禁（Physics 不得越界依赖）
make test-workflow        # 新 Solver 单函数开发、双场耦合、时间序列、ParaView 读取
make test-external        # 仓库外 Solver 构建、1/2/4 进程、负向 API 编译、无 MPI 结果读取器
make test-mpi             # MPI 网格、halo、算子、线性求解、SIMPLE 与标量输运
make validate-cavity      # Re=100、二阶迎风的 Ghia 方腔快速回归
make validate-poiseuille  # 收敛的 Poiseuille 解析解比较
```

方腔验证覆盖 Re=100…10000、网格与格式灵敏度，结果与 Ghia 等基准数据逐点对照：

![Ghia 方腔中心线对比](reports/images/cavity/cavity-ghia-centrelines.png)

完整数据、算例与复现步骤见 [Ghia 方腔验证报告](reports/cavity-ghia-validation.md)。

## 参与贡献

新增求解器不需要改框架：自己的一个 C++ 源文件、一行注册、通用 main 调用 `runApplication`，
只链接公开头与预编译库。

```cpp
#include "babelsim/application.h"   // SolverRegistration / runApplication / SolverResult
#include "babelsim/case.h"          // Case
#include "babelsim/equ.h"           // equ::Equation 与离散项
using namespace babelsim;

SolverResult solveMyCase(Case& problem) {
    ScalarField& U = problem.scalarField("U");
    auto equation = equ::createEquation(U);
    equ::laplacian(equation, 1.0, -1.0);                  // 左端 -div(grad U)
    equ::source(equation, 1.0);                           // 右端体源
    if (!equ::solve(equation).converged()) return SolverResult::notConverged();
    problem.output(U);
    return SolverResult::completed();
}

const SolverRegistration registration("mySolver", solveMyCase);
int main(int argc, char** argv) { return runApplication(argc, argv); }
```

想成为内置命令，把同一份 `main.cpp` 放进 `src/physics/<name>/`，Makefile 自动收集，无需修改
启动器名单、注册宏或 Solver 基类。提交前请至少跑通 `make test-architecture` 与 `make test`，
改动框架分层或公开接口时按 [validation.md](validation.md) 补对应验证；改动离散或后端时请对照
[architecture.md](architecture.md) 的分层约束。

## 文档约定

- 文档里的语义以 `include/babelsim/*.h` 为准；文档与头文件冲突时以头文件为准。
- 每份验证/性能报告的格式：结论 → 算例与命令 → 数据表 → 复现步骤；数值结论必须可复现。
- 归档报告只增不改：反映提交当时的接口与数字，是历史证据，不是当前契约。

## 许可证

MIT License，见 [LICENSE](../LICENSE)。