# PETSc CPU → AMD ROCm 最小迁移可行性研究

> 历史开发记录：保留原报告的接口、命令、数据与结论，只用于追溯对应时期的工作。
> 不作为当前实现或性能的使用说明；现行文档见 [文档中心](../../README.md)，归档范围见 [历史索引](../README.md)。

检查日期：2026-09-26。代码基线：`fa7d7c3`。本次完成环境检查、代码审查和迁移规划；没有安装驱动/ROCm，没有改求解器，没有运行 GPU 性能测试。

## 结论

**有条件可行，建议保留 PETSc，以 HIP 矩阵/向量和 HIP 版 Hypre 作为 GPU 实现；先验证单进程、单 GPU，再选择性迁移压力修正方程。** 当前环境尚不能直接执行该方案。架构适配成本较低，软件栈兼容和性能收益仍需实验确认。

ROCm 是计算平台，不能直接替代 PETSc 的 KSP/PC 求解能力。直接用 rocSPARSE 重写完整数值后端会额外承担 Krylov、预条件器、收敛控制、MPI 和测试维护，超出最小迁移范围。

## 1. 实测环境

| 项目 | 检查结果 | 含义 |
|---|---|---|
| CPU | AMD RYZEN AI MAX+ 392 w/ Radeon 8060S；WSL 可见 12 核/24 线程 | CPU 基线本身较强，GPU 收益不能先验保证 |
| GPU | Windows CIM：AMD Radeon(TM) 8060S Graphics；PCI DEV_1586 | 属于需要核对 Ryzen/Strix Halo 路线的平台 |
| Windows | 10.0.26200，Windows 11 | 有 WSL GPU 宿主环境 |
| 显示驱动 | 32.0.22032.16001 | 未确认对应 Adrenalin 发行版及 ROCDXG 兼容性，不能认定可用或过旧 |
| Linux | Ubuntu 22.04.5 LTS，WSL2，6.6.87.2-microsoft-standard-WSL2 | 应使用 WSL 安装路径 |
| GPU 设备 | `/dev/dxg` 存在；未发现 `/dev/kfd`、`/dev/dri` | DXCore 入口存在；WSL 中缺少 kfd 本身不证明 GPU 不受支持 |
| ROCm 工具链 | PATH 中无 hipcc/rocminfo/rocm-smi；无 `/opt/rocm*`；包清单未发现运行时 | 当前 WSL 未发现可用 ROCm 开发/运行环境；不代表所有其他 WSL 发行版也没有 |
| 内存 | 本次 WSL 可见约 11 GiB RAM，3 GiB swap | 不等于宿主总内存或 GPU 可用容量，需另用 HIP 查询 |
| PETSc | 3.25.5，real double，启用 Hypre，无 HIP/CUDA 宏 | 当前确为 CPU FP64 后端 |
| Hypre / MPI | HYPRE 3.1.0，无 HIP/GPU 宏；Open MPI 4.1.2 | GPU 版 Hypre 必须另建，MPI 的设备指针能力尚未验证 |
| 动态链接 | 求解器链接 `/home/midway/opt/petsc-3.25.5-openmpi-opt/lib/libpetsc.so.3.25` | 系统默认 pkg-config 没找到 PETSc 不代表项目没有 PETSc |

原始探测见 [environment.txt](rocm-feasibility-2026-09-26-evidence/environment.txt)。现有后端测试二进制在 1/2 rank 下退出码均为 0，见 [cpu-smoke.txt](rocm-feasibility-2026-09-26-evidence/cpu-smoke.txt)。这是当前 CPU 二进制的冒烟检查，不是重新编译验证或完整回归。

AMD 文档明确提供 ROCDXG 支持 Strix/Strix Halo 的 WSL 路线，介绍于 Adrenalin 26.2.2 + ROCm 7.2.1，覆盖 Ubuntu 22.04/24.04。[AMD WSL 指南](https://rocm.docs.amd.com/projects/radeon-ryzen/en/docs-7.2.1/docs/install/installrad/wsl/howto_wsl.html)

不能用旧版仅列 Radeon 独显的 WSL 表直接判定 8060S 不支持；也不能把家族级支持当作本机驱动、精确 SKU 和 PETSc FP64 的完整认证。最新版本选择应核对 [ROCm 7.14.1 兼容矩阵](https://rocm.docs.amd.com/en/docs-7.14.1/compatibility/compatibility-matrix.html)。预计目标架构为 gfx1151，安装后必须以 rocminfo/HIP 实测为准。

[ROCDXG 官方仓库](https://github.com/ROCm/librocdxg) 已声明源码移至 rocm-systems；安装时按目标发行版使用匹配的集成组件，不照搬旧仓库构建步骤。旧运行时可能需要 `HSA_ENABLE_DXG_DETECTION=1`，新版本应遵从对应说明。不要在 WSL 中套用原生 Linux 的内核驱动安装路线。

## 2. 代码可迁移性与最小修改位置

| 位置 | 现状 | 建议 |
|---|---|---|
| `src/backend/petsc_linear_system.cpp:73` | 固定 `MATAIJ` | GPU 模式选择 `MATAIJHIPSPARSE`，必须在 COO 预分配之前确定 |
| 同文件：89–94 | `VecCreateMPI` 固定创建 CPU Vec | 用 `MatCreateVecs` 创建相容的解/RHS，再 Duplicate；确保 residual 也为 HIP Vec |
| 同文件：78、更新矩阵路径 | 固定 COO 模式、`MatSetValuesCOO`、矩阵复用 | 保留，HIPSPARSE 支持此路径 |
| 同文件：setVectorValues / warm start | `VecGetArray` 在 host 写值 | 改用 `VecGetArrayWrite`/Restore，避免覆盖写之前无意义读回 |
| 同文件：返回 solution | 使用可写 `VecGetArray` 读回 | 改用 `VecGetArrayRead`/Restore，准确表达只读访问 |
| 同文件：setupPreconditioner | None/Jacobi/Hypre/GAMG/ICC/BJacobi/ASM | 首先验证 None/Jacobi，再 HIP Hypre；其他 PC 单独验证 |
| `src/backend/petsc_session.cpp` | `PetscInitializeNoArguments` | 保持应用参数边界；增加受控设备/诊断选项入口 |
| `src/backend/petsc_backend.cpp` | 按 equation identity 缓存系统 | 可在创建系统时按 identity 选择 CPU/HIP，设备选择在系统生命周期内固定 |
| `Makefile` | 支持 `PETSC_DIR` 与 `BUILD` 覆盖 | 独立安装 GPU PETSc，独立 build 目录，避免混用旧对象文件 |

`ComputeBackend` 接口、FVM 装配、物理模型、场 halo、网格和输出都可以在第一阶段保留 CPU 实现。现有长期保存 Mat/KSP/PC、三分量共用矩阵以及 RHS-only 复用，有利于摊薄 GPU 上传和 AMG setup 成本。

本地 PETSc 3.25.5 的 `src/mat/impls/aij/seq/seqhipsparse/aijhipsparse.hip.cxx` 已核对：存在 HIP COO 实现，host values 输入会复制至设备，HIP 矩阵默认向量类型为 VECHIP。应用可先继续提供 host COO 数值，不必立即写 HIP 装配 kernel。[PETSc COO API](https://petsc.org/release/manualpages/Mat/MatSetPreallocationCOO/)

**仅更换 PETSc 库、添加 `-mat_type aijhipsparse -vec_type hip` 不会完成本项目迁移。** 应用不向 PETSc 传命令行，且对象创建没有调用 Mat/VecSetFromOptions。建议显式、受校验的 CPU/HIP 模式；原型可用专用环境变量，正式接口再决定是否进入配置文件。不要为开放类型选择而无约束地调用 KSPSetFromOptions：这可能覆盖现有配置和收敛约定。若要设置 Hypre 的设备友好参数，使用隔离的 options/prefix 和 PCSetFromOptions，并验证生效值。

## 3. 预条件器与数值约束

| 路径 | 第一阶段决策 |
|---|---|
| None / Jacobi | GPU 正确性冒烟；不作为生产压力性能结论 |
| Hypre BoomerAMG | 压力优先目标；PETSc 和 Hypre 均编译 HIP，确认 GPU 执行 |
| Block Jacobi | 当前非对称方程默认依赖局部 ILU；CG 分支显式设 ICC。不能视为已完成 GPU 迁移 |
| ICC / ASM | 因子分解、子矩阵/子求解器及 CPU 回退需逐项确认，先保留 CPU |
| GAMG | 可作为后续对照，但粗化、平滑及粗网格求解的设备驻留需要实测 |

PETSc 官方说明 `--with-hip --download-hypre` 配合 HIP Vec 可启用 AMD GPU BoomerAMG。[PCHYPRE GPU 说明](https://petsc.org/release/manualpages/PC/PCHYPRE/)

GPU AMG 参数集合和默认值可能不同，压力当前使用 CG，因此要检查矩阵及预条件器的对称正定条件；不能随意换非对称 smoother。保持 real double、原始绝对/相对容差、warm start 和最终 `b-Ax` 真残差检查。性能不足时再研究混合精度，它不属于最小迁移。APU 共用物理内存不意味着 PETSc host/device 缓冲零复制；不以 AI/FP16 吞吐量推断 FP64 稀疏求解性能。

## 4. 分阶段实施与通过条件

### P0：验证 ROCm 平台

确认准确 CPU/GPU SKU、Windows 驱动发行版、WSL 版本与目标 ROCm 组合，安装匹配用户态组件。验收：rocminfo 枚举 GPU，hipcc 编译运行 FP64 kernel，HIP 查询设备架构/可用内存，双精度 BLAS 和稀疏 SpMV 冒烟成功。若失败，先解决环境，不进入应用改造。不采用伪造 GPU 架构作为默认解法。

### P1：隔离构建 PETSc/Hypre HIP

保留现有安装。下面仅是待环境就绪后的配置草案，未执行，ROCm 路径、架构和 BLAS 设置需按实际安装确认：

```bash
# 在独立 PETSc 3.25.5 源码副本中运行
./configure PETSC_ARCH=arch-hip-opt \
  --prefix=/home/midway/opt/petsc-3.25.5-openmpi-hip-opt \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpifort \
  --with-debugging=0 --with-scalar-type=real --with-precision=double \
  --with-hip=1 --with-hip-dir=/opt/rocm --with-hip-arch=gfx1151 \
  --download-hypre \
  COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3'
```

构建完成后先运行 PETSc 单进程 HIP Mat/Vec/KSP 示例与 Hypre 示例；检查配置宏、实际库路径和对象类型。保留当前索引宽度，不在同一次迁移修改精度或 MPI 实现。[PETSc 构建指南](https://petsc.org/release/install/install/)

**MPI 特别限制：**本地 `config/BuildSystem/config/packages/hypre.py:105` 在 HIP 分支添加 `--enable-gpu-aware-mpi`。现有 Open MPI 4.1.2 是否支持 ROCm 设备指针尚未证实。首轮单 rank；可用 PETSc 初始化选项 `-use_gpu_aware_mpi 0` 关闭其设备指针通信假设，但必须独立处理 Hypre，不能假定该开关控制 Hypre。已安装 HYPRE 头文件提供 `HYPRE_SetGpuAwareMPI(0)`，实施时需确认调用时序和版本行为，或构建禁用该能力的 Hypre。多 rank 后续比较 host staging 与 ROCm-aware MPI/UCX。PETSc 文档给出的后者构建路线包括 HIP + Open MPI + ROCm-enabled UCX。

### P2：最小应用原型

先实现全 CPU/全 HIP 的单 rank 小型线性系统开关，只用 None/Jacobi 验证矩阵、向量及收敛。对照已有后端测试中的矩阵数值变化、RHS-only、warm start、标量/向量三 RHS、非收敛与真残差约定；针对 HIP 路径补充有意义的测试。

随后为压力 equation identity 开启 HIP + GPU Hypre，其余方程保留 CPU Mat/Vec。在 HIP 版 Hypre 共用进程时确认 CPU 方程的 PC 行为；不得仅从 Mat 类型推断整个 PC 的执行位置。显式 GPU 请求在 CPU-only 编译或无设备时给出清晰错误；不允许静默回退后继续报告 GPU 成绩。

```bash
# 代码适配完成后的独立应用构建，当前执行此命令本身不会切换对象类型
make BUILD=build-petsc-hip \
  PETSC_DIR=/home/midway/opt/petsc-3.25.5-openmpi-hip-opt -j4
```

WSL 当前仅约 11 GiB 可见内存，编译先限并行度，是否调整 WSL 内存根据实际峰值决定。

### P3：数值验收与性能去留

先 heat/Poisson、小 cavity，再 NACA0012 127,013 单元的五步短算。相同网格、时间步、输入、容差、输出频率；对照真实残差、失败次数、守恒、有限值和所有字段。NACA 可沿用仓库 `compare_parallel_results.py --atol 5e-6 --rtol 5e-6` 作为既有场比较门槛，同时记录逐字段 L2/最大绝对差与质量误差。它是启动短程回归标准，不替代长期物理验证；失败时分析数值原因，不能为通过而自动放宽容差。

性能至少比较：CPU 1 rank、HIP 1 rank；再比较经过合理调优的 CPU 2/4 rank。同算法同容差比较与 GPU 最优参数比较分开报告。预热至少一次、重复至少三次，记录中位数和波动，分别报告初始化、CPU 装配、矩阵上传、PC setup、KSP、真实残差、解读回及整步/端到端时间。GPU 异步执行下现有 chrono 阶段计时可能把成本归到后续同步点；解读回也不在当前 KSP 计时里，应以包含完成同步的完整 solve/整步墙钟为主，并辅以可用的 PETSc GPU 事件/日志。

确认真实 HIP Mat/Vec、Hypre 构建/执行模式及 GPU 活动，不能只看 binary 链接了 HIP。压力局部加速但整步没有稳定收益时保留 opt-in，不切默认。建议整步中位数至少改善 15%、超过噪声且数值验收通过才讨论默认启用；15% 是工程决策门槛，不是预测。

### P4：收益成立后的扩展

再研究其他方程 GPU PC、设备驻留的 RHS/解、装配 kernel 和 MPI。单 APU 多 rank 共享 GPU 可能增加上下文、同步和内存成本，不能套用 CPU rank 数。不要首轮引入 Kokkos、全场 GPU 化、异步流水线或混合精度。

## 5. 性能判断边界

仓库现有 [NACA CPU 计时](../../../cases/naca0012/validation/petsc_one_step_timing.md) 中，旧容差下 CG+Hypre 单步计算约 1.745 s，其中压力 KSP 约 0.755 s。用这组历史样本作解释性 Amdahl 估算，若只优化该部分，比例约 43.3%，即使它无限快也仅约 1.76 倍整步；局部 3 倍时约 1.40 倍。这里没计迁移成本，也没有假定 GPU PC setup 收益。

这些数值来自历史短程、旧容差、2 rank 记录，不是本次实测预测；当前正式容差已收紧，必须重测。同记录中原 BJacobi 9.058 s 与 CPU Hypre 1.745 s 的差异是预条件器选择收益，不能作为 GPU 加速基线。

实际收益受 FP64 吞吐、稀疏带宽、AMG 层次、小网格 launch 成本、同步、host/device 复制、CPU/GPU 共享带宽及功耗影响。首先测压力迁移，能够以最小改动回答这台机器是否值得继续投入。

## 6. 投入估计与决策

以下为工程估计，不含未知驱动兼容问题：平台和 PETSc 冒烟约 0.5–2 人日；设备开关/压力原型约 1–2 人日；数值和性能对照约 1–2 人日。总体约 2.5–6 人日，MPI/其他 PC/长期验证另外评估。

当前判断：**架构可行；平台待就绪；GPU 正确性待测；性能收益未知。下一步应完成 P0/P1，再实施压力方程 HIP 原型。** CPU 版本继续作为参考和回退，独立前缀与构建目录可避免迁移试验影响现有使用。
