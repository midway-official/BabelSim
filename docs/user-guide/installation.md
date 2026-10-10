# 02 安装与构建

[用户手册](README.md) · 上一章：[产品与工作流程](overview.md) · 下一章：[首次运行](quick-start.md)

## 系统要求

正式构建使用 GNU Make、支持 C++17 的 GCC、MPI C/C++/Fortran 开发环境、Python 3.9 或更新版本、pkg-config、tar/gzip、BLAS 和 LAPACK。Fortran 编译器用于依赖构建，BabelSim 应用代码使用 C++。MPI 编译器包装器必须来自同一实现与版本。

自动安装脚本支持 Ubuntu/Debian 的 apt-get 环境。当前实际验证环境是 Ubuntu 22.04、GCC 11.4、OpenMPI 4.1.2、Python 3.10、GNU Make 4.3 与 OpenBLAS。其他 Linux 环境需要先准备这些系统依赖，再执行源码构建。Windows 用户可以在具备上述工具的 WSL Linux 环境中使用；没有提供原生 Windows 二进制安装器。

## 推荐的一键入口

进入交付的 BabelSim 目录，执行：

```bash
./dependencies/install.sh
```

脚本实际编译、链接和运行检查程序，核对 C++17、MPI 及 BLAS/LAPACK。缺少系统包时，在 apt-get 环境安装所需包，然后重新检查并执行 `make`。普通用户只在系统安装阶段需要 sudo 权限；依赖与 BabelSim 的本地构建不需要 root 权限。

| 命令 | 行为 |
| --- | --- |
| `./dependencies/install.sh --check` | 只检查，缺少依赖返回非零 |
| `./dependencies/install.sh --deps-only` | 检查并安装缺少的系统依赖，不构建 BabelSim |
| `./dependencies/install.sh --jobs 4` | 完整检查并使用 4 个构建作业 |
| `./dependencies/install.sh --clean --jobs 4` | 清除默认 `build/` 构建产物后重新编译 |
| `./dependencies/install.sh --help` | 显示安装脚本帮助 |

`--check`、`--deps-only`、`--clean` 互斥；作业数必须为正整数。脚本默认保留已有构建。缺少包的安装需要系统软件源可用；已准备好系统依赖后，PETSc、Eigen 和 Hypre 的构建使用交付目录内的源码。

## 使用 make

系统依赖齐全时，一条命令构建正式软件：

```bash
make
```

首次构建自动从 `dependencies/petsc` 配置并编译 PETSc，使用 `dependencies/distfiles` 中经过 SHA-256 校验的 Hypre 源码包。Eigen 为头文件库，直接由 `dependencies/eigen` 使用。依赖源码版本、来源、许可证和精简范围见[依赖源码说明](../../dependencies/README.md)。

```bash
make -j4 PETSC_JOBS=4
make petsc-check
```

`-j4` 控制 BabelSim 的构建作业；`PETSC_JOBS=4` 控制 PETSc/Hypre 构建。`petsc-check` 在 1、2、4 个进程下检查实际 PETSc 后端，并不是重新下载或运行已经移除的 PETSc 教程。

## 产物与目录

```text
dependencies/                 交付的必要依赖源码与安装脚本
build/
  petsc-build/                PETSc/Hypre 配置和编译中间产物
  petsc/                     本地 PETSc 安装目录及 Hypre 库
  libbabelsim.a               BabelSim 框架静态库
  babelsim-solve              正式求解程序
  babelsim-post               正式后处理程序
  core/, io/, backend/, ...   框架和求解器目标文件
```

PETSc 源目录可能包含指向构建目录的相对链接和少量配置状态；库与主要编译产物生成在 `build/`。不要将 `build/petsc-build` 当作需要交付给使用者的依赖源码。

重新构建默认目录：

```bash
make clean
make -j4 PETSC_JOBS=4
```

`make clean` 会删除对应构建目录，先将需要长期保留的测试数据、日志或用户算例移到构建目录外。它不清理 `cases/` 下的结果，也不卸载系统包。

## 可移植与调试构建

默认优化针对当前 CPU，包含本机指令集和浮点优化。若交付预编译程序到其他 CPU，需单独构建可移植版本，并重新运行相关验证：

```bash
make BUILD=build-portable USE_METIS=0 \
  OPTFLAGS='-O3 -DNDEBUG -ffp-contract=off' \
  PETSC_OPTFLAGS='-O3' -j4 PETSC_JOBS=4
```

该命令同时去掉 BabelSim 与 PETSc 的 `-march=native`；仅更改其中一方不足以保证跨机器兼容。动态 MPI、BLAS/LAPACK 和操作系统库仍需与运行机器兼容。默认浮点优化和上述构建可能有不同的舍入路径，应根据允许误差比较结果。

```bash
make debug
```

调试构建默认写入 `build-debug/`，启用地址与未定义行为检查，不替换正式构建。调试框架所链接的 PETSc 默认仍采用依赖脚本的优化配置。

## 可选设置与已有 PETSc

| Make 变量 | 默认/作用 |
| --- | --- |
| `BUILD` | `build`，选择本地构建与依赖安装根目录 |
| `PETSC_JOBS` | `4`，依赖构建并发数 |
| `CXX`、`PETSC_CC`、`PETSC_FC` | `mpic++`、`mpicc`、`mpifort` |
| `USE_METIS` | `auto`；找到 METIS 则使用，否则图生长分区回退 |
| `USE_METIS=0` | 禁用可选 METIS，减少外部交付依赖 |
| `USE_METIS=1` | 要求 METIS 可用，缺少时构建报错 |
| `EIGEN_DIR` | 默认交付目录中的 `dependencies/eigen` |
| `PETSC_DIR` | 默认对应 `BUILD/petsc`；显式指定时使用已有安装 |

使用已有 PETSc 是开发者的可选路径：

```bash
make BUILD=build-external PETSC_DIR=/absolute/path/to/petsc -j4
```

该安装必须提供 `lib/pkgconfig/PETSc.pc`，具备与框架一致的实数双精度/MPI 配置，并提供所选预条件器。交付时记录 PETSc、MPI、编译器和可选 METIS 信息。修改编译器、优化参数或外部库后，使用新的 `BUILD` 或清理对应目录，避免复用旧目标文件。

## 安装完成后的检查

确认 `build/babelsim-solve`、`build/babelsim-post` 和 `build/libbabelsim.a` 存在，并完成[首次运行](quick-start.md)。启动失败时先检查 `ldd build/babelsim-solve`；MPI 或动态库问题见[故障排查](troubleshooting.md)。源码目录就是当前的使用位置，没有另行提供 `make install` 或 `make uninstall` 目标。
