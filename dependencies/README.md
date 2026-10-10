# 依赖源码

PETSc、Eigen 和所需 Hypre 的必要源码集中放在本目录，随 BabelSim 一起交付。
版本、上游地址、提交和源码包 SHA-256 见 [sources.json](sources.json)，
精简范围与前后大小见 [pruning.json](pruning.json)。

| 路径 | 内容与使用方式 | 许可证 |
| --- | --- | --- |
| [petsc](petsc/) | PETSc 3.25.5，提交 `a0b507cc1cc9eac739c34d6e71521f4854e84e36`；保留 CPU/MPI/real/double/Hypre 库的 974 个编译单元、传递头文件和构建脚本，编译到 `build/petsc-build` 并安装到 `build/petsc` | [BSD-2-Clause](petsc/LICENSE) |
| [eigen](eigen/) | Eigen 3.4.0；保留 Core、Cholesky 及其传递头文件，包括 Jacobi 和 CPU 相关分支。直接编入 BabelSim，无独立二进制库 | [MPL-2.0](eigen/COPYING.MPL2)，另保留相关许可证 |
| [distfiles](distfiles/) | Hypre 3.1.0 加 PETSc 要求的修复，提交 `85b779557005b2eb94c231c1b516e988b87f4e53`；只保留 GNU configure 与 CPU 库源码的压缩包 | COPYRIGHT、LICENSE-MIT 和 LICENSE-APACHE 随源码包保留 |
| [build_petsc.py](build_petsc.py) | 本地源码配置、编译、安装与构建目录清理脚本；根目录 `make` 自动调用 | BabelSim MIT |
| [install.sh](install.sh) | 系统依赖检查、缺少包的安装和一键源码构建入口；Ubuntu/Debian 自动使用 apt-get | BabelSim MIT |

Hypre 的版本和提交按 PETSc 3.25.5 的 `config/BuildSystem/config/packages/hypre.py`
固定，构建时先校验源码包 SHA-256，再交给 PETSc 的 `--download-hypre=<本地路径>`。
依赖构建无需从网络获取 PETSc、Hypre 或 Eigen。

PETSc 源码来自官方 tag 的干净 Git 快照，按其 `gmakegen.py` 的构建选择和头文件依赖
精简，去掉文档、教程、测试、测试数据、开发工具、未启用的 GPU/外部包实现和 Python/
Fortran 绑定。保留的库文件内容不改，`.gitignore` 增加了交付与构建产物规则。
PETSc 的单库构建包含包初始化和注册依赖，因此仍保留所选 Sys/Vec/Mat/DM/KSP/
SNES/TS/TAO/ML 库源码，不能仅按 BabelSim 的直接 API 调用删除其中某个包。

Eigen 同样去掉文档、测试、示例、unsupported 和未使用模块；Hypre 去掉文档、测试、
示例、CMake 和 CUDA 源码；不保留完整源码包备份。没有嵌套 Git 仓库。
当前源码集针对上述 CPU/MPI 配置，启用 GPU、其他 PETSc 外部包或增加 Eigen 模块时
需要补齐相应源码。

首次准备机器时，在 BabelSim 根目录执行：

```bash
./dependencies/install.sh
```

脚本检查 GNU Make、GCC C/C++/Fortran、MPI 开发环境、Python >= 3.9、pkg-config、
tar/gzip 和 BLAS/LAPACK，并实际编译、链接和运行 MPI/C++17 与 BLAS/LAPACK 小程序，
核对 C/C++/Fortran 使用同一 MPI 实现和版本。
Ubuntu/Debian 缺少依赖时自动调用 apt-get 安装对应系统包，安装后再次检查；普通用户
仅在系统安装阶段使用 sudo。已有可用 MPI 环境时不安装另一套 MPI。
PETSc、Eigen、Hypre 均使用本目录中的源码，产物仍生成在 `build/`。

```bash
./dependencies/install.sh --check       # 仅检查，缺少依赖时返回非零
./dependencies/install.sh --deps-only   # 只准备系统依赖
./dependencies/install.sh --clean --jobs 4  # 清空 build 后从头编译
```

默认保留已有构建。依赖齐全时可继续直接使用根目录 Makefile：

```bash
make
make petsc-check
```

自动安装适用于 apt-get 环境，建议 Ubuntu 22.04+ 或 Debian 11+；其他系统可在准备好
上述构建环境后直接运行脚本。可选 METIS 可以使用系统安装，也可以使用现有图生长回退。
`petsc-check` 使用 BabelSim 的 PETSc/Hypre 后端回归检查 1/2/4 rank，不依赖被移除的
PETSc 示例目录。目录、清理、可移植编译和外部 PETSc 用法见 [构建说明](../docs/user-guide/installation.md)。
