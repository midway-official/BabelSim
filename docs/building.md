# 构建、依赖与运行

本页对应 [Makefile](../Makefile)。依赖或编译参数变化时使用新的 `BUILD` 目录；Make
不跟踪命令行参数变化，旧对象不能证明新配置生效。

## 依赖和安装路径

| 依赖 | 用途 |
| --- | --- |
| GNU Make、GCC C++17、`gcc-ar` | 构建与 LTO 静态库 |
| `mpic++`、`mpirun` | MPI 应用；与 PETSc 使用同一 MPI 实现 |
| Eigen 3 | 局部最小二乘，默认 `/usr/include/eigen3` |
| PETSc 头、库、`PETSc.pc` | 分布式 Mat/Vec/KSP/PC |
| Hypre（经 PETSc） | 选择 `pcType hypre` 的方程 |
| METIS（可选） | 多约束图划分 |
| Python 3 | 架构、工作流、结果比较、性能脚本 |

本机证据使用 PETSc 3.25.5，前缀 `/home/midway/opt/petsc-3.25.5-openmpi-opt`。
这也是 Makefile 默认 `PETSC_DIR`，其他机器需设置实际前缀。
Makefile 查询 `<PETSC_DIR>/lib/pkgconfig/PETSc.pc`；源码树多架构安装需指定实际包含
该目录的前缀，或显式提供匹配的 flags。

```bash
export PETSC_DIR=/path/to/petsc-install
PKG_CONFIG_PATH="$PETSC_DIR/lib/pkgconfig" pkg-config --modversion PETSc
PKG_CONFIG_PATH="$PETSC_DIR/lib/pkgconfig" pkg-config --cflags --libs PETSc
make -j4
```

PETSc 必须含示例所需的 Hypre。默认产出 `build-petsc/libbabelsim.a`、
`build-petsc/babelsim-solve`、`build-petsc/babelsim-post`；默认目标不执行测试。

## 编译配置

默认优化为：

```text
-O3 -march=native -mtune=native -flto=auto -ffat-lto-objects
-ffast-math -fno-finite-math-only -ffp-contract=off -DNDEBUG
```

使用 double，保留显式 NaN/Inf 检查；快速数学仍可改变舍入、极小数和溢出行为。
异构集群需选共同 CPU 基线。fat LTO 保留普通对象码，外部可用 `-fno-lto` 链接。
不启用快速数学、使用编译器默认 CPU 基线的独立构建：

```bash
make -j4 BUILD=build-portable OPTFLAGS='-O3 -flto=auto -ffat-lto-objects'
```

`CXX`、`AR` 可覆盖，MPI/PETSc、C++ ABI、LTO 工具需匹配。清理只删除所选构建目录：

```bash
make BUILD=build-portable clean
```

`make -j4 debug` 使用独立 `build-debug/`、`-O1 -g3`、ASan/UBSan、保留栈帧，无 LTO 和
快速数学。它只构建。需要实际检查时，测试使用相同的调试 flags：

```bash
make -j4 debug
make BUILD=build-debug \
  OPTFLAGS='-O1 -g3 -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=all' \
  LDFLAGS='-fsanitize=address,undefined' test
```

## 可选 METIS

| 参数 | 含义 |
| --- | --- |
| `USE_METIS=auto` | 默认，找到链接配置时启用，否则图生长 |
| `USE_METIS=1` | 要求依赖可用 |
| `USE_METIS=0` | 显式图生长 |
| `METIS_DIR` | 安装前缀，默认 `$HOME/opt/metis-5.1.0` |
| `METIS_CFLAGS`、`METIS_LIBS` | 覆盖头和链接 flags |

```bash
make -j4 BUILD=build-metis USE_METIS=1 METIS_DIR=/path/to/metis-install
make -j4 BUILD=build-growth USE_METIS=0
```

头和库需 ABI 匹配。非系统位置的 `libmetis.so` 还需可被动态加载器找到。
策略说明见 [网格图模型](performance/partition-model.md)。

## 外部 Solver

源码只需公开头。静态库链接仍需 MPI/PETSc；启用 METIS 时需同一 `METIS_LIBS`。
下例针对 `USE_METIS=0` 的库，使用已测试的仓库示例：

```bash
mpic++ -std=c++17 -O2 -fno-lto -Iinclude \
  tests/examples/coupled_scalar.cpp build-petsc/libbabelsim.a \
  $(PKG_CONFIG_PATH="$PETSC_DIR/lib/pkgconfig" pkg-config --libs PETSc) \
  -Wl,-rpath,"$PETSC_DIR/lib" -o /tmp/babelsim-coupled
```

启用 METIS 时在 `-o` 前加实际 flags，例如 `-L/path/to/metis-install/lib -lmetis`，并
提供其运行时路径。`make test-external` 按当前配置带入依赖，检查 1/2/4 rank 和负向 API。
内置扩展放入 `src/physics/<name>/main.cpp` 并注册，Makefile 自动收集。

## 运行与退出状态

```text
babelsim-solve -case <case-directory> [-time <run-name>] [-performance <directory>]
babelsim-post -case <case-directory> [-time <selection>] -format <vtk|tecplot> [...]
```

运行名是单个非数字目录名，不能包含路径分隔符或为 `.`/`..`，不能使用 `all`/`latest`。
`demo`、`mpi4` 均可。`-time` 不选择重启时刻；重启需配置初场，见 [文件格式](mesh-and-results.md)。
退出码 `0` 表示 Solver 完成，`2` 表示迭代耗尽或数值失败，`1` 表示配置/应用或结束阶段
错误。固定迭代吞吐只在达到目标工作量、状态 `maxIterations` 且退出码 `2` 时接纳。

`COMPUTE_BACKEND_SOURCES` 在构建期替换内部后端工厂；不承诺动态插件或稳定后端 ABI。
测试见 [验证](validation.md)，持续运行规则见 [性能](performance/README.md)。
