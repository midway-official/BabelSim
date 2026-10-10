# 可运行示例

[文档首页](../README.md) · [用户教程](../user-guide/tutorials.md) · [开发教程](../developer-guide/extending-solvers.md)

所有命令从 BabelSim 根目录执行。示例不依赖已保存的研究结果。生成器仅创建新目录；各教程用自己的运行名，避免混写。

| 文件 | 作用 |
| --- | --- |
| [create_case.py](create_case.py) | 创建 heat/transport/simple/piso/custom-heat/coupled 的完整小型算例 |
| [external_heat.cpp](external_heat.cpp) | 完整热方程、Euler/BDF2、错误状态与 customHeat 注册 |
| [external_main.cpp](external_main.cpp) | 独立应用入口与错误报告 |
| [coupled_scalar.cpp](coupled_scalar.cpp) | 双场耦合、显式时间历史、派生 total 场 |
| [Makefile](Makefile) | 与当前框架依赖一致地构建示例程序 |
| [external_project.mk](external_project.mk) | 只有公开 SDK 的独立工程 Makefile |
| [one-cell.mesh](one-cell.mesh) | 完整合法的 BABELSIM_MESH 3 单单元文件 |
| [verify_examples.py](verify_examples.py) | 构建并运行手册示例，保留日志与 JSON 证据 |

## 构建与调用

```bash
make -f docs/examples/Makefile -j4
python3 docs/examples/create_case.py --template custom-heat --output work/customHeat
build/examples/custom-heat -case work/customHeat -time external
python3 docs/examples/create_case.py --template coupled --output work/coupled
build/examples/coupled-scalar -case work/coupled -time serial
```

产物为 build/examples/custom-heat、coupled-scalar、registered-heat。最后一个同时链接正式求解器和示例注册项，可用于验证同一启动器按 case 分派。

如果已经按开发教程把热示例安装到 `src/physics/custom_heat/main.cpp`，示例 Makefile 在 registered-heat 中排除该同名模块对象，再加入示例对象，避免重复注册。若在其他目录又注册 customHeat，应消除重名后构建。

## 独立 SDK 工程

```bash
make -f docs/examples/Makefile sdk-example
build/examples/sdk-project/custom-heat -case work/customHeat -time sdk
```

生成的 sdk-project 包含 `sdk/include/`、`sdk/lib/libbabelsim.a`、两个源文件与独立 Makefile；不包含 src/internal、Eigen 或 MPI/PETSc 头。源文件使用标准 g++ 编译，MPI/PETSc 等仅在链接阶段加入。仍需要运行库，不代表可脱离 PETSc/MPI 独立执行。

搬到其他目录后，可按实际安装设置：

```text
make PETSC_DIR=/absolute/path/to/petsc PARTITION_LIBS='<框架启用的 METIS 链接参数>'
```

未启用 METIS 时 PARTITION_LIBS 留空。SDK ABI 与 CPU 指令集需匹配构建环境。生成 SDK 工程不是系统安装过程，不需要 sudo。

## 自动验收

```bash
python3 docs/examples/verify_examples.py --output build/docs-refresh
```

报告为 build/docs-refresh/validation.json，独立输入与逐命令日志保存在 tutorials-* 子目录。验证覆盖六类算例 × 1/2/4 进程，以及 SDK、注册分派、Euler/BDF2、输出序列、继续计算、网格格式、拒绝覆盖和未知 solver 错误。

范围是软件教程与 API/I/O 验证。小型网格的功能演示不是模型精度、强扩展或工程物理认证。
