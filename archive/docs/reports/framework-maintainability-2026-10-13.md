# 2026-10-13 框架可维护性重构与性能对照

本次整理 `src/physics` 以外的通用框架实现，目标是让掌握 C++17 基础、MPI 和基本数值方法的开发者，
能够沿实际调用链理解、修改和验证代码。重点是职责、数据含义和控制流程，而不是压缩行数。

基线提交为 `c23d6d46aeaff379b2e2d2db8f01f4240bc90127`。修改位于本地工作区，
涉及 8 个原有源码/内部头、7 个新增源码/私有头，以及构建清单和开发文档。
公开头的接口、Physics 源码、数值方法、线性容差、PC 复用预算和 MPI 通信流程保持原有契约。
任务开始前已有改动的 8 个文件经 SHA-256 复核保留。

## 结构调整

| 维护问题 | 调整后的入口 | 保留的约束 |
| --- | --- | --- |
| 一个离散文件混合显式计算、梯度重建和隐式矩阵装配 | [operators.cpp](../../../src/discretization/operators.cpp)、[gradients.cpp](../../../src/discretization/gradients.cpp)、[equation_assembly.cpp](../../../src/discretization/equation_assembly.cpp) | 面方向、积分系数、边界规则、遍历和浮点运算顺序 |
| 显式计算和装配各自包含面采样规则 | [operator_support.h](../../../src/discretization/operator_support.h) | 真正共用的短函数保留 inline；通信由后端完成 |
| 标量/向量梯度重复构造最小二乘边界样本，常量/场容量重复时间项公式 | `leastSquaresSample`、`addCellTimeDerivative` | 右端和结果写回清晰可见；Euler/BDF2 公式及各自容量检查保留 |
| 网格构造和长串校验交织 | [mesh.cpp](../../../src/core/mesh.cpp)、[mesh_validation.cpp](../../../src/core/mesh_validation.cpp) | 几何计算分阶段；长度检查先于数组访问，失败次序与校验临时容器的生命周期保留 |
| Case 初始化检查嵌套较深，字段名字和初始化来源分开存储 | [case.cpp](../../../src/io/case.cpp) | `FieldRecord` 将名字和来源放在一起；场仍由稳定的类型容器拥有 |
| 时钟推进混入 Case 的配置、场和输出实现 | [time_step.cpp](../../../src/io/time_step.cpp) | 推进、接受、拒绝和历史恢复位于同一模块 |
| 模块注册、校验、描述输出和教程生成混在一个文件 | [model_description.cpp](../../../src/io/model_description.cpp)、[model_description_output.cpp](../../../src/io/model_description_output.cpp)、[case_template.cpp](../../../src/io/case_template.cpp) | 共享词汇和校验，描述和模板保持一致 |
| PETSc 求解函数同时负责初值、残差尺度、Krylov 续算和陈旧 PC 恢复 | [petsc_linear_system.cpp](../../../src/backend/petsc_linear_system.cpp) | 暖启动仍批量归约两个范数；三次 Krylov 求解共享原预算，PC 最多再尝试一次 |
| 分区候选使用位置含义不直观的五元组，缓存失效与迁移流程难以定位 | [mesh_partition.cpp](../../../src/parallel/mesh_partition.cpp) | 候选字段命名，比较次序和选择策略保留；队列版本与连通性缓存版本区分 |

原来 1351 行的离散文件现在由 473 行的显式算子、242 行的梯度实现、460 行的装配实现和一个私有共用头承担。
模块描述主文件从 842 行变为 480 行，输出和模板分别为 278 和 119 行。
源码总行数略有增加，增加的是边界、命名和必要说明。

| 主流程 | 重构前行数 | 重构后行数 |
| --- | ---: | ---: |
| `Mesh::polyhedral` 主构造流程 | 278 | 191 |
| `Mesh::validate` 校验入口 | 189 | 16 |
| `Case::finalizeInitialization` | 77 | 10 |
| `PetscLinearSystem::solveRhs` | 145 | 82 |
| `QuotaBalance::apply` | 80 | 14 |

这里的行数只说明主流程和职责分离；具体几何、校验和求解逻辑仍存在于命名明确的小函数中。
修改导航、数学约定和验证方法见[框架源码阅读与修改](../../../docs/developer-guide/framework-maintenance.md)。

## 编译与行为验证

重构前后均执行正式构建、`test-full`、独立 `test-piso-coupling`、`petsc-check` 和 `test-descriptions`，全部通过。
`test-full` 包括串行契约、工作流、外部接口、RANS、生产标量时间阶、SIMPLE 串并行、MPI 和 Poiseuille 组。
最后保留网格校验临时容器原有生命周期后，再次通过正式构建、网格/几何测试、PETSc 串并行检查及 2/3/4 进程分区测试。
另外用 `USE_METIS=0` 构建分区测试，并在 2/3/4 个进程下验证，全部通过。

架构检查覆盖 94 个源码/头文件，include 无环且分层约束通过。
clang-format 14 检查全部 50 个非 Physics C++ 文件通过；正式框架构建没有编译警告。
文档链接检查和 `git diff --check` 通过。

直接运行保留的基线和新程序进行额外逐字节对照：8 个模块各自的 Markdown/JSON 输出一致；
heat/transport 生成的教程文件一致；simple/transientSimple/piso 不支持通用模板时的错误行为一致。
共完成 21 组描述或模板行为比较。

## 性能验证口径

保留重构前后程序，用同一组编译优化、MPI/PETSc/METIS、线程环境和输入进行顺序交替测量。
机器信息、二进制和输入哈希、原始逐 rank JSON、stdout 和最终场均保存于任务证据目录。
六项编译变量已直接比对：`CXX`、`AR`、`OPTFLAGS`、`CXXFLAGS`、`PETSC_OPTFLAGS`、`CPPFLAGS` 一致。

两组 54000 单元输入分别固定为 PISO 的 300 个接受时间步和 SIMPLE 的 600 次外迭代。
使用 1/4 个进程，每种程序/算例/进程数组合有 1 次预热和 3 次正式重复，共 32 次运行。
主要指标为 `max_rank(iterationComputeSeconds)`：从第一次完成的迭代之后计时，扣除窗口内输出，
排除启动、读取、网格构造、分区和第一次求解。另保留算例初始化和应用墙钟时间供核对。
SIMPLE 的固定迭代耗尽是吞吐测试，不表示稳态物理收敛。

**性能对照尚在运行，结果和字段/工作量核对完成后更新此段。**

## 证据与复现

完整工作目录为 `build/maintainability-20261013`，其中保存基线源码、前后程序、编译和回归日志、
输入副本、原始测量、最终场和哈希清单。该目录不会改变或清理原生产算例的运行状态。
大型网格、程序和逐单元字段保留在本机证据目录，不加入版本化报告。

运行回归：

```bash
make -j4 all
make -j4 test-full
make -j4 test-piso-coupling
make -j4 petsc-check
make -j4 test-descriptions
python3 tests/architecture_test.py
python3 tools/check_docs.py
```

本次重构保留原有数值语义。回归、字段一致和固定工作量计时各有自己的证据范围，
不代替独立物理验证，也不覆盖所有网格、规模、机器或编译器下的性能。
