# BabelSim 验证与维护检查

本文是测试与验证入口的维护清单。当前测试源码以 `Makefile` 的 `TEST_SOURCES` 和专用
`test-*` 目标为准；历史性能和物理研究报告保存在 [reports/](reports/)，仅作为对应提交的证据。
接口与用法以 [DSL 与运行时用户手册](dsl-runtime-manual.md) 和
[架构与维护边界](architecture.md) 为准。

## 测试分组

| 入口 | 覆盖范围 |
| --- | --- |
| `make test` | 架构门禁和本地 C++ 回归：网格与几何、边界、算子与数学、配置与方程契约、后端、Case/字段 I/O、标量方程、SIMPLE 与腔体回归 |
| `make test-workflow` | 内置 Heat/Transport 和自定义耦合 Solver 的完整运行工作流、输出时间序列与失败路径 |
| `make test-external` | 仓库外只依赖公开头和库构建 Solver；检查 1/2/4 rank、错误 API 用法和 MPI-free 结果读取 |
| `make test-rans` | SA、k–omega、k–epsilon 方程项、时间推进、并行和失败处理 |
| `make test-scalar-time` | 生产 Solver 的标量时间精度、Euler/BDF2 与末步步长 |
| `make test-simple-parallel` | cavity/channel/cube/warped 的 steady/Euler/BDF2 在 1/2/4 rank 间的一致性 |
| `make test-mpi` | 分区、halo、并行数学/算子、PETSc、SIMPLE 和标量输运；包含 Heat 的 1/2 rank 对照 |
| `make test-mpi-poiseuille` | Poiseuille 案例的 1/2/4 rank 启动、结果对照和后处理 |

Heat 与 Transport 的一单元方程装配检查合并在 `scalar_equation_test` 中；完整生产求解流程仍由
`test-workflow` 覆盖。`time_history_test` 单独检查 History 生命周期和重复内迭代，生产级时间阶由
`test-scalar-time` 检查。

## 常用命令

快速本地回归：

~~~bash
make -j4 test
~~~

完整测试矩阵按组顺序执行，避免不同工作流并行写入相同算例结果目录：

~~~bash
make test-full
~~~

`test-full` 按顺序包含上表全部入口。物理基准快速检查另由下列目标执行：

~~~bash
make validate-cavity
make validate-poiseuille
make validate
~~~

`make validate` 运行 `make test`、Ghia 方腔检查和 Poiseuille 解析解比较；它不替代完整测试矩阵。
MPI 数值一致性不代表逐 bit 可复现，也不构成大规模扩展性证明。不同进程数不要复用同一结果目录；
Case 会拒绝分区数不一致的旧目录。

## 架构门禁

~~~bash
python3 tests/architecture_test.py
git diff --check
~~~

门禁检查生产 include 闭包、分层边界、Physics 禁止依赖和公共头的独立包含。它还防止已移除的临时
装配/求解包装、旧配置别名和隐藏 `Case::loop()` 重新出现。

## 新 Solver 的验证最低线

新增物理求解器时至少增加：

1. typed 参数缺失、重复、错误范围和未消费配置测试；
2. load/create/existing 的生命周期测试，确认初值没有被静默忽略；
3. 一个解析解或制造解；
4. 稳态和瞬态（若适用）的结构化 SolverResult；
5. 线性不收敛、NaN/Inf、外迭代耗尽的失败路径；
6. 串行和 1/2/4 rank 的收敛、输出字段和时间目录检查；
7. 若使用 RANS 或非正交格式，模型/修正子循环的独立回归。

新增长期运行或多 rank 脚本时，为其选择明确的 `test-*` Make 目标，并同步本页入口表。新增小型
C++ 回归时，将源文件加入 `TEST_SOURCES`，并按功能归入对应测试分组。不要把未收敛结果作为物理精度
结果；报告中区分源码存在、调用链接入、回归通过、短跑稳定和物理验证五种状态。
