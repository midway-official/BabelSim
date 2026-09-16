# BabelSim 验证与维护检查

本文记录当前源码应执行的验证入口。历史性能和物理研究报告保存在 docs/reports/，它们是带日期
的证据，不是当前 DSL API 的规范；新开发者以本目录的架构、Case、DSL 和 Solver 指南为准。

## 1. 快速回归

~~~bash
make -j4 all
make -j4 test
~~~

当前构建实测通过：

- architecture_test：生产 include 闭包、层次边界和 Physics 禁止依赖检查通过；
- procedural_equation、operators、field、mesh、backend、numerical contract、monitor
  等基础测试通过；
- case_lifecycle、case_io、field_writer：Case 稳定引用、load/create/existing、typed
  参数、单场线性覆盖和 output.bs 字段筛选通过；
- heat_solver、transport_solver、time_history：标量 PDE、Transport、Euler/BDF2 和末步
  时间元数据通过；
- simple_solver、cavity、3D cavity、二维/三维 non-orthogonal cavity：SIMPLE 外循环、
  质量平衡、Rhie–Chow 和压力修正回归通过。

代表性本次输出：

| 测试 | 结果 |
| --- | --- |
| simple_solver_test | 62 次迭代，质量相对不平衡 3.14e-17 |
| cavity_regression_test | 153 次，质量 6.99e-17 |
| cavity_3d_test | 66 次，质量 5.61e-15 |
| nonorthogonal_cavity_test | 178 次，质量 1.51e-15 |
| nonorthogonal_cavity_3d_test | 116 次，质量 3.03e-15 |

这些数字验证当前实现没有因 API 重构改变数值路径；它们不代表所有网格和雷诺数的精度。

## 2. 工作流、外部 Solver 和 RANS

~~~bash
python3 tests/solver_workflow_test.py
python3 tests/external_solver_test.py
python3 tests/rans_validation_test.py
~~~

工作流测试覆盖 Heat/Transport/耦合时间序列、结果读取和失败路径。外部 Solver 测试在仓库
外只使用 include/ 和预编译库，确保 Physics 不需要 MPI/Eigen/内部头。RANS 脚本覆盖 SA、
k–omega、k–epsilon 的 Euler/BDF2 和稳态/瞬态配置。

## 3. MPI 和并行结果

~~~bash
make test-mpi
make test-mpi-poiseuille
~~~

并行检查至少使用 1、2、4 rank，比较 global cell ID 对齐后的场值、质量平衡和收敛状态。
不同进程数不要复用同一结果目录；Case 会拒绝分区数不一致的旧目录。MPI 数值一致性不等于
逐 bit 可复现，也不等于大规模扩展性证明。

## 4. 架构门禁

~~~bash
python3 tests/architecture_test.py
git diff --check
~~~

门禁应拒绝 Physics 对 src/internal、runtime、MPI、Eigen、CSR、原始 Field 存储和其它
Solver 私有头的依赖。它还检查公共头可独立包含以及旧的含糊 API 没有重新出现。

## 5. 新 Solver 的验证最低线

新增物理求解器时至少增加：

1. typed 参数缺失、重复、错误范围和未消费配置测试；
2. load/create/existing 的明确生命周期测试，确认初值没有被静默忽略；
3. 一个解析解或制造解；
4. 稳态和瞬态（若适用）的结构化 SolverResult；
5. 线性不收敛、NaN/Inf、外迭代耗尽的失败路径；
6. 串行和 1/2/4 rank 的收敛、输出字段和时间目录检查；
7. 若使用 RANS 或非正交格式，模型/修正子循环的独立回归。

不要把未收敛结果作为物理精度结果。报告中区分源码存在、调用链接入、回归通过、短跑稳定
和物理验证五种状态。
