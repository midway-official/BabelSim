# 验证与维护检查

测试入口以 [Makefile](../Makefile) 与测试源码为准。现行用法见 [文档中心](README.md)；
旧审查和验证报告在 [历史记录](history/README.md)，仅证明对应时期的工作。

## 选择测试组

| 命令 | 覆盖 |
| --- | --- |
| `make test` | 架构门禁、几何/边界/算子/配置/方程/后端/IO 与本地 C++ 回归 |
| `make test-workflow` | 生产 Heat/Transport、自定义双场 Solver、时间序列与失败路径 |
| `make test-external` | 只用公开头和静态库、1/2/4 rank、负向 API、无 MPI 结果读取器 |
| `make test-rans` | SA/k–ω/k–ε 方程、时间推进、并行及失败处理 |
| `make test-scalar-time` | 生产标量 Euler/BDF2 时间阶、显式修正和缩短末步 |
| `make test-simple-parallel` | cavity/channel/cube/warped，steady/Euler/BDF2，1/2/4 rank |
| `make test-piso-coupling` | PISO 校正与冻结动量缺陷、耦合接受和失败路径；独立专用组 |
| `make test-mpi` | 分区/halo、数学与算子、PETSc、数值与 Physics 契约、流动/标量；含 Heat 对照 |
| `make test-mpi-heat` | 生产 Heat 的 1/2 rank 字段对照；已被 test-mpi 调用 |
| `make test-mpi-poiseuille` | 生产 Poiseuille 的 1/2/4 rank、字段与后处理 |
| `make test-architecture` | 项目头依赖、公开接口和 Physics 分层门禁 |

`make test-full` 严格按以下顺序执行：test → test-workflow → test-external → test-rans →
test-scalar-time → test-simple-parallel → test-mpi → test-mpi-poiseuille。
**PISO 专用 test-piso-coupling 不在该列表里，需另行执行。**
Heat/Transport 的局部装配手算在 `scalar_equation_test`，History 生命周期与重复内迭代
在 `time_history_test`；它们不替代生产 Solver 的时间阶测试。

## 执行方式

```bash
make -j4 all
make test-full
make test-piso-coupling
python3 tools/check_docs.py
python3 tests/architecture_test.py
git diff --check
```

编译可用多个作业，测试组和求解任务按顺序执行，不同时启动另一个求解器或性能驱动。
切换编译选项/依赖需独立 BUILD；调试配置见 [构建页](building.md)。

物理基准快速检查逐项执行：

```bash
make validate-cavity
make validate-poiseuille
```

`make validate` 的依赖为 test、validate-cavity、validate-poiseuille，不含完整测试矩阵。
不要用 `make -j validate` 并行运行数值验证。Ghia 快速目标针对当前 Re=100 方腔；多 Re
和网格研究是 [历史报告](history/reports/cavity-ghia-validation.md) 的特定实验，不能把它们
全部算作一次快速目标已验证的内容。

## 证据分别说明

| 证据 | 能证明什么 |
| --- | --- |
| 源码/API 核对 | 定义和调用路径存在 |
| 构建、回归、失败路径 | 测试覆盖的执行和契约成立 |
| 解析解、制造解、阶数研究 | 指定问题与细化范围的数值精度 |
| 同输入 MPI 字段对照 | 给定网格、时刻和容差下的一致性 |
| 持续运行和强扩展 | 相同工作量下的实测速度与效率 |
| 独立物理基准 | 指定物理范围的定量适用性 |

有限值、线性收敛、质量误差小或一次数值推进完成，均不能单独证明工程物理正确。
MPI 对照不要求逐 bit 一致，也不证明强扩展。固定工作量 SIMPLE 允许以明确的迭代耗尽
状态测吞吐，必须单独说明尚未达到稳态判据，不能作为收敛时间或稳态精度数据。

## 改动与最低检查

配置/IO 变更检查合法输入、缺失/重复/未知/未消费键、初场读取、输出与重启。
公开接口变更检查 architecture、external、workflow 和全部调用方。离散/后端/分区
变更增加对应局部契约、守恒、1/2/4 rank 对照；涉及 PISO 执行专用耦合组。
时间格式变更用实际生产 Solver 检查 Euler/BDF2、启动与末步，不只检查一个时间系数函数。

新 Solver 至少覆盖 typed 参数、load/create/existing 生命周期、一个解析/制造解、
结构化成功/失败状态、线性不收敛/NaN/外迭代耗尽、同一时刻的串并行字段与输出。
RANS 或非正交子循环还需独立验证。长期物理结论需适当的网格/时间步研究与采样。

最近工程优化的完整回归日志、源码哈希、字段及测量审计见
[优化报告](reports/flow-optimization-2026-10-07.md)。文档变更不能改写历史测试日志；
源码快照相符才可引用其结果，最新文档检查单独记录。
