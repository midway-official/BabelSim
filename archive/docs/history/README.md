# 历史开发记录

本目录只保存历史开发文件：过时接口、重构计划、执行提示、审查、修复记录和当时的
验证/性能报告。**这些文件不作为当前框架的使用说明、配置规范或性能承诺。**
现行入口为 [文档中心](../../../docs/README.md)，最近实验见 [版本化报告](../reports/README.md)。

## 归档范围

2026-10-08 将旧 `docs/reports/` 中的早期报告、附属证据和图像迁入 `history/reports/`，
将前提过时且已执行的 polyhedral 重构计划和提示迁入 `history/plans/`。
原位置和新位置见 [迁移清单](archive-map.json)。保留原数字、结论、命令和原始数据；
只增加历史标识并修正阅读链接。冻结源码、日志和机器可读证据保持原样。

| 历史主题 | 记录 | 当前入口 |
| --- | --- | --- |
| polyhedral 计划与执行提示 | [计划](plans/polyhedral-fvm-refactor-plan.md)、[提示](plans/polyhedral-fvm-refactor-prompts.md) | [架构](../../../docs/developer-guide/architecture.md)、[网格格式](../../../docs/reference/file-formats.md) |
| 早期框架审查与 F1–F7 | [审查](reports/framework-review-2026-09-09.md)、[实现](reports/f1-f7-implementation.md) | [架构](../../../docs/developer-guide/architecture.md)、[验证](../../../docs/developer-guide/verification.md) |
| polyhedral 阶段门禁 | [状态](reports/polyhedral-fvm/2026-09-22-status.md)、[性能](reports/polyhedral-fvm/2026-09-22-performance-gate.md)、[并行](reports/polyhedral-fvm/2026-09-23-parallel-consistency.md) | [验证](../../../docs/developer-guide/verification.md) |
| 方腔与 Reynolds 数研究 | [Ghia](reports/cavity-ghia-validation.md)、[原 PDF](reports/cavity-ghia-validation.pdf)、[Re=10000](reports/cavity-re10000.md) | [求解器](../../../docs/user-guide/models.md)；PDF 保留原版 |
| 旧后端与线性配置比较 | [后端](reports/backend-performance-optimization.md)、[线性基准](reports/re1000-linear-solver-benchmark.md)、[n128](reports/re1000-n128-mpi-performance.md)、[n256](reports/re1000-n256-linear-solver-mpi.md) | [当前配置](../../../docs/reference/numerical-controls.md)、[测量](../../../docs/developer-guide/performance.md) |
| SIMPLE 串并行 | [2026-09-10](reports/simple-parallel-consistency-2026-09-10.md) | `make test-simple-parallel` |
| RANS 方程核对 | [模型报告](reports/rans-equation-verification.md) | [模型](../../../docs/user-guide/models.md)、`make test-rans` |
| 物理审计与修复 | [审计](reports/physics-audit-2026-09-26.md)、[修复](reports/physics-fixes-2026-09-26.md) | [验证](../../../docs/developer-guide/verification.md) |
| MPI 审计与扩展测量 | [审计](reports/mpi-audit-2026-09-26.md)、[测量](reports/mpi-scaling-2026-09-26.md) | [工程规模报告](../reports/flow-optimization-2026-10-07.md) |
| ROCm 可行性 | [调研](reports/rocm-feasibility-2026-09-26.md) | 当前生产证据为 CPU PETSc/MPI |
| 2026-09-28 扩展实验 | [分析](reports/strong-scaling-analysis-2026-09-28.md)、[CSV](reports/strong-scaling-analysis-2026-09-28.csv)、[实验](reports/solver-scaling-experiments-2026-09-28/README.md) | [当前测量方法](../../../docs/developer-guide/performance.md) |

相应 `*-evidence/`、实验目录及图像与报告一起保留。历史命令中的旧 API、绝对路径、
临时目录和分支反映当时环境，不能直接当作现行复现命令。完整旧 CFD 场可能未纳入仓库；
报告存在不表示全部原运行产物仍可获取，结论应按其版本、输入、状态和日志理解。
