# 文档整理与交付核对：2026-10-08

按当前公开头、IO、Physics、PETSc/MPI 和 Makefile 更新 11 份现行入口/手册，其中新增
构建页和网格/场/结果页。早期开发文件移入 `docs/history/`，只作为历史记录。
现行入口为 [文档中心](../../../docs/README.md)，迁移范围见 [历史索引](../history/README.md) 和
[31 项迁移清单](../history/archive-map.json)。

## 修正内容

- 根 README 与文档中心：更新 PETSc KSP/PC 支持、五个注册求解器、依赖安装路径和导航。
- 构建页：PETSC_DIR、可选 METIS、独立 BUILD、快速数学、本机指令、调试与外部库链接。
- DSL 与架构：声明/验证/启动生命周期、Physics 边界、长期方程及 Mat/KSP/PC、矩阵和 RHS 分离。
- 数值页：方程/项/算子层级、真实残差、复用与重试，以及 PISO 总压力 p 与 SIMPLE pPrime 的区别。
- 求解器页：PISO 动量松弛默认 1、耦合与湍流接受、NACA 当前配置；删除固定旧日志数字的承诺。
- 文件页：v3 面拓扑、初场和重启历史、命名/未命名序列、最近快照与后处理。
- 验证页：精确列出 test-full 顺序，明确 PISO 专用组另行执行，区分数值与物理证据。
- 性能页：48 个工程测点、同版串行分母、有效循环窗口、全生命周期 RSS/PSS 和限制。

归档保留旧数字和证据，添加历史标记、修正链接。295 个原已跟踪非 Markdown 文件
逐字节不变，冻结实验源码/清单保持原样。历史源码位置只作为对应时期的引用，不冒充
当前行号。NACA 文档的已不可用大场路径改为明确的历史路径说明。

生成 CSV 的 CRLF、Matplotlib SVG 的路径格式和原始日志保持原字节；
`.gitattributes` 仅将这些生成证据的格式排除出空白风格检查，不启用换行转换。
源码、手册和分析脚本仍执行正常检查。

## 实际检查

| 检查 | 结果与证据 |
| --- | --- |
| `python3 tools/check_docs.py` | 现行及历史文件链接、锚点、Make 目标；[诊断](documentation-refresh-2026-10-08-evidence/docs-check.json) |
| `python3 tests/architecture_test.py` | 78 个源码/头的依赖闭包和分层通过 |
| DSL 最小 Heat 示例 | 语法、链接、1/2 MPI 求解和字段对照通过；[记录](documentation-refresh-2026-10-08-evidence/minimal-heat-validation.json)、[源码](documentation-refresh-2026-10-08-evidence/minimal-heat.cpp) |
| `make test-external` | 外部公开接口、1/2/4 rank、负向 API 与错误路径通过；[日志](documentation-refresh-2026-10-08-evidence/external-solver.log) |
| `make test-piso-coupling` | 保守但未耦合场的拒绝、适应校正、1/2/4 rank 通过；[日志](documentation-refresh-2026-10-08-evidence/piso-coupling.log) |
| 当前求解源码和二进制 | 25 个冻结生产/测试修改文件与已通过 full 回归快照相符，当前生产树较检查点无改动 |
| 48 个工程测点 | 重新读取原 rank JSON、工作量和状态，重算端点中位数及自身串行效率分母；[审计](documentation-refresh-2026-10-08-evidence/source-and-measurement-audit.json) |
| 性能图 | 人工检查 PISO/SIMPLE 强扩展及 PSS 图，图例、工作量和指标范围可读 |

检查和数值任务按顺序执行。完整 full 回归是前一轮生产源码的实际执行结果，其日志哈希
与当前生产源码匹配；本次文档更新没有把它改称为新执行的完整回归。
工程报告的 Git 检查点反映当时状态，文档推送后的 HEAD 以 Git 和远端为准。

实际性能结论见 [工程优化报告](flow-optimization-2026-10-07.md)：PISO/SIMPLE 计算速度
明显提高，12 rank 自身效率改善仍有限，SIMPLE 8 rank 效率没有改善；不承诺理想线性
或多节点扩展。工程 SIMPLE 采用固定工作量，尚未达到稳态判据，不能当作稳态精度证据。
