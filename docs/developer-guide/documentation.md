# 05 文档维护规范

[开发指南](README.md) · 上一章：[验证与回归](verification.md) · [文档首页](../README.md)

## 文档层次

| 目录 | 面向读者 | 内容边界 |
| --- | --- | --- |
| user-guide | 使用者 | 按任务解释从安装到结果与交付的全过程 |
| reference | 参数查询者 | 精确语法、枚举、默认值、约束和格式 |
| developer-guide | 扩展作者/维护者 | API、数学契约、架构、开发与验证 |
| examples | 学习与验收 | 完整可编译代码、生成输入、自动执行 |
| developer-guide 性能章节 | 高级使用者 | 性能测量与分区方法 |
| archive/docs/reports | 研究证据读者 | 版本化结论及冻结证据，位于产品 docs 之外 |
| archive/docs/history | 追溯读者 | 历史开发记录，位于产品 docs 之外 |

根 README 是产品入口；docs/README 是导航；用户手册有连续章节，配置参考有集中表格，开发指南提供可执行代码。避免再把安装、字典、API 和实验日志混进一篇大文档。

## 本次重写范围

本次手册按 2026-10-10 的当前源码接口核对。旧平铺文档内容已按下表重新组织并扩写，旧路径已删除；维护阅读链接同步迁移，冻结实验源码副本保持原样。docs 根目录只保留 README 与四个正式目录，不放兼容空壳、研究日志或历史计划。

| 原文档 | 当前位置 |
| --- | --- |
| building.md | [安装](../user-guide/installation.md)、[扩展链接](extending-solvers.md) |
| solvers.md | [模型使用](../user-guide/models.md)、[物理参数](../reference/physical-models.md)、[算法控制](../reference/numerical-controls.md) |
| numerical-configuration.md | [数值使用](../user-guide/numerics.md)、[完整配置](../reference/numerical-controls.md) |
| mesh-and-results.md | [网格与边界](../user-guide/mesh-and-boundaries.md)、[结果/重启](../user-guide/results-and-restart.md)、[文件格式](../reference/file-formats.md) |
| dsl-runtime-manual.md | [公开接口](public-api.md)、[详细开发示例](extending-solvers.md) |
| architecture.md | [架构](architecture.md) |
| validation.md | [验证](verification.md)、[交付验收](../user-guide/delivery.md) |
| performance/README.md、partition-model.md | [性能测量](performance.md)、[分区模型](partitioning.md) |
| reports/、history/ | `archive/docs/`，位于产品 docs 之外 |

## 内容依据

| 契约 | 主要源码 |
| --- | --- |
| 程序参数、注册、状态 | [application.cpp](../../src/runtime/application.cpp)、[postprocess.cpp](../../src/io/postprocess.cpp) |
| case/control/output 与生命周期 | [case_reader.cpp](../../src/io/case_reader.cpp)、[case.cpp](../../src/io/case.cpp)、[runtime.cpp](../../src/runtime/runtime.cpp) |
| 数值配置 | [numerics_reader.cpp](../../src/io/numerics_reader.cpp)、[methods.h](../../include/babelsim/methods.h)、[solver_control.cpp](../../src/backend/solver_control.cpp) |
| 网格/场/结果 | [mesh_reader.cpp](../../src/io/mesh_reader.cpp)、[field_reader.cpp](../../src/io/field_reader.cpp)、[parallel_writer.cpp](../../src/parallel/parallel_writer.cpp)、[result_reader.cpp](../../src/io/result_reader.cpp) |
| 正式算法与模型 | [src/physics](../../src/physics/) |
| 安装与构建 | [Makefile](../../Makefile)、[dependencies](../../dependencies/README.md) |
| 示例运行 | [verify_examples.py](../examples/verify_examples.py) 与实际 JSON/日志 |

不要依据旧报告、其他 CFD 软件同名配置或默认记忆填写当前键。检查公开函数重载与大小写，尤其是逐方程配置、PC/AMG 控制、Field 生命周期和输出目录。

## 写作规则

教程应给出问题、完整输入生成方式、可复制命令、产物、成功条件和错误定位。代码块中的占位符须明确说明，主阅读路径优先提供可以直接执行的命令。初学者流程不包含没有产品意义的内部布局细节；开发指南再解释积分量、边界、并行和后端契约。

各项能力区分已实现、执行验证与物理验证。模型限制、默认接受阈值和重启范围应在相应使用位置讲清楚。章节包含入口导航和相邻章节链接；参考表统一参数正式名称，旧兼容别名集中说明。

## 检查流程

```bash
python3 tools/check_docs.py --json build/docs-refresh/docs-check.json
python3 tools/check_docs.py --delivery-only --json build/docs-refresh/delivery-check.json
python3 docs/examples/verify_examples.py --output build/docs-refresh
git diff --check
```

默认文档检查覆盖根 README、当前手册/参考/开发/示例、依赖说明、归档/报告入口、历史 Markdown 和受版本管理的 case README，并检查每篇现行手册都能从 docs 首页导航到达。--delivery-only 检查交付说明而不要求 Git/可选归档，现行 docs 不通过阅读链接依赖可选归档。检查不递归将冻结源码里的旧文档当成现行手册。历史缺失链接单独报告，不豁免现行导航。

修改配置或 API 后先更新表格和可执行示例，按改动运行[验证矩阵](verification.md)。不通过的命令不能继续标成当前可用。检查脚本验证语法链接和执行结果，源码语义与工程物理仍需要维护者独立核对。
