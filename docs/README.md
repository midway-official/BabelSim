# BabelSim 产品文档

[项目首页](../README.md) · [用户手册](user-guide/README.md) · [配置参考](reference/README.md) · [开发指南](developer-guide/README.md) · [可运行示例](examples/README.md)

本目录是随 BabelSim 源码交付的现行产品文档，按使用、查询和开发任务组织。内容与当前源码一起维护；本轮接口核对日期为 2026-10-10。终端命令默认在项目根目录执行。

## 选择阅读路线

| 你的任务 | 建议路线 |
| --- | --- |
| 第一次安装与使用 | [安装](user-guide/installation.md) → [首次运行](user-guide/quick-start.md) → [结果处理](user-guide/results-and-restart.md) |
| 建立自己的计算 | [算例管理](user-guide/case-management.md) → [网格/边界](user-guide/mesh-and-boundaries.md) → [模型](user-guide/models.md) → [数值](user-guide/numerics.md) → [运行](user-guide/running-and-parallel.md) |
| 编写或接入新求解器 | [架构](developer-guide/architecture.md) → [公开 API](developer-guide/public-api.md) → [完整开发示例](developer-guide/extending-solvers.md) → [验证](developer-guide/verification.md) |
| 准备软件交付 | [安装](user-guide/installation.md) → [交付与验收](user-guide/delivery.md) → [文档维护](developer-guide/documentation.md) |

## 目录结构

```text
docs/
  README.md                  产品文档入口
  user-guide/                完整用户手册，按 01—12 连续阅读
  reference/                 命令行、字典、模型、数值与文件格式
  developer-guide/           架构、API、扩展、验证和性能维护
  examples/                  完整源码、输入生成器与自动验证
```

历史计划、审计、研究报告和冻结证据已单独归档到 `archive/docs/`（可选归档），交付手册不依赖这些旧记录。

## 用户手册

| 章节 | 主题 |
| --- | --- |
| 01 | [产品与工作流程](user-guide/overview.md) |
| 02 | [安装与一键构建](user-guide/installation.md) |
| 03 | [首次运行：热传导](user-guide/quick-start.md) |
| 04 | [算例管理](user-guide/case-management.md) |
| 05 | [网格与边界条件](user-guide/mesh-and-boundaries.md) |
| 06 | [物理模型](user-guide/models.md) |
| 07 | [数值设置与收敛](user-guide/numerics.md) |
| 08 | [运行与并行计算](user-guide/running-and-parallel.md) |
| 09 | [结果与继续计算](user-guide/results-and-restart.md) |
| 10 | [详细操作示例](user-guide/tutorials.md) |
| 11 | [故障排查](user-guide/troubleshooting.md) |
| 12 | [交付、维护与验收](user-guide/delivery.md) |

## 查询与开发

参数查询集中在[参考目录](reference/README.md)：[命令行](reference/command-line.md)、[算例字典](reference/dictionaries.md)、[物理参数](reference/physical-models.md)、[数值配置](reference/numerical-controls.md)、[文件格式](reference/file-formats.md)。

开发内容集中在[开发指南](developer-guide/README.md)：[架构](developer-guide/architecture.md)、[公开 API](developer-guide/public-api.md)、[新增求解器](developer-guide/extending-solvers.md)、[验证](developer-guide/verification.md)、[性能测量](developer-guide/performance.md)、[分区模型](developer-guide/partitioning.md)、[文档维护](developer-guide/documentation.md)。

[示例目录](examples/README.md)提供完整可执行文件：小型热/输运/流动算例生成、公开接口热求解器、双场耦合、独立 SDK 工程、单单元原生网格与自动验收。

## 文档使用约定

操作教程提供直接可复制命令与预期产物；参考表给出准确参数与约束；开发指南解释数学和生命周期。示例通过表示其软件执行范围通过，物理精度仍需对应研究。模型与重启限制分别在相关章节明确说明。

接口以当前[公开头](../include/babelsim/)和实现为依据，文档检查覆盖链接、章节入口与 Make 目标，执行检查覆盖真实教程。依赖来源与安装脚本见[dependencies](../dependencies/README.md)，项目许可见[LICENSE](../LICENSE)。
