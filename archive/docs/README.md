# 文档与验证记录归档

[现行产品文档](../../docs/README.md) · [用户手册](../../docs/user-guide/README.md)

本目录保存从产品 docs 目录迁出的历史开发和版本化研究资料。2026-10-10 将交付手册统一为用户、参考、开发和示例四个目录；历史资料仍用于追溯对应版本，不能直接作为现行配置或产品性能承诺。

| 归档 | 内容 |
| --- | --- |
| [history](history/README.md) | 历史计划、接口重构、审计与旧数值/性能报告 |
| [reports](reports/README.md) | 版本化实验报告、读取索引及冻结证据 |
| [旧性能手册](manuals/performance.md)、[旧分区手册](manuals/partitioning.md) | 产品重写前的说明与当时研究数字，保留供追溯 |

原始数据、图像、冻结源码、日志、哈希和机器可读证据保留字节内容。本次迁移仅调整维护中的 Markdown 阅读链接；不改写原始测量、数值结果或冻结源码。迁移核对共验证 6,804 个原始证据文件的 SHA-256 保持一致，本地明细见 build/docs-refresh/migration.json。

这里是项目追溯材料，不是产品手册组成部分；仅交付现行 docs 时，可以独立留存本目录。冻结副本内部旧路径与旧命令反映原运行环境，不自动更新。现行性能使用方法已并入[开发指南](../../docs/developer-guide/performance.md)。
