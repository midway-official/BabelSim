# 配置与接口参考

[文档首页](../README.md) · [用户手册](../user-guide/README.md) · [开发指南](../developer-guide/README.md)

本目录用于查询精确语法、必需项、默认值与约束。操作步骤从用户手册进入；不要将不同求解器的完整参数块混在同一个算例里。

| 参考 | 内容 |
| --- | --- |
| [命令行](command-line.md) | 求解、后处理及示例工具的参数、退出码与选择规则 |
| [算例字典](dictionaries.md) | case/control/output、文本语法和路径 |
| [物理模型参数](physical-models.md) | 热传导、输运、流动与 RANS 的输入和默认系数 |
| [数值配置](numerical-controls.md) | 时间与空间离散、逐方程线性控制、SIMPLE/PISO 控制 |
| [文件格式](file-formats.md) | 网格 3、场文件、非均匀初值与结果格式 |

公开 C++ API 的类型、生命周期、积分量与算子符号在[开发接口说明](../developer-guide/public-api.md)中解释。参考与当前源码一起维护，版本依据及更新流程见[文档维护规范](../developer-guide/documentation.md)。
