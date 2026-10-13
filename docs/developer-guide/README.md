# BabelSim 开发指南

[文档首页](../README.md) · [用户手册](../user-guide/README.md) · [配置参考](../reference/README.md)

本指南面向用 BabelSim 编写新方程、耦合算法或独立计算程序的开发者。建议先完成用户手册的热传导教程，再编译公开接口示例。框架维护者还应阅读架构与验证章节。

| 阅读顺序 | 内容 |
| --- | --- |
| 01 [架构与维护边界](architecture.md) | 模块职责、数据所有权、公开接口与后端生命周期 |
| 02 [公开 C++ 接口](public-api.md) | Case、Field、math、equ、时间历史、诊断和错误状态 |
| 03 [新增与调用求解器](extending-solvers.md) | 完整外部热求解器、内置注册与双场耦合示例 |
| 04 [验证与回归](verification.md) | 测试目标、失败路径、并行一致性和物理证据范围 |
| 05 [文档维护](documentation.md) | 源码依据、目录规则、迁移关系和示例检查 |
| 高级 [时间推进与关联历史](time-integration.md) | 公共时间权重、预测、历史同步与可选步骤回滚 |
| 高级 [模块描述与模板](module-descriptions.md) | 单份字段/参数声明、前端检查、生成文档与算例模板 |
| 记录 [模块描述与时间工具验证](dsl-time-validation.md) | 源码精简、代码量、回归与性能对照 |
| 维护 [框架源码阅读与修改](framework-maintenance.md) | 非 Physics 实现的调用链、数学符号、MPI 数据流和维护规则 |
| 高级 [性能测量](performance.md) | 固定工作量、计时、并行与内存证据 |
| 高级 [图分区模型](partitioning.md) | 合法网格图、分区候选与负载比较 |

[示例目录](../examples/README.md)包含直接构建的 C++ 源文件、生成完整算例的 Python 脚本和自动验证入口。代码不是延迟表达式配置：开发者通过普通 C++ 循环明确决定装配、求解、耦合、接受及写出时机。
