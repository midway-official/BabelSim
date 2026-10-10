# BabelSim 用户手册

[文档首页](../README.md) · [配置参考](../reference/README.md) · [开发指南](../developer-guide/README.md)

本手册面向安装、配置和运行 BabelSim 的使用者。第一次使用时，按章节顺序完成热传导教程；准备自己的工程计算时，重点阅读网格、边界、物理模型、数值设置和结果检查。所有终端命令默认在 BabelSim 源码根目录执行。

| 章节 | 解决的问题 |
| --- | --- |
| 01 [产品与工作流程](overview.md) | BabelSim 的用途、能力和适用范围是什么？ |
| 02 [安装与构建](installation.md) | 如何检查环境、安装依赖并用 make 编译？ |
| 03 [首次运行](quick-start.md) | 如何完成第一个计算并在 ParaView 中打开结果？ |
| 04 [算例管理](case-management.md) | 输入文件放在哪里，修改哪个文件，怎样保留计算记录？ |
| 05 [网格与边界条件](mesh-and-boundaries.md) | 如何生成或导入网格，设置场和边界？ |
| 06 [选择物理模型](models.md) | 热传导、标量输运、SIMPLE、PISO 和 RANS 如何选择？ |
| 07 [数值设置与收敛](numerics.md) | 如何选择离散格式、时间步和线性求解器？ |
| 08 [运行与并行计算](running-and-parallel.md) | 如何启动串行/MPI、保存日志并判断运行状态？ |
| 09 [结果处理与继续计算](results-and-restart.md) | 如何导出结果、选择时刻及从已有场继续？ |
| 10 [详细操作示例](tutorials.md) | 如何完成输运、通道流、PISO 和求解器调用示例？ |
| 11 [故障排查](troubleshooting.md) | 输入错误、构建失败或未收敛怎样定位？ |
| 12 [交付、维护与验收](delivery.md) | 如何打包源码、重建软件并记录验收证据？ |

手册中的参数示例是可运行的入门配置。实际工程计算还需根据几何、流动尺度和关注量完成网格、时间步及模型验证。可查询的完整字段定义集中在[配置参考](../reference/README.md)，新增求解器和外部程序链接步骤集中在[开发指南](../developer-guide/README.md)。
