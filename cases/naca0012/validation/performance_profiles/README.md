# 性能短测原始数据

`summary.json` 汇总了 2/4/8 rank 扩展、非正交修正次数对照、压力方程 ILUT/AMG 对照、PISO corrector 对照，以及全部/仅压力方程关闭预条件器的失败单步试验。`fullrun_dt0.05_100steps/` 和 `fullrun_dt0.5_100steps/` 分别保存两组完整 100 步运行的 rank 计数、输入配置、运行日志与墙钟记录。

每个子目录保存该试验的 rank 性能 JSON、运行日志、运行摘要和 `solution.bs` 配置副本。短测均使用独立临时算例、同一 127,013 cells 网格、相同初始场和 `deltaT=0.005`；成功组执行 5 步。完整运行分别用 `deltaT=0.05` 与 `0.5`，覆盖物理时间 5 s 与 50 s，不将它们与短测混作同一时间步的速度 A/B，也不替代时间步精度或物理验证。

性能报告中的计数器单位和字段定义见 BabelSim runtime 的 `-performance` 报告；多 rank 子目录保留各 rank 的 owned/ghost/communication-face 数据，可检查负载分区差异。
