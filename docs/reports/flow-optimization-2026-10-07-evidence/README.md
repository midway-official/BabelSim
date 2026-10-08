# 通用框架优化证据

实验日期 2026-10-07—08，正文见 [优化报告](../flow-optimization-2026-10-07.md)。
本目录是冻结实验记录，现行用法见 [性能手册](../../performance/README.md)。

## 读取顺序

1. [PISO 汇总](strong-summary.json)、[SIMPLE 汇总](simple-summary.json)：端点三次采样
   中位数、各版本自身串行分母和同 rank 版本对照。
2. 各 `study.json`、`records.json`、rank JSON 和日志：输入、哈希、绑定、工作量、
   退出码、循环时间、内存采样。
3. [字段对照](strong-field-comparisons.json)、[配对清单](strong-field-pairs.json)：
   14 组 PISO 全部 owned 单元、时刻、坐标和字段。
4. [完整回归](tests/test-full-amg-krylov-policy-cf.log)、[交付审计](delivery-source/manifest.json)、
   [审计结果](audit-delivery-launcher.log)：源码匹配、48 次有效测量、分母重算、16 个字段目录复核。
5. [通用性与抽象边界](framework-scope-review.md)、[退出码检查](benchmark-exit-contract-audit.json)、
   [第一次推送记录](checkpoint-git-push.json)。

`final-source/`、`final-source-amg/`、`delivery-source/` 的源文件、README、文档和清单
属于源码快照，本次文档整理不改这些文件。历史审计 Git HEAD 是审计时刻的 HEAD；最新
提交以仓库和远端为准。

计算时间取各 rank 的 `iterationComputeSeconds` 最大值，排除首迭代及此前初始化和
窗口内输出。阶段计数覆盖整个回调。RSS/PSS 峰值覆盖全生命周期，不能称循环内存峰值。
SIMPLE 固定执行 384 次外迭代，状态 `maxIterations`、退出码 `2`，不作为稳态精度或
收敛时间证据。PISO 完成规定物理步数。完整字段保存在实验机清单指定的位置。
