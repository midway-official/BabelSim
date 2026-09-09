# 本轮验收证据

本目录对应 F1–F7 修复后工作区；不覆盖或改写 framework-review-2026-09-09-evidence。

- acceptance.log：最终完整 Make 目标日志（并行执行的目标输出可能交错），整体退出 0。
- physical-validation.log：已收敛的方腔和 Poiseuille，整体退出 0。
- numerical-contract-*.log：最后新增的向量 Sp/原方程残差断言与 30 组预条件/尺度测试。
- rans-summary.json、rans-logs/：最终模型验收脚本的全部运行日志；case 由 tests/rans_validation_test.py 确定生成，临时目录位置写在 summary 中。
- fixed-*.log：使用旧审查探针的原始源文件，链接当前框架重新执行。
- channel-old-tolerance.log / channel-refined.log：新原方程残差拒绝旧精度及收紧内层容差后通过的证据。
- source-sha256.txt：当前受测源码、测试、工具及仓库 case 输入的 SHA-256；不包含生成结果。
- status.json：已核对终端退出状态。

复现主命令见上级 f1-f7-implementation.md。独立旧探针源位于 ../framework-review-2026-09-09-evidence/{scale,inletoutlet,transient}.cpp。
它们访问维护级 Eigen/RunTime 接口，需使用项目相同编译配置，包括 -march=native、-DNDEBUG，避免 Eigen 对齐 ABI 不匹配；普通 Solver SDK 外部构建另由 test-external 验证。

本次第一次用 -O2 且未匹配 -march=native/-DNDEBUG 编译 Eigen 维护探针，发生对齐 ABI 相关的内存释放异常；按项目 Makefile 配置重编后三个探针全部通过。
这不是可忽略的通用 ABI 保证：维护级接口跨编译选项混用仍须避免。仓库日志中保留的正确运行使用与库相同配置。

负向测试会有 healthy=0、converged=false 或 mpirun 退出 2；脚本验证这些预期结果后整体成功。
文件名中的 fixed 表示修复后运行，不表示编译器或数值格式被手工改成更容易通过的路径。
