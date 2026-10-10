# 11 故障排查

[用户手册](README.md) · 上一章：[详细示例](tutorials.md) · 下一章：[交付与验收](delivery.md)

## 按失败阶段定位

先保留完整日志与退出码，确认问题发生在环境检查、依赖配置、框架编译、输入读取、线性求解、耦合推进还是结果导出。最后一条明确错误通常比 MPI 启动器随后打印的通用提示更有用。

| 现象/错误 | 检查与处理 |
| --- | --- |
| 安装脚本检查失败 | 用 `--check` 查看缺项；核对编译器、MPI、Python、BLAS/LAPACK，而不只看命令是否存在 |
| apt-get/sudo 无法安装 | 检查系统源和系统安装权限；可由管理员准备包后再执行 make |
| C/C++/Fortran MPI 版本不一致 | 修正 PATH 和 MPI 包装器，使用同一 MPI 安装重新检查 |
| Hypre SHA-256 不符 | 核对交付包与 dependencies/sources.json；恢复正确本地源码包 |
| 缺少 PETSc.pc | 默认省略 PETSC_DIR 使用本地构建；外部安装检查 lib/pkgconfig 路径 |
| Eigen 头文件缺失 | 检查 dependencies/eigen；增加新模块时补齐源码，不能依赖未交付的系统头 |
| 程序缺少动态库 | `ldd` 检查 MPI/PETSc/Hypre/可选 METIS；确认构建目录和安装路径未被搬迁 |
| Illegal instruction | 可能运行在不同 CPU；按安装手册同时重建框架和 PETSc 的可移植版本 |
| unknown BabelSim solver | 核对 case 的大小写与程序实际链接的注册项 |
| usage / needs -case | 使用支持的 `-case`、`-time`、`-performance`，求解器没有 --help |
| 字典 missing/duplicate/unused 项 | 修正文件与键；不同算法的配置块不能混用 |
| incomplete spatial configuration | 每个方程提供四个基础选项；不要只给实际使用的某一项 |
| BDF2 uniform time steps 错误 | 将 `(endTime-startTime)/deltaT` 设置为整数，或采用 Euler |
| field patch 缺失/重复 | 每个输入场配置网格所有物理 patch，名称完全一致 |
| internal file 覆盖不完整 | 核对 global_id 是否覆盖全部单元、是否重复、分量数是否正确 |
| create fields before starting | 在第一次 advance/Case::start 前声明所有 Case 命名场 |
| icc 并行失败 | 改用可用 MPI 预条件器，如 hypre/gamg/bjacobi，并重新检查系统性质 |
| 退出码 2、线性不收敛 | 核对系数/边界相容性、矩阵与 KSP 的匹配、PC、迭代上限和真实残差 |
| SIMPLE 达上限 | 检查 mass、dU、rU、dP 与模型指标，再有针对性调整松弛/配置 |
| PISO 时间步不接受 | 检查线性、mass、rCoupling、校正次数、dt 与湍流耦合设置 |
| 后处理 rank/mesh 不一致 | 使用完整原始快照和匹配网格；不要混合不同运行的数据 |
| Tecplot 不支持网格 | 一般多面体改用 VTK；当前 Tecplot 使用 FEBRICK |
| 找不到预期最后时刻 | 看 writeInterval、物理 endTime、接受状态和实际最新保存时刻 |
| 输出越积越多 | 用新运行名，控制字段与间隔，并预算数据量 |

## 区分线性与外层失败

先观察 linU/linP 或对应线性结果；它们失败时外层指标没有足够意义。线性已成功但 dU/mass/rCoupling 未满足，问题在耦合、边界或时间推进层面。提高线性上限不是处理所有未收敛的通用办法。

如果残差、场值或矩阵出现非法数值，检查物理参数、初始场、网格几何、模型下限和边界。以更小算例、较小时间步或基础格式复现，再定位造成异常的改动。需要排查内存或未定义行为时使用 `make debug`，并保留独立调试日志。

## 后处理与重启区别

后处理校验配套元数据和网格；场重启按全局编号读取值。后者能够读取数据，不代表已验证新旧几何完全相同。复制一个时刻时，应一起保存 rank 数据和原始 case 输入；从 BDF2 单快照继续时会重新建立历史。

## 反馈所需的最小材料

提供可复现的小型 case、完整日志、启动命令、退出码、MPI 进程数、编译器/MPI/PETSc 版本和构建选项；说明预期行为、实际行为及最后成功时刻。若有版本管理，附源码版本；若为源码包，附交付包校验值。

```bash
./dependencies/install.sh --check
g++ --version
mpirun --version
ldd build/babelsim-solve
```

避免只提供最后一张图或一个 converged 标签；它们无法重建配置与算法状态。
