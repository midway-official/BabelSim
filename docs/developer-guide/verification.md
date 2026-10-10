# 04 验证与回归

[开发指南](README.md) · 上一章：[新增求解器](extending-solvers.md) · 下一章：[文档维护](documentation.md)

## 证据层次

| 证据 | 结论范围 |
| --- | --- |
| 源码/接口核对 | 功能、定义与调用路径存在 |
| 构建与执行 | 指定配置可以编译、链接并启动 |
| 局部契约与失败路径 | 测试覆盖的边界、状态和数值契约成立 |
| 同时刻串并行比较 | 指定网格、版本、容差与进程数下字段一致 |
| 解析/制造解及细化 | 指定问题与细化范围内的误差与阶数 |
| 独立物理参考与敏感性 | 特定物理条件和关注量的适用性 |
| 固定工作量测量 | 指定机器与口径的时间、效率和内存 |

上述结论分别记录。有限场、线性收敛、极小质量误差、一次完成的 pilot 或文件名都不能独立证明工程物理准确。

## Make 测试目标

| 命令 | 内容 |
| --- | --- |
| `make test-architecture` | include 闭包、公开头与 Physics 分层边界 |
| `make petsc-check` | 实际 PETSc/Hypre 后端的 1/2/4 rank 检查 |
| `make test` | 架构与串行单元/契约测试 |
| `make test-workflow` | 新求解器工作流、Case 编程、时间序列、后处理 |
| `make test-external` | 仓库外公开 include 编译/链接与负向 API 检查 |
| `make test-rans` | RANS 方程、模型常数和执行 |
| `make test-scalar-time` | 实际热/输运 Euler/BDF2 时间阶 |
| `make test-simple-parallel` | SIMPLE 系列案例串并行一致性 |
| `make test-mpi` | 分区、halo、数学、契约、后端和各类 MPI 案例 |
| `make test-mpi-heat` | 正式 heat 串并行字段比较 |
| `make test-mpi-poiseuille` | 正式 SIMPLE 1/2/4 rank 与几何导出 |
| `make test-piso-coupling` | PISO 耦合误差、校正与接受判据专用检查 |
| `make validate-cavity` | 方腔参考验证 |
| `make validate-poiseuille` | 通道解析参考验证 |

`make test-full` 顺序执行 test、workflow、external、rans、scalar-time、simple-parallel、mpi、mpi-poiseuille。PISO 专用 test-piso-coupling **不包含**在此汇总中，需要单独执行。`make validate` 组合 test 与方腔/通道入口，不等于完整 test-full。

生产测试与求解任务逐组顺序运行，避免不同组写同一个 case 的结果。构建可以用 -j，不应将多个共享结果目录的测试组并行执行。命令与覆盖以[Makefile](../../Makefile)和[tests](../../tests/)为准。

## 按改动选择检查

| 改动 | 至少核对 |
| --- | --- |
| 文档 | 链接/锚点/目标、实际教程命令、源码参数一致性 |
| 配置/IO | 合法与缺失/重复/未知/未使用键、初值、输出、重启 |
| 公开 API | architecture、external、workflow、正式调用方 |
| 空间离散 | 局部恒等式、边界、非正交、守恒、串并行 |
| 时间格式 | 生产求解器启动、末步、历史、时间阶 |
| 后端/PC | 当前矩阵真实残差、复用/重建/失败重试、1/2/4 rank |
| PISO | 专用 coupling 组及时间步/网格研究 |
| RANS | 模型系数、独立方程、限幅、耦合及物理案例 |
| 图分区/MPI | 合法划分、完整 halo、跨分区通量、owned 输出 |

新算法至少包括一个可解释参考问题，以及线性不收敛、非法输入、外迭代耗尽、适用的非法数值处理、同一时刻串并行和实际保存文件。直接镜像实现的测试不能代替独立参考。

## 手册示例验收

```bash
python3 tools/check_docs.py --json build/docs-refresh/docs-check.json
python3 docs/examples/verify_examples.py --output build/docs-refresh
make test-architecture
```

示例验证使用新建工作目录保留输入与日志，检查真实程序退出码、唯一拥有单元、有限场、时间序列、API 外部链接、注册分派、Euler/BDF2 一致性、跨进程数场继续计算和 VTK 几何。它不会运行所有生产测试，也没有给小型通道/PISO 教程增加工程物理认证。

## 尚需独立关注的契约

向量 symmetry 的隐式扩散切向分量需要专项回归，不能依据标量 symmetry 或平面案例推广任意方向滑移条件。纯 Neumann 参考值处理没有检查源项/边界通量相容性；需用独立兼容与不兼容输入分别验证。保守通量下 mass 很小，也不保证 PISO splitting/coupling 误差已小。

## 报告与性能

每份报告标明源码/二进制、输入、机器、MPI/线程/绑定、退出状态、工作量、计时边界和重复次数。SIMPLE 固定迭代耗尽测吞吐应明确尚未达稳态条件，不能标作收敛时间。强扩展基准须采用相同版本和物理工作量；MPI 字段一致不证明加速。

版本化研究位于可选归档 `archive/docs/reports`，高级口径见[性能指南](performance.md)。冻结源码、原始日志、哈希与测量保持原样；新文档或新测试不能改写旧证据。
