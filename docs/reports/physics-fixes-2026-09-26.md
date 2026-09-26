# 数学物理审查：最小修复与验证（2026-09-26）

对应 [原始审查](physics-audit-2026-09-26.md) 的 F1–F5。原审查及其证据保留修复前状态。
本次保持公开 DSL API、Field 值赋值/边界约束语义、RANS 模型接口及公式，
Physics 继续只通过公开 DSL 组织方程，未引入跨层后端访问。

## 修复

| 项目 | 实现 | 契约 |
|---|---|---|
| F1 | scalar laplacian 保存实际显式面通量；faceFlux 只重算正交项 | 限幅结果与矩阵在同一装配状态冻结；copy/scale/add 保持一致 |
| F2 | vector flux 读取输入场的边界流向；首次以 owner 速度初始化 | 输出缓冲区不决定 inletOutlet 分支；求值后继续更新输入场流向上下文 |
| F3 | 单次 PISO 对模型输运做同时间层内迭代 | 物理历史每步只保存一次；保留默认 0.7 欠松弛及模型接口 |
| F4 | 湍流线性收敛状态加入 PISO 步接受条件 | MaxIterations 或模型未收敛时不写出成功结果、不继续推进 |
| F5 | heat/transport 以变步长二阶外推状态装配显式修正 | BDF2 时间导数仍读取真实历史；Euler 启动、物理边界约束不变 |

新增可选 PISO 设置 `turbulenceMaxIterations`，默认 1000，仅在启用湍流且
`maxIterations=1` 时使用。收敛使用现有 `turbulenceTolerance`；达到上限则返回
`notConverged`。这会增加过去仅求解一次模型的算例成本，避免将未收敛的欠松弛状态
作为真实物理时间步接受。多次外层 PISO 保持每个外层一次模型更新。

外推公式 `q*=(1+r)q_old-r q_older`，其中 `r=dt/dt_previous`。
这是显式修正的二阶时间一致处理，仍有显式稳定性限制。

## 回归证据

- `physics_contract_test`：1/2/4 进程，在扭曲网格上覆盖 orthogonal、corrected、
  limitedCorrected，常/变系数及方程复制、缩放、相加。求解后和再次改变场值后，
  `V*div(faceFlux)-(apply-rhs)` 的 L2 范数均小于 `1e-10`。
  同时验证已有回流上下文、首次流向初始化和污染的输出缓冲区。
- 原有 procedural_equation、numerical_contract、math_runtime、operators、time_history、
  field_boundary 回归通过；parallel_math 在 1/2/4 进程下最大误差为 0；
  numerical_contract 在 2/4 进程下通过。见 [日志](physics-fixes-2026-09-26-evidence/contracts.log)。
- `make test-rans` 对 SA、kOmega、kEpsilon 的公式、Euler/BDF2 衰减与稳态场回归通过。
  新增单次 PISO 与步内收敛的 transient SIMPLE 比较，所有比较差值小于 `2e-8`；
  中间步长覆盖 1/2/4 进程。裁剪导致未收敛与线性迭代上限的 PISO 反例，
  在 1/2/4 进程下均返回退出码 2，未生成成功结果。见 [日志](physics-fixes-2026-09-26-evidence/rans.log)。
- `make test-scalar-time` 调用真实生产入口，固定空间网格，20/40/80/160 步推进至 0.2。
  下表是相邻时间分辨率最终场差值的比值，二阶期望约 4。

| 求解器 | 等步长两组比值 | 缩短末步两组比值 |
|---|---|---|
| heat：corrected diffusion | 4.11685 / 4.04799 | 4.02742 / 4.00508 |
| transport：linearUpwind | 4.04556 / 4.02764 | 3.92358 / 3.96601 |

这两种入口的 2/4 进程结果与串行差值小于 `1e-9`。
见 [标量时间阶日志](physics-fixes-2026-09-26-evidence/scalar-time-order.log)。

`make test-workflow` 通过：heat/transport/coupled、1/2/4 进程、物理时刻、
失败路径，以及 ParaView 的 PVD 时间序列和 VTU 单元值。
架构检查通过（76 个源文件/头文件、无环包含与层间边界），`git diff --check` 通过。
见 [工作流日志](physics-fixes-2026-09-26-evidence/workflow.log)。

验证范围是所修复的离散契约、均匀湍流衰减与标量时间一致性，
不等价于任意湍流流场或完整 PISO 耦合算法的二阶精度证明。
