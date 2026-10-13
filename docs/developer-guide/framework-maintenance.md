# 框架源码阅读与修改

[开发指南](README.md) · [架构边界](architecture.md) · [验证入口](verification.md)

本文面向掌握 C++17、MPI 和基本数值方法、准备修改 `src/physics` 以外实现的开发者。
先沿一条真实调用链阅读，再进入具体算法。普通函数和数据结构应足以说明流程；
不要为了缩短代码增加运行时注册器、表达式解释器或多层策略类。

## 从一条方程走到结果

以标量扩散方程为例：

1. [application.cpp](../../src/runtime/application.cpp) 管理应用入口和 MPI 会话，选择已注册的求解器。
2. [case.cpp](../../src/io/case.cpp) 读取配置、拥有网格和命名场，完成初始化检查。
3. [procedural_equation.cpp](../../src/discretization/procedural_equation.cpp) 将 `equ::ddt`、`equ::laplacian`、源项等立即加入 `DiscreteEquation`。
4. [petsc_backend.cpp](../../src/backend/petsc_backend.cpp) 根据方程身份取得长期存在的线性系统，处理同步、装配和求解。
5. [petsc_assembly_plan.cpp](../../src/backend/petsc_assembly_plan.cpp) 把面系数映射到预先确定的 COO 位置；[petsc_index_map.cpp](../../src/backend/petsc_index_map.cpp) 负责局部网格编号到 PETSc 自由度编号的转换。
6. [petsc_linear_system.cpp](../../src/backend/petsc_linear_system.cpp) 更新矩阵与右端、控制预条件器复用、求解并检查真实残差。
7. 求解值回到原来绑定的场。显式场运算经 [fvm_execution.cpp](../../src/discretization/fvm_execution.cpp) 校验输入并同步；几何离散公式在 [operators.cpp](../../src/discretization/operators.cpp)。

`Equation<T>` 的模板实现留在 `.cpp`，文件末尾只为 `double` 和 `Vec3` 显式实例化。
这不是动态类型分派；`if constexpr` 在编译时选择标量或向量分支。
`Storage` 隐藏装配状态，`EquationAccess` 仅提供内部访问，公开接口不用知道 PETSc 类型。

## 修改应落在哪里

| 要修改的内容 | 首先阅读 | 需要保留的约定 |
| --- | --- | --- |
| 几何与连接 | [mesh.cpp](../../src/core/mesh.cpp)、[mesh_access.h](../../src/internal/mesh_access.h) | 面方向、owner/neighbour、体积及局部/全局编号 |
| 梯度、插值、对流、扩散 | [operators.cpp](../../src/discretization/operators.cpp) | 共享面的符号、边界贡献、显式修正 |
| `equ` 装配和矩阵操作 | [procedural_equation.cpp](../../src/discretization/procedural_equation.cpp) | 体积分量、绑定未知场、冻结通量、时间历史 |
| `math` 执行及全局统计 | [fvm_execution.cpp](../../src/discretization/fvm_execution.cpp)、[field_math.cpp](../../src/discretization/field_math.cpp) | 输入输出位置、同步顺序、只归约 owned 数据 |
| MPI 基本归约 | [parallel_context.cpp](../../src/parallel/parallel_context.cpp) | communicator、所有 rank 的 collective 顺序 |
| 分区读取与分发 | [distributed_mesh.cpp](../../src/parallel/distributed_mesh.cpp) | root 上的全局网格、局部分区打包格式、集体错误传播 |
| 图划分策略 | [mesh_partition.cpp](../../src/parallel/mesh_partition.cpp) | 图邻接与负载指标；不能改变几何定义 |
| 重复 halo 交换 | [halo_exchange.cpp](../../src/parallel/halo_exchange.cpp) | 已建计划、缓冲复用、begin/finish 对应关系 |
| 线性求解性能 | [petsc_linear_system.cpp](../../src/backend/petsc_linear_system.cpp) | 当前矩阵与复用 PC 分离；真实残差决定接受 |
| 数值配置解析 | [numerics_reader.cpp](../../src/io/numerics_reader.cpp)、[config.cpp](../../src/io/config.cpp) | 完整校验、清晰错误位置、命名选项优先级 |
| 模块声明与生成 | [model_description.cpp](../../src/io/model_description.cpp) | 字段/参数声明、初始化检查、模板和描述输出一致 |
| 结果读写 | [result_reader.cpp](../../src/io/result_reader.cpp)、[parallel_writer.cpp](../../src/parallel/parallel_writer.cpp)、[postprocess.cpp](../../src/io/postprocess.cpp) | 全局 ID、分区元数据、owned 输出 |

## 数学符号与存储

`DiscreteEquation` 保存的是积分方程 `A*x = b`。`diagonal`、`source` 按单元索引，
`upper`、`lower` 按面索引。内部面 `f` 的 owner 为 `o`、neighbour 为 `n` 时：

```text
(A*x)[o] += upper[f] * x[n]
(A*x)[n] += lower[f] * x[o]
```

一个 rank 可以读取 ghost 值，但只更新属于自己的矩阵行。`owner/neighbour` 表示几何方向，
`owned/ghost` 表示进程归属，两组概念不同。边界面的 neighbour 是 `invalid_index`。

底层 `addDiffusion` 装配 `-div(k grad(x))`；公开 `equ::laplacian` 表示
`+div(k grad(x))`，因此加入矩阵时乘 `-factor`。不能按函数名字猜符号。
`equ::source`、`reaction`、时间项会乘体积；`addRhs`、`addDiagonal` 接收已积分的量。

时间项采用 `V * capacity * (a0*x - a1*previous - a2*older)`：`a0` 进入对角，
已知时间层进入右端。BDF2 中 `a2` 为负；步长比和初始 Euler 步由既有时间契约决定。
装配函数只读历史，不能推进或覆盖历史场。

标量扩散装配保存系数、边界和实际显式通量修正快照。`faceFlux` 用新解计算正交部分，
再加装配时冻结的修正。限制器依赖旧场，不能在解更新后重新计算限制器来替代快照。
方程复制共享不可变快照；方程相加先复制贡献列表，以支持自加而不使迭代器失效。

## MPI 数据流和缓存

`distributed_mesh.cpp` 中 `PackedPartition` 的计数槽有名称，打包和解包必须遵守相同顺序。
拓扑数组依次包含面顶点偏移、顶点 ID、owner/neighbour/patch、面 ID/rank、单元 ID/rank/depth。
修改协议时必须同时检查两端，以及 `collectiveMeshError` 的调用位置，避免单 rank 抛错后其余 rank 卡住。

`HaloExchange` 构造阶段发送需要的全局 ID，并将收到的请求转换为本地索引。
计划的 `receive_counts` 是本 rank 向邻居请求的实体数；`send_counts` 是邻居向本 rank 请求的实体数。
稳态交换阶段只做打包、非阻塞通信、等待和解包。分量布局与缓冲区可复用；
标量、向量、张量分别使用 1、3、9 个连续 double。
不要把全局 ID 查找、全局网格遍历或重复分配移回交换热循环。

PETSc 稀疏模式、矩阵数值和预条件器是三个不同层次。
网格/编号/隐式 stencil 决定模式；当前系数决定矩阵；PC 可以在既有预算内复用，
但残差始终针对当前 `b-A*x`。结构拆分不能放宽容差、减少迭代工作或改变复用策略来补偿性能。

## 代码写法

- 一条语句表达一个动作。校验、同步、局部计算、结果同步之间保留清晰边界。
- 用 `field`、`equation_state`、`face_diffusivity` 等名字说明角色；局部数学公式可用 `a0/a1/a2`，并解释公式。
- 同一数学规则的标量/向量分量可以共用小函数；不同同步需求不能用复杂通用回调掩盖。
- 注释解释符号、所有权、冻结状态和复用条件，避免逐行翻译 C++。
- 添加全局指标时更新 `maximumPerformance` 中相应成员表；打包和还原使用同一份表，仍执行一次批量 MPI 归约。
- `src/.clang-format` 固定四空格、100 列、不重排 include。使用 clang-format 14；只格式化本次修改范围。

对框架做一次完整格式检查时可使用：

```bash
find src -path src/physics -prune -o \( -name '*.cpp' -o -name '*.h' \) -print0 \
  | xargs -0 clang-format --dry-run --Werror
```

格式化只是阅读辅助。拆职责、去掉重复的数值规则、保留可定位的失败状态才是维护工作的核心。

## 验证修改

先执行 `make test-architecture`，再按[验证指南](verification.md)选择测试。
跨离散/MPI/后端的修改至少覆盖单元测试、`test-mpi`、`petsc-check` 和相关模型工作流。
`test-full` 不包含 `test-piso-coupling`，涉及通量、边界或装配时需单独执行。

性能对照保留修改前后程序，使用相同编译参数、输入、rank 绑定和数值设置，交替重复运行。
同时比对最终场、接受步数、线性迭代与真实残差。编译通过、数值测试通过、性能通过分别记录。
短算例只能反映对应路径，不能据此声称任意规模、任意机器都无性能损失。
