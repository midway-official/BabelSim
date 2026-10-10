# 框架优化范围核对

本记录是源码边界核对，不代表性能、数值或物理验收。比较基线为
`83d7b72757f7f9af2db1a01a097389f86c296117`，优化源码检查点为
`78b712ef24704d14fc64cbcba3b4f9f7cbd2c438`。

| 用户约束 | 本次直接检查的源码证据 | 证据范围 |
|---|---|---|
| 保持顶层 SIMPLE/PISO 语义和抽象层次 | `git diff <基线> -- src/physics` 为空；性能计时由 `RunTime::recordIteration` 和通用 Reporter 生命周期记录 | 证明顶层 Physics 源码未修改；不代替数值对照 |
| 配置解析保留在 Case/IO 层 | `src/io/case.cpp` 的 `equationControl` 解析新增 PC 复用参数；`LinearSolverConfig` 传递已解析的配置并校验范围 | 新增参数为后端执行控制，原有容差、迭代上限及 warm-start 默认值未修改 |
| 求解优化保留当前矩阵和真实残差标准 | `PetscLinearSystem::updateMatrix` 集体判断矩阵数值是否相同；变化时调用 `MatSetValuesCOO`；`solveRhs` 每次计算当前矩阵的 `b-Ax`，滞后 PC 失败时恢复原始初值并重建重试 | 复用 PC 不等于复用旧矩阵或放宽残差目标；失败恢复也保持原始目标 |
| 划分与负载模型适用于一般合法网格 | `partitionOwners` 由 `cellNeighbours` 构建 CSR，以 `[1, cellFaces.size()]` 为双约束；连通性由图本身判断；无 METIS 时使用图距离和图生长 | 选择由拓扑和工作量决定；适用范围仍受网格契约、索引宽度与可用内存约束 |
| 共享 FVM/Field 路径优化 | `Field::updateBoundaryTrace` 与内部 field access 遍历已有 patch 面列表；`Mesh::patchFaces` 提供只读面列表视图 | 优化对象是任意合法 patch 的面访问，不依赖测量算例的网格尺寸或名字 |

核对时也搜索了 PETSc、网格划分、算子和并行上下文实现，未发现本次测量算例的名字或
单元总数常量用于选择执行路径。这个搜索只是辅助证据；上表的实际控制逻辑才是通用性的依据。

矩阵更新、PC 生命周期及完整 stencil/halo 的数值行为仍须结合已保存回归日志和工程字段
对照判断。速度、强扩展效率及内存收益须以相同输入、相同工作量的完整测量记录判断，
不能由本记录推导。最终交付审计会另外核验源码快照、计时分母和场数据哈希。
