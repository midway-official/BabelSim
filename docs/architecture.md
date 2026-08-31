# BabelSim 架构：作者与维护者的两条入口

## 本轮修正了什么

旧代码的主要问题不是数学表达式不够漂亮，而是每个新 Solver 还要重复编写
Case reader、Field/RunTime 构造、并行输出和日志；SIMPLE 的公开头还列出了全部算法场和工作区。
这使阅读一个 PDE 之前必须先理解许多实现对象。

本轮把重复执行工作收回框架：

- `Case` 统一命名配置、读场、对象生存期和结果工作流；
- `solver.h` 只提供 solve、fvc 和 diagnostics，不包含 RunTime/MPI/代数实现头；
- SIMPLE 公开头只声明算法操作，状态在 `src/internal/simple_state.h`；
- 删除三个专用 Case reader 和三个 Case 运行适配器；
- 新 Solver 一个 `src/physics/<name>/main.cpp` 加显式选择表的一项；
- 时间步自动输出；PVD 使用真实物理时间并引用 XML VTU；
- 时间历史与 solve 次数解耦，允许一个时间步内进行耦合迭代。

没有引入 Solver 基类、Factory、Registry 或通用算法引擎。只有 Case 和 Parameters 两个新的
职责对象：前者是输入/运行结果的所有者，后者是带诊断的命名配置，不是数学或物理管理器。

## 1. 面向 Solver 作者的架构

```text
Case：命名场、物性、方法/时间/输出设置
  └── 普通 Solver 函数或已有 Algorithm
        ├── Field / Boundary / 物性
        ├── fvm / fvc / Equation / solve
        └── 数学收敛与守恒 diagnostics
```

作者描述求解什么、方程如何耦合、何时收敛。普通源码只需要 case.h、solver.h；
使用现成 SIMPLE 再包含 simple.h。无需知道运行对象、分区、消息、矩阵和内部存储。

学习流程和可复制示例见 [solver-development.md](solver-development.md)。
可见 Field 仍是统一的数学/存储类型，而不是另一套包装对象；底层数据接口尚未由 C++ 访问控制
完全隔离，但它们不是 Solver 编程契约。源码和头文件闭包测试防止主 Solver 重新引入实现依赖。

## 2. 面向框架维护者的实际依赖

```text
应用：参数 / MPI 生命周期 / 显式 solver_selection
  └── Case（IO 与对象所有权） ────────┐
        └── Physics / Algorithm       │
              └── solver.h / fvm / fvc│
                    └── RunTime ◀─────┘
                          ├── 离散 / 装配 / DiscreteEquation
                          ├── 线性代数
                          └── Parallel / Halo / 全局归约
Mesh / Field / Boundary 为上述层共享的基础数学对象，不反向依赖物理或 MPI。
```

这是组合依赖，不是 Manager 调用链。Case 不计算梯度，不组装方程，不判断 SIMPLE 收敛；
RunTime 不读热物性，不认识 heat/simple/transport；代数不认识具体物理变量。

| 模块 | 可以知道 | 不应该知道 |
| --- | --- | --- |
| apps/solver_selection | Solver 函数名称 | 方程、Field 数据、工作区 |
| Case / IO | 路径、命名配置、Field 生命周期、时间输出 | SIMPLE 内迭代和离散系数 |
| Physics / Algorithm | PDE、物性、边界、数学修正和收敛 | MPI、CSR/LDU、Eigen、原始场存储 |
| fvm/fvc 描述 | 运算种类、常数、Field 引用、数学正负号 | 具体矩阵分配或通信 |
| RunTime | 当前时间、历史场、方法、数学操作的执行 | 物理模型或 Solver 名称 |
| 通用离散 | Mesh、Field、BC、Method、离散方程 | SIMPLE/热模型 |
| SIMPLE 专用离散 | rAU、压力参考、动量插值 | Case 路径、输出策略 |
| Algebra | 数值矩阵、向量、范数、预条件器 | 温度/压力的物理含义 |
| Parallel | 局部分布、halo、collective、MPI 生命周期 | PDE 或外迭代停止规则 |

`simple_discretization.cpp` 是明确的领域专用数值实现，不应被理解为通用
梯度/扩散核。它复用通用场、离散和线性后端，并没有第二套通信结构。

## 3. 所有权和初始化契约

Case 的私有实现依次拥有配置、并行绑定、局部 Mesh、命名 Field 和 RunTime。
析构顺序反向：RunTime/历史/后端先释放，源 Field 再释放，Mesh 最后释放。
Field 实体地址稳定，容器扩展不会移动已返回的 Field。命名中间场也由同一个 Case 拥有。

Solver 的引用不得逃出 Case 的生存期。所有物理/中间场应在循环开始前声明；同名同类型再次读取
返回已有对象，错误地重用类型/位置会被拒绝。程序正常运行时一个线程只能有一个活动 RunTime；
当前不支持多个同时活动的 Case，多区域需另行设计，不能用全局切换偷偷模拟。

低层内存型数值测试仍可显式构造 Mesh、Field、RunTime，并用
`simple_control.h` 的控制与 iterate API。这个入口服务框架维护和嵌入，
不是新 Solver 必须学习的编程模型。调用者此时负责让源场覆盖整个时间积分生命周期。

SIMPLE 保留一个固定分配的私有 State：

| 状态 | 内容 | 所属职责 |
| --- | --- | --- |
| 物理场引用 | U、p、phi | Case/物理对象 |
| 算法场 | pPrime、rAU、phiHbyA、UPrevious | 压力速度耦合 |
| 数值工作场 | gradP、rAUGradP、面插值、散度 | 避免每次外迭代分配大数组 |
| 小型控制状态 | 当前步骤、迭代次数、各类收敛结果 | 防止漏步/错序；不含 MPI 状态 |

公开 SimpleSolver 只持有一个内部所有权句柄。新增此句柄的成本是构造时一次分配；
外迭代仍复用全部工作场。隐藏状态不是删除状态，更不是每步重新构造状态。

## 4. 数学表达与离散方程

`lhs == rhs` 形成轻量描述符列表；不会在加减表达式时复制 LDU 或大型 Field。
只有 solve 才解释这些项，选择 Method、处理边界并生成局部 owned 行。
ScalarEquation/VectorEquation 是内部离散系数表示，不是作者看到的数学 Equation。

fvm::source 是显式已知体源的方程贡献，不是对源 Field 的隐式线性化。
fvm::source(a,C) 直接复用描述符的系数；不用创建 a*C 的中间大场。
详见 [fvm-fvc.md](fvm-fvc.md)。

底层已有稀疏结构、线性工作区和几何缓存继续复用。表达式列表仍有小型容器成本，
不能声称零分配；本轮没有新增每次 solve 的大场复制。历史场复制从每次 solve
降为每个时间步一次，多次内迭代尤其受益。

## 5. 时间与并行契约

- RunTime 负责时间推进；Case 只调用它并安排写出，不另建 TimeStepper。
- 同一物理时间步重复 solve 不改变旧时间层；进入下一步才推进历史。
- Euler 可缩短末步精确到 endTime；BDF2 首步 Euler，随后等步长 BDF2。
- 当前不支持非均匀步长 BDF2，配置时拒绝非整数步数区间。
- MPI 由应用初始化/结束；Case 和 RunTime 检查 initialized/finalized 后绑定。
- 所有线性残差、场变化和守恒诊断按全局值控制行为。
- Krylov 若在有限精度下提前停止但真实残差未达标，可在原总迭代预算内用真实残差重启。
  不放宽容差，也不通过增加 SIMPLE 外迭代来补偿。
- Case 写出复用 owned-only 并行 writer；独立后处理才聚集全局结果。

## 6. 数据结构与扩展边界

Mesh 仍是三维结构化六面体：二维取 nz=1，整数索引描述 owner/neighbour、cell-face、
patch、owned/ghost/global-ID；连续数组缓存几何与非正交量。Field<T> 共用连续数组，
支持 scalar/vector/tensor，容量由 Mesh/位置确定，不允许任意 resize。

这些实现足以让新方程复用现有 FVM，不等于已经支持任意 PDE、非结构网格、FDM/FEM、
自动单位检查、restart、一般材料模型或隐式块耦合。未来新增后端应在数学描述到离散系统之间扩展，
不要让 Solver 自己选择 LDU、halo 或矩阵类型。

维护者修改后至少运行：make test、make test-workflow、make test-mpi 和
make test-mpi-poiseuille。二次开发测试不仅看代码长度，还检查公共头文件依赖、
新双场算法、多个时间层、失败路径和真实 ParaView 读取。
