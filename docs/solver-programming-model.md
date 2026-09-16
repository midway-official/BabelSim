# 过程式数值求解器编程模型

求解器是使用公共数值 API 的普通 C++ 函数。作者依次读取场和参数、计算物性、
清空矩阵、逐项离散、求解、修正场、检查收敛和写出结果。
不要求作者定义结构体、继承算法类或拆分公共算法文件。

## 独立的算法程序

- Heat 和 Transport 各自在一个 `main.cpp` 中完成读取、时间循环、组装和输出。
- 稳态 SIMPLE 在自己的 `main.cpp` 中使用一个主要迭代循环。
- 瞬态 SIMPLE 在自己的 `main.cpp` 中使用时间循环和内部 SIMPLE 循环。
- 两个 SIMPLE 不调用共享的 SIMPLE 算法；压力非正交修正可使用显式子循环。
- SA、k-omega、k-epsilon 各自实现输运方程、闭合和历史状态。
  `RANS/api.h` 只提供交互接口与模型选择，不提供共享数值算法。

具体入口见 [稳态 SIMPLE](../src/physics/simple/main.cpp)、
[瞬态 SIMPLE](../src/physics/transient_simple/main.cpp) 和
[Heat](../src/physics/heat/main.cpp)。

## API 的执行语义

`equ::createEquation(field)` 创建绑定未知场的线性系统。`clear` 清空贡献，
`ddt/div/laplacian/reaction/source` 立即组装；`solve` 只求解已组装的矩阵。
`math::grad`、`math::div`、插值和场代数直接产生已计算的场。
组装后的矩阵与已计算的临时场不因输入后续变化而重新求值。

`loadMethods(problem)` 显式加载离散配置。`enableTime(problem)` 创建时间服务；
普通 while 循环调用 `advance(time)` 推进。服务负责末步截断和时间元数据，
不会保存历史、运行物理算法、判断收敛或写文件。
`math::saveOld` 在物理时间步开始时显式保存历史，内迭代不得推进历史。
`equ::ddt` 根据历史和格式配置处理 BDF2 起步与变步长系数。

## 责任边界

| 层 | 负责内容 |
| --- | --- |
| 物理求解器 | 方程、更新次序、循环、松弛、物理收敛与所有运行信息打印 |
| 模型 | 独立闭合与输运，向求解器返回物性和诊断数据 |
| 公共 DSL | 场、矩阵、直接数学运算与历史 |
| 离散层 | 网格几何、边界贡献、离散格式与一致通量 |
| 代数及并行后端 | 线性求解、线性收敛、通信及全局归约 |
| 运行时与 IO | MPI 生命周期、配置读取、显式时间服务和显式结果写入 |

运行时不打印、不推断物理收敛、不根据成功返回码自动输出。
应用可向 runApplication 提供错误打印回调。求解器使用 primaryProcess()
控制单进程打印，使用全局 diagnostics 保证各 rank 的停止决策一致。

`Case` 管理命名场；局部数学场可用 auto 创建。赋值复制数值并保留目标的边界约束。
中间派生场需要计算边界迹时，直接接收 math 返回值，或对复用目标显式设置
useCalculatedBoundary()。不向物理层暴露 LDU、CSR、Eigen 或 MPI。

完整 API 契约见 [过程式 DSL](procedural-dsl.md)，编译与算例准备见
[开发指南](solver-development.md)。旧 eqn 表达式和 loop 接口只保留兼容回归，
不作为新求解器的开发模式。
