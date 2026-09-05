# 方程驱动与算法驱动：同一框架的两种组织方式

两种方式不是两套架构，也不需要共同的 Solver 基类。普通作者都从 Case、命名 Field、
eqn/math、solve 和数学 diagnostics 出发；RunTime、矩阵和 MPI 由框架维护者负责。

## 公共的是框架语言，不是具体求解器

`include/babelsim/` 只发布 Case、Field、边界、eqn、math、solve、diagnostics 和应用注册等
通用接口，不发布 Heat、Transport、稳态 SIMPLE 或瞬态 SIMPLE 的类。具体求解器是使用这门
“框架语言”编写的应用，而不是供另一个 Solver 继承或调用的 SDK。

- Heat、Transport 足够简单，各自在一个 `main.cpp` 中完成方程与注册；
- 稳态和瞬态 SIMPLE 需要跨步骤共享数学状态，只在自己的 `src/physics/<solver>/` 中使用
  私有 `algorithm.h/.cpp`、`momentum.cpp`、`pressure.cpp` 和 `state.h`；
- 普通 Solver 作者不能包含这些私有头，也不需要理解它们；应直接组合公共数学 API；
- Case 用户只写 `solver heat`、`solver transport`、`solver simple` 或
  `solver transientSimple`，不会构造具体 Solver 对象。

这种边界避免“内置算法的当前 C++ 组织”变成必须长期兼容的公共抽象，同时保证所有新 Solver
仍可使用同一套 Field、Equation、离散、代数和 MPI 执行设施。

## 方程驱动

源码直接描述一个或几个 PDE。完整实际入口见
[heat/main.cpp](../src/physics/heat/main.cpp) 和
[transport/main.cpp](../src/physics/transport/main.cpp)。
一个函数已经包含读物性、取得场、时间循环和方程，不再另建 Case reader 或运行适配文件。

```cpp
while (problem.loop()) {
    if (!solve(eqn::ddt(C) + eqn::div(phi, C) ==
               eqn::laplacian(D, C) + eqn::source(S)).converged()) return 2;
}
```

## 算法驱动

本质是在同一时间层或稳态迭代里组织方程、修正与收敛。
[双场例子](../tests/examples/coupled_scalar.cpp) 用一个普通函数完成交替耦合，
没有新增类、框架对象、解析器或通信文件；SIMPLE 较复杂，才在自身模块内使用私有算法对象：

```cpp
while (simple.loop()) {
    simple.solveMomentum();
    simple.solvePressure();
    simple.correctVelocity();
    simple.correctFlux();
    simple.checkContinuity();
}
```

不要为了形式统一把普通 PDE 也做成继承式算法类；也不要要求新耦合算法照搬 SIMPLE 的全部文件。

## 两种方式共享的契约

1. Case 拥有物理场和命名中间场，引用在运行期间稳定。
2. Field 的赋值、缩放和积是数学操作，不是让作者遍历存储。
3. eqn 表达方程贡献，math 表达显式求值，真正执行集中在 solve/evaluate。
4. 时间历史按物理步推进，内迭代不会改变旧时间层。
5. diagnostics 返回全局量；每个耦合未知量都要被正确纳入停止判据。
6. Case 按时间步统一保存已完成结果；失败提前返回，不能写成收敛结果。
7. 新算子才需要修改离散层；组合已有算子不改 MPI 或线性代数。

## 降低学习曲线，而不是隐藏物理决策

作者仍需决定：隐式/显式项、BC、物性、耦合顺序、松弛、时间步和收敛准则。
这些是数学与数值算法，框架不应该替作者猜测。
被移除的是反复编写 reader、创建执行对象、维护历史/通信/矩阵和组织并行输出的负担。

新增独立求解器的流程只有：自己的函数与一行注册、通用 runApplication 入口 → 一个 Case →
链接已构建框架 → 串行验证 → MPI 对比。各 Solver 在自己的源文件注册，不再有内置应用选择表。
完整可执行步骤见 [开发指南](solver-development.md)；维护者的责任边界见
[架构](architecture.md)。make test-external 将真实 Solver 复制到仓库外构建，
覆盖方程驱动、双场算法驱动和矢量响应；同时确认底层存储访问不能编译。

这遵循 OpenFOAM 的场、方程、算法与 Case 分离思想，而不是复制其类型体系。
[官方 simpleFoam 主程序](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/applications/solvers/incompressible/simpleFoam/simpleFoam.C)
展示了这种以算法步骤组织代码的方式；BabelSim 用普通函数和少量对象实现自己的边界。
