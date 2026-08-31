# 瞬态热传导：最小方程驱动 Solver

BabelSim 的求解器名称是 `heat`，不是 heatFoam。完整用户入口只有
`src/physics/heat/main.cpp` 一个函数：

```cpp
int runHeat(Case& problem) {
    ScalarField& T = problem.scalarField("T");
    const double rho = problem.properties().positive("density");
    const double cp = problem.properties().positive("heatCapacity");
    const double k = problem.properties().nonnegative("conductivity");
    const double Q = problem.properties().number("source");

    while (problem.loop()) {
        if (!solve(fvm::ddt(rho * cp, T) ==
                   fvm::laplacian(k, T) + fvm::source(Q)).converged()) return 2;
    }
    return 0;
}
```

直接对应
\[
\rho c_p\partial_tT=\nabla\cdot(k\nabla T)+Q.
\]

作者无需写 Case reader、构造 RunTime、配置线性对象、维护历史或编写输出。
Case 在下一次 loop 前写出已完成时间步，正常退出保证最终时刻保存。
不收敛立即返回，不把失败步伪装成有效最终结果。

## 工作流

```bash
make -j4
mpirun -np 4 build/babelsim-solve -case cases/heat -time mpi4
build/babelsim-post -case cases/heat -time mpi4/all -format vtk tecplot
```

数据在 `results/mpi4/<物理时间>/rank-*/`，
ParaView 打开 `post/mpi4/series.pvd`。
默认不指定 -time 时，序列直接在 results/<物理时间>。

output.bs 的 writeInterval 是步数间隔，默认 1。例如 endTime=0.05、deltaT=0.01、
writeInterval=2，输出 0.02、0.04、0.05，而不是仅写 final 或把 1/2/3 当作时间。

## 变系数与嵌入式接口

将 k 或 rhoCp 替换为从 Case 读入的 scalar Field，fvm 写法不变。
材料模型可在每一步前更新物性场；它仍属于物理数学代码，不应操作通信或矩阵。

thermal.h 的 ThermalProperties、ThermalFieldProperties、solveHeatStep 和
solveTransientHeat 保留给内存型数值测试/嵌入式调用。它们使用同一后端，但不承担 Case 自动输出。
新 Solver 作者优先使用上面的 Case 入口，不需要先学习这组显式 RunTime 接口。

## 与 OpenFOAM 对照

OpenFOAM-8 的 laplacianFoam 在主循环中构造 ddt/laplacian 方程，并包含非正交修正及写出步骤；
它短小是因为创建网格/场、边界离散、方法选择和执行工作已由框架提供。
BabelSim 学习这个职责分离，而不复制其头文件片段包含方式、fvMatrix、IOobject 或注册机制。
参见 [OpenFOAM-8 laplacianFoam 官方源码](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/applications/solvers/basic/laplacianFoam/laplacianFoam.C)。

| 对照点 | BabelSim |
| --- | --- |
| Field 从 Case 创建 | problem.scalarField("T") |
| 方程与数学同形 | solve(ddt == laplacian + source) |
| 数值格式来自 Case | methods.bs / solution.bs |
| 时间、历史和输出下沉 | Case::loop 调用内部 RunTime |
| 普通作者新增代码 | 一个普通函数 + 一项显式选择 |

不同点：当前 heat 的每步只求一次方程；强非线性物性或需要多次显式非正交修正时，应像
双场例子一样在该时间步内组织收敛迭代。框架已保证重复 solve 不推进历史，
但不会替物理作者猜测非线性收敛准则。
