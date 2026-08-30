# 瞬态热传导 Solver

`heatFoam` 是 BabelSim 的最小 PDE Solver 模板。它验证 Field、边界、时间、方程表达、
FVM 装配、线性求解和 MPI 都可以在不向物理代码暴露底层细节的情况下组合。

## 控制方程

\[
\rho c_p\frac{\partial T}{\partial t}
= \nabla\cdot(k\nabla T) + Q.
\]

当前 `ThermalProperties` 支持常数 \(\rho,c_p,k\) 与常数体源 \(Q\)。通用方程 API 同时
允许把 \(\rho c_p\)、\(k\) 和 \(Q\) 作为 Field 传入，因此后续温度相关物性不需要改变
Runtime、代数或 MPI 层。

## 核心源码

实际实现位于 `src/physics/heat/transient_heat_solver.cpp`。核心循环是：

```cpp
while (run_time.loop()) {
    result.linear = solve(
        run_time,
        fvm::ddt(material.volumetricHeatCapacity(), temperature) ==
            fvm::laplacian(material.conductivity, temperature) +
            fvm::source(volumetric_source));
}
```

它直接对应控制方程：左端是蓄热项，右端是扩散与体源。此文件没有 MPI、rank、ghost、
buffer、稀疏矩阵、LDU、Eigen、`mutableData()` 或手写 cell/face 循环。

如果系数是 Field，写法保持不变：

```cpp
solve(
    run_time,
    fvm::ddt(rho_cp, T)
        == fvm::laplacian(conductivity, T) + fvm::source(source));
```

其中 `rho_cp` 和 `conductivity` 是 cell scalar Field。Runtime 会同步它们并为扩散自动
产生所需面系数。

## 案例

内置例子位于 `cases/heat`：

```bash
make -j
build/babelsim-solve -case cases/heat -time serial
mpirun -np 2 build/babelsim-solve -case cases/heat -time mpi2
python3 tools/compare_parallel_results.py \
  cases/heat/results/serial cases/heat/results/mpi2 \
  --atol 1e-9 --rtol 1e-9
```

案例字典含义：

```text
case.bs                 选择 heatFoam 与各资源路径
mesh/heat.mesh          3D 结构化网格；例子是 nz=1 的二维退化网格
fields/initial/T.field  温度初值与 hot/cold/symmetry 边界
physics/thermal.bs      density、heatCapacity、conductivity、source
numerics/heat.bs        方法、时间区间、步长和 scalarSolver
output.bs               结果目录与默认时刻名
```

## 与 OpenFOAM laplacianFoam 的对照

OpenFOAM 的 `laplacianFoam`/标量输运 Solver 之所以短，是因为网格、Field、边界、离散格式、
线性求解和并行控制都由框架管理；Solver 只保留 `fvm::ddt`、`fvm::laplacian`、`solve` 和
时间循环。BabelSim 的热 Solver 采用同样的职责划分：

| OpenFOAM 思想 | BabelSim 对应 |
| --- | --- |
| `runTime.loop()` | `RunTime::loop()` |
| `fvm::ddt(T)` | `fvm::ddt(..., T)` |
| `fvm::laplacian(k,T)` | `fvm::laplacian(k, T)` |
| `TEqn.solve()` | `solve(run_time, equation)` |
| `fvSchemes/fvSolution` | `numerics/heat.bs` |
| decomposition/MPI/IO | RunTime、启动器、并行写出与后处理器 |

不同点是 BabelSim 不使用 `GeometricField`、`fvMatrix`、`tmp<>`、`IOobject` 或运行时对象
注册表。表达式只描述方程，实际 LDU 与分布式线性系统始终留在 Runtime 内部。

可参考 OpenFOAM 的 [scalarTransportFoam 源码](https://github.com/OpenFOAM/OpenFOAM-8/blob/master/applications/solvers/basic/scalarTransportFoam/scalarTransportFoam.C)
与 [laplacianFoam 文档](https://doc.openfoam.com/2306/tools/processing/solvers/rtm/basic/laplacianFoam/)。
