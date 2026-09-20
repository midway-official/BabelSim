# BabelSim DSL 与运行时用户手册

面向二次开发：只用一个 C++ 源文件、公开头文件和预编译库，就能写出新的物理求解器并跑
MPI。本文覆盖当前 DSL 的全部算子、运行时 API、配置文件键和文件格式，语义以
`include/babelsim/*.h` 为准（本文与头文件不一致时以头文件为准）。

阅读顺序建议：第 1 节跑通最小例子 → 第 2 节建立心智模型 → 写自己的方程时查
第 6/7 节（math/equ）→ 配置与文件格式查第 3/4 节 → 出问题查第 15 节附录。

---

## 1. 最小可运行 Solver

一个求解 `ddt(T) = laplacian(k, T) + Q` 的完整 Solver（外部构建，不修改框架）：

```cpp
// solver.cpp
#include "babelsim/application.h"   // SolverRegistration / runApplication / SolverResult
#include "babelsim/case.h"          // Case / readLinearControl / write
#include "babelsim/equ.h"           // equ::Equation 与离散项
#include "babelsim/time.h"          // time::start / time::history
#include "babelsim/monitor.h"       // monitor::Reporter
using namespace babelsim;

SolverResult heat(Case& problem) {
    const monitor::Reporter reporter("heat");
    ScalarField& T = problem.scalarField("T");          // 读 fields/initial/T.field
    const double k = problem.physics().nonnegative("conductivity");
    const double Q = problem.physics().number("source");
    const auto linear = readLinearControl(problem, T);  // solution.bs 的配置（含 T 覆盖）
    const int writeInterval = readWriteInterval(problem);

    auto time = time::start(problem);                   // TimeStepper
    auto history = time::history(T);                    // 显式历史场
    auto equation = equ::createEquation(T);             // 绑定 T 的空方程
    problem.validate();                                 // 校验配置，关闭声明阶段

    while (time.value() < time.end()) {
        time.advance();
        history.save(T, time.dt());
        equation.reset();                               // 复用分配，清空项
        equ::ddt(equation, 1.0, history);               // ddt(T)
        equ::laplacian(equation, k, -1.0);              // -div(k grad T) 在左端
        equ::source(equation, Q);                       // 右端体源
        const auto result = equ::solve(equation, linear);
        reporter.record({{"time", time.value()}, {"residual", result.relative_residual}});
        if (!result.converged()) return SolverResult::notConverged();
        if (time.step() % writeInterval == 0 || time.finished()) write(problem, time);
    }
    return SolverResult::completed();
}

const SolverRegistration registration("myHeat", heat);

int main(int argc, char* argv[]) {
    return runApplication(argc, argv);
}
```

构建与运行（仓库外也能做：`include/` + `build/libbabelsim.a` 即可）：

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Iinclude -c solver.cpp -o solver.o
mpic++ solver.o build/libbabelsim.a -o my-solver
./my-solver -case cases/heat          # 或 mpirun -np 4 ./my-solver -case cases/heat
```

想成为内置命令，把同一个文件放到 `src/physics/<name>/main.cpp`：Makefile 自动收集并
链接进 `build/babelsim-solve`，不需要改启动器名单、注册宏或基类。

---

## 2. 心智模型：谁能碰什么

| 层 | 头文件 | 职责 | Solver 可否直接使用 |
|---|---|---|---|
| Case | `case.h` | 命名场、物性、配置、结果写出、时间元数据 | 是（唯一入口对象） |
| Field | `field.h` | 连续存储的 scalar/vector/tensor 场、边界条件 | 是 |
| geometry | `geometry.h` | 体积/面积/中心/法向等几何量，以普通 Field 暴露 | 是 |
| math | `math.h` | 已知场的显式数学（梯度、散度、通量、插值…） | 是 |
| equ | `equ.h` | 绑定未知量的有限体积方程（离散项 + 线性求解），唯一的方程入口 | 是（内置求解器使用） |
| time | `time.h` / `history.h` | 显式时间推进与历史层 | 是 |
| monitor | `monitor.h` | 通用观测值打印，不做收敛判断 | 是 |
| diagnostics | `solver.h` / `equ.h` | 相对残差、相对变化、通量守恒、全局归约 | 是 |
| solver | `solver.h` | `numericalMethods()`、`primaryProcess()` 与诊断入口（求解是 `equ::solve`） | 是 |
| Runtime | `runtime.h` | 时间推进、运行域与后端生命周期；内部对象 | 否（普通 Solver 不构造） |
| Backend | 内部 | halo、全局归约、稀疏装配、线性求解 | 否 |

两条硬边界：

1. **Physics 不接触实现细节。** 内置求解器源码不出现 MPI、halo、CSR/LDU、Eigen、
   `src/internal/*`、field 底层存储指针或单元索引。`make test-architecture` 与
   `make test-external` 自动检查这两条（后者在仓库外用公开头重新编译一个 Solver）。
2. **整场语义。** `math::` 的每个操作都是一个完整的整场计算：内部会自动同步输入、
   处理边界、归约 MPI；调用者永远不接触 rank、通信器或本地索引。

求解器只有两种组织方式，都共享同一套 Field/equ/math/Runtime：

- **Equation-driven**：核心就是一个或少数几个 PDE（heat、transport、Poisson）。
- **Algorithm-driven**：多个方程 + 修正步骤表达算法流程（SIMPLE、PISO、耦合）。
  这类算法把循环、松弛、收敛判据都显式写在 `main.cpp` 里，没有隐藏的 `problem.loop()`。

---

## 3. Case 与配置

### 3.1 case.bs

```text
solver   piso                  # 注册名，必须是已注册的字符串字面量
mesh     mesh/planeJet.mesh    # 相对 case 目录
fields   fields/initial        # 初值场目录（读 *.field）
physics  physics/kOmega.bs     # 物性/模型字典
methods  numerics/methods.bs   # 离散格式字典
solution numerics/solution.bs  # 线性求解与算法数值控制
control  control.bs            # 时间区间与步长
output   output.bs             # 结果目录与写出策略
ghostLayers 3                  # 可选，MPI 分区 halo 层数，>= 3
```

未知键、重复键、绝对路径、缺项都在构造时拒绝。路径必须是相对路径且不含 `.`/`..`。

### 3.2 全部配置文件键

**physics/*.bs**（物性由求解器自己声明，框架只提供读取接口）：

| 键 | 类型 | 使用者 |
|---|---|---|
| `density` | positive | heat、transport、simple/piso、RANS |
| `dynamicViscosity` | positive | simple/piso、RANS 工厂 |
| `heatCapacity` / `conductivity` / `source` | positive / nonnegative / number | heat |
| `storage` / `diffusivity` / `source` | number | transport |
| `turbulenceModel` | word：`none`(=`laminar`/`off`)、`SA`(`spalartAllmaras`)、`kOmega`(`wilcox1988`)、`kEpsilon`(`standardKepsilon`) | 动量方程求解器 |

**numerics/methods.bs**（离散格式，Case 构造时读一次）：

| 键 | 取值（默认值加粗） |
|---|---|
| `interpolation` | `linear`、**`corrected`** |
| `gradient` | `greenGauss`、**`leastSquares`** |
| `convection` | **`upwind`**、`linearUpwind`(= `secondOrderUpwind`)、`central` |
| `diffusion` | `orthogonal`、**`corrected`**、`limitedCorrected` |
| `time` | **`steady`**、`euler`、`bdf2` |

默认值来自 `Methods` 结构体初值，只对未覆盖的场生效；`time` 没有按场覆盖。

每个键都可写成两列（默认）或三列（按场名覆盖），覆盖按被离散的场名查找：

```text
gradient T greenGauss
diffusion p orthogonal
convection U central
```

调用者可用 `problem.methods()` 读取只读结果；`methods.convectionFor("T")` 之类按场查询。

**numerics/solution.bs**：

```text
scalarSolver cg     incompleteCholesky 1e-14 1e-10 1000
vectorSolver bicgstab amg             1e-12 1e-8  500 amgMaxLevels=10 amgCoarseSize=64
scalarSolver.T bicgstab ilut 1e-13 2e-9 321 ilutDropTolerance=1e-4 ilutFillFactor=3
```

线性求解行的语法固定为：

```text
<scalarSolver|vectorSolver[.field]> <cg|bicgstab> <none|incompleteCholesky|ilut|amg>
    <atol> <rtol> <maxIterations> [name=value ...]
```

可选 `name=value`：`amgMaxLevels`、`amgCoarseSize`、`amgSmoothingSteps`、
`amgRefreshInterval`、`ilutDropTolerance`、`ilutFillFactor`。`scalarSolver` 与
`vectorSolver` 必须同时存在；`scalarSolver.<field>`/`vectorSolver.<field>` 是按场覆盖。

其余键由各求解器自己读取，例如 SIMPLE/PISO 家族：`maxIterations`、
`velocityRelaxation`、`pressureRelaxation`、`nonOrthogonalCorrections`、
`continuityTolerance`、`velocityTolerance`、`momentumTolerance`、
`pressureCorrectionTolerance`；PISO 另有 `nCorrectors`；RANS 输运读取
`turbulenceRelaxation`、`turbulenceTolerance`。

**control.bs**：`startTime`、`endTime`、`deltaT` 三项必填。BDF2 要求等步长且
`endTime` 是 `deltaT` 的整数倍。

**output.bs**：

```text
directory results        # 相对 case 目录，必须相对
timeName final           # 最终时刻目录名（给 -time <name> 时用 name）
writeInterval 100        # 每多少步写一次；省略时每步都写；末步总写
writeFields U p k        # 可选：白名单
excludeFields p          # 可选：黑名单
```

未指定白名单时，默认写出“文件加载的 cell 场”；程序创建的场必须被 `output.bs` 列名
或 `problem.output(field)` 显式选入。面向量/面张量等 face 场不进入结果文件。

### 3.3 读取参数（Parameters）

`problem.physics()` 与 `problem.solution()` 返回只读 `Parameters`：

| 方法 | 语义 |
|---|---|
| `contains(key)` | 是否配置 |
| `word(key[, fallback])` | 字符串 |
| `boolean(key[, fallback])` | 布尔 |
| `number(key[, fallback])` | 有限实数 |
| `positive(key[, fallback])` | > 0 |
| `nonnegative(key[, fallback])` | >= 0 |
| `fraction(key, fallback)` | [0, 1] |
| `integer(key[, fallback[, min, max]])` | 整数（可带范围） |
| `inspect(key)` | `{configured, consumed, line}`，调试生效配置 |
| `requireAllUsed()` | 拒绝拼错的键；由 `problem.validate()` 调用 |

`problem.validate()` 只校验（不推进时间、不写结果、不关闭声明阶段）；它同时校验
Case 的所有配置字典里没有“配了但没人读”的键，改完配置立刻能发现拼写错误。

### 3.4 场文件格式 `fields/initial/<name>.field`

```text
field U
{
    type vector            # scalar | vector | tensor
    location cell          # 目前只支持 cell
    internal uniform (1 0 0)     # 初值：只能是均匀值
    boundary
    {
        inlet  { type fixedValue  value (1 0 0) }
        outlet { type zeroGradient }
        wall   { type fixedValue  value (0 0 0) }
        top    { type fixedGradient value (0 0 1) }
        side   { type symmetry }
        far    { type inletOutlet value (0 0 0) }
    }
}
```

边界类型（`field.h::BoundaryType`）：

| 类型 | 别名 | 数学语义 |
|---|---|---|
| `fixedValue` | `dirichlet` | 面值固定为 `value` |
| `fixedGradient` | `neumann` | 外法向导数固定为 `value` |
| `zeroGradient` | — | 面值取 owner 单元值（零法向梯度） |
| `symmetry` | `mirror` | 镜像：标量取 owner 值，向量去掉法向分量 |
| `inletOutlet` | — | 依赖面通量符号：流出用零梯度，流入用 `value` |

两个容易踩的点：

- 每个 patch 都必须显式配置（processor patch 例外，由 halo 提供）。
- `inletOutlet` 是**依赖通量符号**的混合条件。它需要“当前面通量”作为上下文：
  方程层的对流项会自动提供；独立用 `math::` 求值时，调用者必须先
  `field.setBoundaryFlux(phi)` 传入面通量，否则抛异常。
- 逐 patch 常值意味着**不能表达空间非均匀入口剖面**。需要 tanh/阶梯型入口时，把它
  切成多条与网格线对齐的 patch，每段给一个常值（见 `essay/plane-jet-komega-trip/`）。

### 3.5 网格文件格式 `BABELSIM_MESH 2`

```text
BABELSIM_MESH 2
vertices <n>
<x y z>            # 每行一个顶点，17 位有效数字
cells <m>
<8 个顶点索引>      # 六面体，槽位顺序见下
patches <p>
patch <name> <kind> <faceCount>
<4 个顶点索引>      # 每个边界面一行，顺序任意（读取时会按外法向自动定向）
end
```

- 单元顶点槽位：`(i,j,k), (i+1,j,k), (i+1,j+1,k), (i,j+1,k), (i,j,k+1), (i+1,j,k+1), (i+1,j+1,k+1), (i,j+1,k+1)`；
  单元体积必须为正（逆时针定向）。
- `kind` 取 `generic`、`wall`、`inlet`、`outlet`、`symmetry`（`mirror`）、`processor`。
  kind 只描述角色，边界条件仍由 `.field` 文件按 patch 名给。
- 每个外部面必须恰好属于一个 patch；非流形面、重复面会被拒绝。
- 二维算例 = z 方向 1 层单元 + front/back 用 `symmetry`（框架不维护独立二维算子）。

---

## 4. Field API

```cpp
ScalarField/VectorField/TensorField    // Field<double>/Field<Vec3>/Field<Tensor3>
```

| 成员 | 语义 |
|---|---|
| `mesh()` / `location()` / `name()` / `size()` | 只读元数据 |
| `fill(value)` | 整场填充（含边界迹），halo 直接有效 |
| `evaluate(function)` | 按**空间位置**定义场：`function(Vec3) -> T`，用于初值/物性/源 |
| `evaluate(source, function)` | 逐点映射，输入输出布局必须相同：`function(T) -> R` |
| `evaluate(a, b, function)` | 双场逐点核：`function(A, B) -> T` |
| `assign(source)` | 数值拷贝（保留本场的 mesh/name/边界约束） |
| `assignScaled(factor, source)` / `addScaled(factor, source)` | 缩放赋值 / 累加 |
| `assignProduct(coef, source)` / `addProduct(factor, coef, source)` | 逐点乘（可为原位：`f.assignProduct(a, f)`） |
| `operator= += -= *= /=` | 值语义；`Field copy = U;` 是快照，不共享存储 |
| `boundary(patch_index 或 patch 名)` | 读写边界条件对象 |
| `setBoundary(patch, condition)` | 设置边界条件 |
| `useCalculatedBoundary()` | 把当前边界迹冻结进本场：之后的逐点代数同时作用于单元值与边界迹（派生场专用） |
| `calculatedBoundary()` | 是否处于 calculated 模式 |
| `setBoundaryFlux(flux)` | 为 `inletOutlet` 提供面通量上下文 |

语义要点：

- **未知场 vs 已知场**：`createEquation` 绑定的未知场值更新会保留其边界约束；
  派生/中间场要参与逐点代数并携带边界迹，必须先 `useCalculatedBoundary()`
  （`math::` 的所有即时算子都会自动这么做）。
- **`evaluate(function)` 只依赖位置**，禁止依赖调用次数或分区，因此并行下天然一致。
- 场的长度由 `(mesh, location)` 唯一决定，任何越界/重绑定都会在 halo、算子、输出入口
  抛 `field storage invariant is violated`。

便捷构造齐次边界场（压力修正这类"同边界类型的零场"）：

```cpp
auto pPrime = field::homogeneousLike(p, "pPrime");   // 固定值->零固定值，symmetry 保留，其余零梯度
```

---

## 5. geometry：几何量的场视图

```cpp
const auto V   = geometry::cellVolumes(problem.mesh());       // ScalarField(cell)
const auto Sf  = geometry::faceAreaVectors(problem.mesh());   // VectorField(face)，外向/owner->neighbour
const auto Af  = geometry::faceAreas(problem.mesh());         // ScalarField(face)
const auto C   = geometry::cellCentres(problem.mesh());       // VectorField(cell)
const auto Cf  = geometry::faceCentres(problem.mesh());       // VectorField(face)
const auto n   = geometry::faceUnitNormals(problem.mesh());   // VectorField(face)
```

它们是只读快照，可以像普通场一样参与 `math::` 运算（例如 `V / equation.diagonal()`
得到 SIMPLE 的 `rAU`）。

---

## 6. math：已知场上的显式数学

全部函数返回**独立的即时结果**（立即求值，不是延迟表达式），返回位置见下表。
`evaluate(op, result)` 可复用调用者的目标场，`math::` 会同步输入并覆盖结果。

| 调用 | 返回位置 | 语义 |
|---|---|---|
| `grad(f)` | cell 向量 / cell 张量 | 单元梯度（标量→向量，向量→张量） |
| `normalGradient(f)` | face 标量 | 面外法向梯度；`normalGradient(f, g)` 复用已有单元梯度 |
| `flux(U)` | face 标量 | 有向面积积分通量 `Sf·U`（cell 向量先插值，face 向量直接用） |
| `flux(k, f)` / `flux(k, f, g)` | face 标量 | 数学正定扩散通量 `k·grad(f)·Sf`，**不是** `-k·grad(f)·Sf` |
| `div(faceFlux)` | cell 标量 | `Σ(outward flux)/V`；**只接受 face 标量**，cell 标量会被拒绝 |
| `div(U)` | cell 标量 | 速度散度 |
| `div(T)` | cell 向量 | 张量散度 `(div T)_i = ∂T_ij/∂x_j` |
| `div(phi, f)` / `div(phi, U)` | cell 标量/向量 | 显式对流 `div(phi·f)`（phi 是 face 通量，f 是已知 cell 数据） |
| `interpolate(f)` | face | 只做插值（含配置的插值格式） |
| `reconstruct(f, g)` | face | 用梯度做面重构 |
| `laplacian(f)` / `laplacian(k, f)` | cell 标量 | `div(k grad f)`；k 可为常数或 cell 场 |
| `dot(U, U)` | cell 标量 | 逐点内积 |
| `cross(U, U)` | cell 向量 | 逐点叉积 |
| `transpose(T)` / `trace(T)` / `isotropic(s)` | cell 张量/标量 | 张量代数 |
| `max(f, bound)` / `sqrt(f)` | cell 标量 | 逐点函数 |
| `map(f, function)` | 同位置 | 任意逐点函数（含边界迹） |
| `add(increment, target[, region])` / `subtract(...)` | 原地 | `target += increment`；`FaceRegion::Interior` 可只更新内部面 |
| `subtract(coef, grad(f), target)` | 原地向量 | `target -= coef·grad(f)`，SIMPLE 压力修正用 |
| `add(flux(U), target, region)` / `subtract(flux(k,f), target, region)` | 原地 face | Rhie–Chow 增量式组合（只动 `FaceRegion::Interior`，保留物理边界通量） |
| 运算符 | 新场 | `a+b`、`a-b`、`-a`、`scalar*field`、`field/scalar`、`scalar/scalar`、`scalar±double` |
| `sum/integral/max/normL2` | 标量 | 全局归约：`integral` 用 cell 体积或 face 面积加权，`normL2` 只对 owned 实体 |

`FaceRegion::All`（默认）包含物理边界面；`Interior` 只含内部面（含跨分区的内部面），
用于"物理边界通量保持不变"的增量组合。几何选择与 MPI 分区无关。

---

## 7. equ：绑定未知量的方程层

```cpp
auto A = equ::createEquation(U);   // 绑定未知场 U 的空方程；不推断任何项
```

| 调用 | 作用 |
|---|---|
| `reset()` | 保留分配与未知绑定，清空所有项 |
| `ddt(A, capacity, history)` | `capacity·d(x)/dt`，时间格式取自 `methods.time`；BDF2 自动用历史两层 |
| `ddt(A, capacity, previous, dt[, method, older, previous_dt])` | 无 History 对象时的显式时间层形式 |
| `div(A, phi[, scale])` | 左端 `scale·div(phi·x)`；`phi` 是**有向积分面通量**，并成为该方程唯一的边界通量上下文 |
| `laplacian(A, coefficient, multiplier)` | `multiplier·div(coefficient·grad x)`；左端扩散通常取 `-1` |
| `reaction(A, coefficient)` / `reaction(A, field, scale)` | 左端局部隐式反应项（如 `β*ω*k` 里的 `β*ω`） |
| `source(A, value)` / `source(A, field[, scale])` | 右端显式源 |
| `relax(A, previous, alpha)` | 标准欠松弛：`aP /= alpha`，`b += (aP_new-aP_old)·previous` |
| `scale(A, factor)` / `add(A, B[, factor])` / `addDiagonal(A, f)` / `addRhs(A, f)` | 方程级代数（已积分系数上操作，不再乘体积） |
| `apply(A, x)` / `residual(A, x)` | 只观察：`A·x` 与 `b-A·x` |
| `A.diagonal()` | 积分后的对角 `aP`（cell 场） |
| `A.rhs()` | 积分后的右端 `b` |
| `A.copy()` | 独立副本（系数相同、绑定同一个未知） |
| `A.reference(cell, value)` / `A.referenceIfUnanchored(value)` | 标量方程零空间规范（Neumann 压力等） |
| `equ::faceFlux(A, x)` | 标量扩散方程的左端面通量（与装配时冻结的梯度/边界一致），用于修正通量 |
| `equ::solve(A[, config])` | 求解；不传 config 时用运行时的默认线性配置 |
| `equ::solve(A, x[, config])` | 显式目标场（必须与方程绑定的未知一致，否则抛异常） |
| `diagnostics::relativeResidual(A, x)` | `‖b-Ax‖ / max(‖Ax‖+‖b‖, 1e-30)`；**必须在欠松弛之前**测，否则不是原方程残差 |

时间离散：

- `methods.time = euler`：一阶后向欧拉。
- `methods.time = bdf2`：二阶，需要两个历史层；第一步自动退化为欧拉。
- `methods.time = steady`：`ddt` 直接不加项（稳态算法应当走 SIMPLE 家族）。

`capacity`（密度×比热、或体积热容）可以是常数或 cell 场。传 `time::History` 时框架
自己取 `previous/older/dt`；历史层只由 `History::save` 推进。

写方程的推荐读法——动量预测方程：

```cpp
momentum.reset();
equ::ddt(momentum, rho, U_old);                 // ρ ∂U/∂t
equ::div(momentum, phi, rho);                   // ρ div(phi U)
equ::laplacian(momentum, muEff, -1.0);          // -div(μ_eff grad U)
equ::source(momentum, -math::grad(p));          // -grad p
auto rU = diagnostics::relativeResidual(momentum, U);
equ::relax(momentum, U_previous, alphaU);       // 欠松弛（在测残差之后）
auto aP = momentum.diagonal();
auto rAU = V / aP;                              // 对角响应，Rhie–Chow 用
equ::scale(momentum, alphaU);                   // 保持 SIMPLE 行归一化
const auto solve = equ::solve(momentum, velocitySolver);
```

## 8. 方程 API 只有一层

方程只有第 7 节这一种写法。表达式式方程（`eqn::lhs == rhs` + `solve()`）、它的解释器与
方程级控制已随该层整体删除：两种写法共享同一批离散辅助函数和同一个线性后端，数值没有差别，
但同一 PDE 会有两条装配与生命周期路径（时间历史、欠松弛参照值、参考点规范、残差口径），
差别只能靠对照表解释。现在 `equ::` 是唯一入口，装配顺序、历史推进和求解时机都在源码里
直接可读。

表达式层顺带提供的便利都有等价的显式写法：

| 表达式层的便利 | 现在的写法 |
|---|---|
| 每步自动推进的时间历史 | `time::History<T>`，每个物理步 `save` 一次（第 9 节） |
| `relaxed(alpha)` / `referenceValue(v)` | `equ::relax(A, previous, alpha)` / `A.reference(cell, v)`、`A.referenceIfUnanchored(v)` |
| `solveWithResponse` 返回 V/aP | `V / A.diagonal()`，例如 `geometry::cellVolumes(mesh) / A.diagonal()` |
| `diagnostics::residual(definition)` | `equ::residual(A, x)`、`diagnostics::relativeResidual(A, x)` |
| `Sp(coefficient, field)` | `equ::reaction(A, coefficient[, scale])` |

需要明确的一点：`equ::` **不做项位置校验**，符号或位置写错不一定会报错——`laplacian` 的
`multiplier` 取 `-1` 才是常见的左端扩散（第 7 节），写反会静默改变方程。唯一的一致性检查
是边界通量上下文：`div` 绑定的面通量成为该方程的边界通量，同一方程再出现另一个通量会抛
`one equation requires one boundary flux context`；`reset()` 会同时解除该绑定。

## 9. 时间与历史

```cpp
auto time = time::start(problem);        // TimeStepper，查询无副作用
time.advance();                          // 唯一推进时间的入口
time.value() / time.end() / time.dt() / time.step() / time.finished()

auto old = time::history(U);             // 显式历史，构造时以当前场为初值
old.save(U, time.dt());                  // 每物理步恰好调用一次
old.previous(); old.older(); old.dt(); old.previousDt(); old.levels();
```

- 收敛判断、循环结构、残差打印全部由算法自己写；框架不提供 `problem.loop()`。
- `write(problem, time)` 或 `write(problem, t, step)` 在算法选定的时刻写出。
- 稳态算法用自己的 `for` 循环，不需要 `TimeStepper`（但 `problem.validate()` 仍要调）。

---

## 10. 线性求解与数值控制

```cpp
const auto config = readLinearControl(problem, U);   // 读 solution.bs（含 U 覆盖）
const auto result = equ::solve(equation, config);
```

- `LinearSolverConfig`：`solver`（`ConjugateGradient`/`BiCGSTAB`）、`preconditioner`
  （`None`/`IncompleteCholesky`/`ILUT`/`AlgebraicMultigrid`）、`absolute_tolerance`、
  `relative_tolerance`、`max_iterations`、`warm_start`、`ilut_*`、`amg_*`。
- 收敛合同（所有后端一致）：原始未预条件残差的全局 L2 范数满足
  `‖b-Ax‖ <= max(atol, rtol·max(‖r0‖, ‖b‖))`；不设隐藏的舍入下限，也不把停滞当收敛。
- `SolveResult`：`status`（`Converged`/`MaxIterations`/`NumericalFailure`）、
  `iterations`、`initial_residual`、`final_residual`、`relative_residual`、`performance`。
  矢量方程公开结果要求各分量全部收敛。
- `readLinearControl(problem, field)` 是**唯一**读取按场覆盖的方式；不传配置直接
  `equ::solve(A)` 会退回运行时默认配置，拿不到 `scalarSolver.T` 这类覆盖。

---

## 11. 诊断与监视

```cpp
// 方程级（solver.h / equ.h）
double r = diagnostics::relativeResidual(equation, x);      // 已装配方程上的 ‖b-Ax‖，测在欠松弛之前
Field<double> residual = equ::residual(equation, x);        // 需要残差场本身时直接取
double dU = diagnostics::relativeChange(U, U_previous);     // ‖ΔU‖/‖U‖
double dP = diagnostics::relativeMagnitude(pPrime, p);      // 修正量相对幅值
FluxBalance mass = diagnostics::fluxBalance(phi);           // 每个控制体的通量守恒误差
bool ok = diagnostics::all(std::isfinite(r) && std::isfinite(dU));  // 全局归约

// 通用监视（不做收敛判断，不保留数据）
const monitor::Reporter reporter("PISO", /*iterationInterval=*/100);
reporter.record({{"time", t}, {"residual", r}});
reporter.iteration(iter + 1, maxIterations, urgent, {
    {"mass", mass.relative}, {"dU", dU}, {"rU", r}, {"converged", converged}});
```

`Metric` 的值可以是 `double`/`int`/`bool`/字符串。`Reporter` 打印首轮、周期轮与末轮；
`urgent` 只影响是否立即打印，**不驱动控制流**——收敛判据必须由求解器自己写在代码里。

---

## 12. 输出、结果与后处理

```cpp
problem.output(mut);                 // 把程序创建的场加入写出白名单
problem.output(U, false);            // 关闭某个输入场的写出
write(problem, time);                // 按 output.bs 写出当前时刻
problem.performance();               // 关键路径（最大 rank）性能快照
problem.localPerformance();          // 本 rank 快照
```

结果目录布局（`results/<run>/`）：

```text
results/<run>/<time>/rank-0000/U.csv        # global_id,x,y,z,value0[,value1,value2]
results/<run>/<time>/rank-0000/metadata.bs  # format/time/rank/ranks/global_cell_count/fields
results/<run>/rank-0000/...                 # 最终时刻（timeName 或 -time 名）
```

- 每个 rank 只写 owned cell；ghost 不写。CSV 第二行起是 cell，`global_id` 全局唯一。
- 用不同进程数重跑同一结果目录会被拒绝；用 `-time <name>` 隔离实验。
- 后处理：`build/babelsim-post -case <dir> [-time <name|latest|all>] -format vtk tecplot`
  产出 `.vtu` / `FEBRICK` / `series.pvd`。
- 性能：`-performance <dir>` 每个 rank 写一份 JSON（`solverSeconds`、`haloSeconds`、
  `krylovIterations`、`partition`、`linearSolveSeconds` 等），聚合用法见
  [性能工具说明](performance/README.md)。

---

## 13. 并行与运行时边界

- 分布式网格读取：`readDistributedMesh(path, parallel, ghostLayers)` 在 rank 0 解析
  网格，再按单元邻接图构造每个 rank 的 owned+ghost 局部 Mesh 与 processor patch。
  `ghostLayers >= 3`（见 case.bs）。库代码从已有全局网格分区时可用 `decompose()`。
- Solver 侧**不需要**做任何 MPI 调用：整场算子内部完成 halo 交换与全局归约；
  `diagnostics::all(...)` 把本地判断归约成全局结果。
- 只有 `primaryProcess()` 用于"是否由主进程打印"这类判断；算法不应据此改变数值路径。
- `numericalMethods()` 返回当前运行域的只读 `Methods`（算法据此选修正次数等），
  不接触后端。
- 结果一致性：同样的 Case 用 1/2/4 rank 跑，场值差异应在求解器容差量级
  （`make test-mpi*` 与 `tools/compare_parallel_results.py` 覆盖这条）。

---

## 14. 内置求解器与 RANS

完整用法（方程、算法步骤、全部配置键与默认值、场与边界、收敛与失败语义、验证证据）
见 [内置求解器手册](solvers.md)；本节只做速查。

注册名 → 入口 `src/physics/<name>/main.cpp`：

| solver | 类型 | 方程/算法 | 读取的配置键 |
|---|---|---|---|
| `heat` | Equation | `ddt(ρc_p T) = div(k grad T) + Q` | `density, heatCapacity, conductivity, source` |
| `transport` | Equation | `ddt(storage·C) + div(phi C) = div(D grad C) + S` | `storage, diffusivity, source` |
| `simple` | Algorithm | 稳态 SIMPLE（动量预测 + 压力修正 + Rhie–Chow，私有实现） | `density, dynamicViscosity, maxIterations, velocityRelaxation, pressureRelaxation, nonOrthogonalCorrections, continuityTolerance, velocityTolerance, momentumTolerance, pressureCorrectionTolerance` |
| `transientSimple` | Algorithm | 瞬态 SIMPLE（步内迭代到收敛） | 同上 |
| `piso` | Algorithm | 瞬态 PISO：一次动量预测 + `nCorrectors` 次压力修正 | 上面的 + `nCorrectors`（`maxIterations 1` 即标准 PISO；`deltaT` 必须 > 0） |
| `hello`/其它 | — | `src/physics/<name>/main.cpp` 里 `SolverRegistration` 即可注册 | — |

RANS 通过 `src/physics/RANS/api.h` 与动量求解器交互（SIMPLE/PISO 内部）：

```cpp
auto turbulence = rans::load(problem, U, phi);      // 依 physics 的 turbulenceModel 构造
const ScalarField& muEff = turbulence.effectiveViscosity();     // mu + mut
equ::laplacian(momentum, muEff, -1);                            // 隐式扩散
equ::source(momentum, math::div(turbulence.deviatoricStressRemainder(U)));  // 显式余项
if (turbulence) { turbulence.saveOld(dt); auto r = turbulence.solveTransport(); }
```

模型名（`turbulenceModel`，大小写与下划线被忽略）：`none`/`laminar`/`off`、
`SA`/`spalartAllmaras`、`kOmega`/`wilcox1988`、`kEpsilon`/`standardKepsilon`。

| 模型 | 需要加载的场 | 模型常数（physics/*.bs） |
|---|---|---|
| SA | `nuTilda`（+ 壁面距离场） | `saCb1, saKappa, saSigma, saNuTildaMin, saWallDistanceMin` |
| k-ω | `k`, `omega` | `kOmegaBetaStar, kOmegaBeta, kOmegaGamma, kOmegaSigmaK, kOmegaSigmaOmega, kMin, omegaMin` |
| k-ε | `k`, `epsilon` | `kEpsilonCmu, kEpsilonC1, kEpsilonC2, kEpsilonSigmaK, kEpsilonSigmaEpsilon, kMin, epsilonMin` |

模型自己创建并登记 `mut`（涡黏）与 `muEffective`（有效黏度）两个输出场；输运方程用
`turbulenceRelaxation` 欠松弛、`turbulenceTolerance` 判稳。k-ω/k-ε 的壁面 ω、ε 需要
按解析渐近式在 `.field` 里给（见 `essay/cube3d-komega/`、`cases/naca0012/`）。

---

## 15. 开发检查清单

1. 只包含公开头（`babelsim/*.h`），不 include `src/internal/*`。
2. 物性/算法参数全部经 `problem.physics()/solution()` 的类型化读取，拼错键会被
   `requireAllUsed()` 拒绝。
3. 场只用 `scalarField/vectorField/tensorField`（加载）、`create*Field`（程序创建）、
   `existing*Field`（复用）三类入口，不在声明阶段结束后再新建场。
4. 方程项顺序按 PDE 顺序书写，`laplacian` 左端用 `-1`，`source` 只放显式项。
5. 残差在欠松弛之前测；每步都检查 `SolveResult::healthy()`，`NumericalFailure` 立即返回。
6. 循环、收敛判据、写出时刻自己写；不用 `urgent` 或打印结果驱动控制流。
7. 时间历史每步恰好 `save` 一次；BDF2 只用在等步长算例。
8. 面向量/面张量等中间场用 `createFace*Field`，不入结果文件。
9. 新 Solver 最低验证：`make test-architecture`、`make test-external`、
   `make test-workflow`，再补一个 case 与 1/2 rank 结果对照。
10. 数值结论必须附算例与命令（`docs/reports/` 的既有报告是格式范例）。

验证入口：

```bash
make -j4                  # 只构建 lib + solve + post
make test                 # 串行单元/集成测试
make test-architecture    # 分层与头依赖门禁
make test-external        # 仓库外用公开头编译 Solver（含负向 API 检查）
make test-workflow        # 新 Solver 工作流、时间序列、ParaView 读取
make test-mpi             # MPI 网格/halo/算子/线性求解
make test-rans            # RANS 方程与常数
make validate-cavity      # Ghia 腔体快速回归
make validate-poiseuille  # Poiseuille 解析解比较
```

---

## 附录 A：旧接口 → 当前接口

| 旧写法 | 当前写法 |
|---|---|
| `math::copy/history/saveOld` | Field 值语义（`const Field prev = U;`）与 `time::History::save` |
| `math::correction(field)` / `createHomogeneousField` | `field::homogeneousLike(field[, name])` |
| `equ::Matrix` | `equ::Equation<T>` |
| `Equation::volumeScaledInverseDiagonal()` | `V / equation.diagonal()` |
| `faceField(...)` 混合入口 | `createFaceScalarField/createFaceVectorField/createFaceTensorField` |
| `parameters.entry().tokens` | 类型化 getter（`positive/number/word/...`） |
| `case.loop()` / `simple.solve()` | 求解器自己的显式 `while/for` 循环 |
| `LinearSystem` / `assemble()` | 内部装配（`SparseAssembly`）；公开侧不再暴露矩阵 |
| `fvm::` / `fvc::` 命名空间 | `equ::`（隐式方程项）与 `math::`（显式数学） |
| `eqn::` 表达式式方程 / `solve(definition)` / `EquationControl` / `relaxed` / `referenceValue` | `equ::` 显式组装，配合 `time::History`、`equ::relax`、`reference`/`referenceIfUnanchored`、`V/A.diagonal()`（第 8 节） |
| `Case::loadMethods()` / `Case::createFaceField()` | `problem.methods()` / `problem.createFace*Field()` |
| `properties()` | `physics()` |
| `gmres` 线性求解器 | 只提供 `cg` / `bicgstab`（+ IC/ILUT/AMG 预条件） |

## 附录 B：常见语义陷阱

| 现象 | 原因 |
|---|---|
| `inletOutlet requires a boundary flux context` | 独立用 `math::` 求值时没先 `setBoundaryFlux(phi)` |
| `one equation requires one boundary flux context` | 同一个方程里 `div` 用了两个不同面通量；一个方程只能有一个边界通量上下文（`reset()` 后重新 `div` 即可换通量） |
| `solve target differs from matrix unknown` | `equ::solve(A, x)` 的 x 不是创建方程时绑定的场 |
| `field storage invariant is violated` | 手工构造/移动场导致长度与 `(mesh, location)` 不符 |
| 压力场解出常数漂移 | 全 Neumann 标量方程没设 `referenceIfUnanchored` / `reference` |
| `relativeResidual` 与观察到的收敛不符 | 残差在 `relax` 之后测了（应在其之前） |
| 瞬态算例第二步发散 | 设了 `amgRefreshInterval>1`，复用了过期 AMG 层级 |
| 结果目录被拒绝 | 同一 `<run>` 换进程数重跑；用 `-time <name>` 隔离 |
| 细网格上第一步 `k` 爆掉 | 入口阶跃剪切层首层应变过大而 ω 不足（`S²/ω > βω`），见 `essay/plane-jet-komega*/` 的启动分析 |