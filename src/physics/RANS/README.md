# RANS 物理模块

本目录是不可压缩 SIMPLE 的私有湍流闭合模块，不属于 BabelSim 公共 Solver API。
它负责湍流输运变量、派生系数和松弛前输运方程残差；SIMPLE 用有效动力黏度构造完整偏应力：

```text
修正后的 U、phi
        ↓
RANS 输运方程
        ↓
mu_t
        ↓
muEffective = mu + mu_t
        ↓
下一轮 SIMPLE 动量方程
```

速度、压力、压力修正和连续性检查仍完全属于 `simple` 或
`transient_simple`。RANS 使用公开 equ API 创建和组装方程，不访问 MPI、halo 或底层矩阵/Field 存储。

## 文件与执行顺序

三个模型各自在一个源文件中完成全部闭合和输运：

| 文件 | 模型与物理变量 |
| --- | --- |
| `spalart_allmaras.cpp` | 标准正变量 SA，含 ft2、无 primary trip term；nuTilda |
| `k_omega.cpp` | Wilcox1988m；k、omega |
| `k_epsilon.cpp` | 标准高雷诺数 k-epsilon；k、epsilon |

每个 `solveTransport()` 可以顺序阅读：保存本轮值 → 计算生成/耗散/扩散系数 →
逐项组装 → 记录松弛前残差 → 松弛与求解 → 限制未知量下限 → 更新 mut 与 muEff。
两方程模型的两条方程采用同一轮冻结的闭合系数，保持原有 Picard 更新策略。
各方程每次调用只组装一次；闭合系数是局部数学场，不注册为 Case 的 ransWork 临时变量。
仅物理场、配置、线性求解设置和历史跨调用保存。

模型彼此不调用，不共享输运实现或数值基类。`Model` 是纯接口，`model.cpp`
仅选择模型并连接输入输出。报告中的汇总函数只处理诊断数据，不执行模型算法。

## 与 SIMPLE 的交互

```cpp
auto turbulence = rans::load(problem, U, phi);
const auto& muEff = turbulence.effectiveViscosity();

// 瞬态：每个物理时间步开始时调用一次，稳态不调用。
turbulence.saveOld(time.dt());

// 在 SIMPLE 显式选定的位置求解模型输运，随后更新 muEff。
const auto report = turbulence.solveTransport();
const double residual = report.initialResidual();
const double change = report.relativeChange();
```

`report.equations` 分别保存每个未知量的名称、线性求解结果、归一化初始残差和
限幅后的相对变化。`initialResidual()`/`relativeChange()` 返回无量纲指标的最大值；
不会将 k 和 omega/epsilon 的有量纲绝对残差相加。`linearConverged()` 仅表示
各线性系统收敛，不代表模型与 SIMPLE 的物理耦合已收敛。后者及所有打印由求解器判断。
模型接口不自动推进时间，不隐式增加迭代。有效黏度的单位为动力黏度：muEff = mu + mut。

## Case 配置

模型名称和模型常数全部写在 `case.bs` 的 `physics` 文件所指向的同一个
`.bs` 字典中。未写 `turbulenceModel` 时，为兼容原有层流算例，等价于
`none`。

```text
density 1.0
dynamicViscosity 1.8e-5
turbulenceModel kOmega
turbulenceRelaxation 0.5
turbulenceTolerance 1e-6
```

可选模型：

- `none`：层流，不创建湍流变量；
- `SA`：标准正变量 Spalart-Allmaras；需要初始场 `nuTilda` 和
  `wallDistance`；
- `kOmega`：Wilcox1988m 高雷诺数 k-omega；需要 `k` 和 `omega`；
- `kEpsilon`：标准高雷诺数 k-epsilon；需要 `k` 和 `epsilon`。

湍流变量的初值和边界条件仍使用普通 Field 文件，放在 Case 的
`fields/initial` 目录。模型参数则只放在 physics 字典：

```text
# 公共控制
turbulenceRelaxation 0.5
turbulenceTolerance 1e-6

# SA 可选常数
saCb1 0.1355
saCb2 0.622
saSigma 0.6666666666666667
saKappa 0.41
saCw2 0.3
saCw3 2.0
saCv1 7.1
saCt3 1.2
saCt4 0.5
saNuTildaMin 1e-14
saWallDistanceMin 1e-12

# k-omega 可选常数
kOmegaBetaStar 0.09
kOmegaBeta 0.075
kOmegaGamma 0.5555555555555556
kOmegaSigmaK 0.5
kOmegaSigmaOmega 0.5
kMin 1e-12
omegaMin 1e-12

# k-epsilon 可选常数
kEpsilonCmu 0.09
kEpsilonC1 1.44
kEpsilonC2 1.92
kEpsilonSigmaK 1.0
kEpsilonSigmaEpsilon 1.3
kMin 1e-12
epsilonMin 1e-12
```

没有写出的常数使用上述默认值。字典仍执行“所有条目必须被消费”检查，
因此拼错模型名或给当前模型填写其他模型的系数会被拒绝。

## 数值与边界条件

模型输运方程含对流项，通常应在 `solution.bs` 为标量方程选择
`bicgstab`，而不是仅适合对称正定系统的 `cg`。破坏项使用上一轮场值进行
Picard 隐式耗散线性化，通过 `equ::reaction(equation, destructionRate)` 表达。
生成项通过 `equ::source` 放在右端。SA 对净反应系数按符号拆分，在固定点保持
同一正变量 PDE。`initialResidual` 在松弛与求解之前采样，变化量在求解和下限限制后采样。
下一次正常组装检查上一次更新后的状态；不能仅凭线性成功或下限裁剪后的零变化宣称收敛。
SA 保留 S-tilde >= 1e-30 和 r in [0,10] 的数值保护，没有改成 SA-neg。

SA 的 `wallDistance` 必须表示单元中心到最近真实壁面几何的最短距离；不能用
沿网格线搜索或最近单元中心距离代替。k-omega 的 `omega` 近壁值和标准高雷诺数
k-epsilon 的壁面处理必须由具体 Case 按所用网格与近壁策略给出。本模块当前不
伪装成自动壁函数系统。

## 添加模型

新增模型只需要在本目录工作：

1. 新建一个 `.cpp`，从 `Model` 派生；
2. 构造时读取模型自己的参数和物理场，只保存必要的长期状态；
3. 在 `solveTransport()` 中用局部 math 场和 equ 操作直接表达全部方程；
4. 更新模型自己的 mut，并令 muEff = mu + mut；返回逐方程的 `TransportResult`；
5. 实现显式 `saveOld(dt)`，只在物理时间步切换时保存模型历史；
6. 在 `model.h` 声明构造函数，在 `model.cpp` 增加选择分支。

不应修改 Mesh、Field、离散、线性代数、Runtime、MPI，也不应把模型头文件放进
`include/babelsim`。这样新增闭合模型不会扩大普通 Solver 作者的公共概念面。

派生系数与梯度工作场显式使用 `useCalculatedBoundary()`，按相同数学函数传播面迹；
未知量的下限裁剪不会把 SA 的固定零壁面改成正值。各向同性 k 应力的压力约定见
[动量说明](/home/midway/BabelSim/docs/simple-solver.md)，公开论文、公式及实际验收见
[模型方程报告](/home/midway/BabelSim/docs/reports/rans-equation-verification.md)。

复现：`make test-rans`。测试包含模型系数、非均匀黏度应力、Euler/BDF2 衰减阶、
完整稳态/瞬态 SIMPLE 的 1/2/4 进程及裁剪拒绝；不等同于自动壁函数或工程壁流验证。

## 本轮重构核对的资料

- [NASA/TMBWG 标准 SA](https://tmbwg.github.io/turbmodels/spalart.html)：
  核对 fv1、fv2、ft2、fw、近壁距离、源项符号和标准模型版本。
- [NASA/TMBWG Wilcox 模型](https://tmbwg.github.io/turbmodels/wilcox.html)：
  核对 1988 系数与无交叉扩散的版本边界。
- [OpenFOAM 10 k-epsilon 实现](https://github.com/OpenFOAM/OpenFOAM-10/blob/master/src/MomentumTransportModels/momentumTransportModels/RAS/kEpsilon/kEpsilon.C)：
  核对标准闭合、扩散系数和隐式耗散；BabelSim 保持自身原有两方程冻结系数顺序。

验证采用独立标量公式构造一次输运更新的参考解，检查实际输出与逐方程残差；
不再读取模型私有临时场。Euler/BDF2 衰减、MPI 和裁剪拒绝检查仍然保留。
这些验证不证明工程壁流精度，也不提供自动壁函数。
