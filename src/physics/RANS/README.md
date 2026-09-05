# RANS 物理模块

本目录是不可压缩 SIMPLE 的私有湍流闭合模块，不属于 BabelSim 公共 Solver API。
它只负责湍流输运变量和有效动力黏度：

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
`transient_simple`。RANS 不访问 MPI、halo、矩阵、装配或 Field 存储。

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
Picard 显式线性化；它不改变模型方程，但稳定性通常比隐式 sink 线性化更依赖
欠松弛。

SA 的 `wallDistance` 必须表示单元中心到最近真实壁面几何的最短距离；不能用
沿网格线搜索或最近单元中心距离代替。k-omega 的 `omega` 近壁值和标准高雷诺数
k-epsilon 的壁面处理必须由具体 Case 按所用网格与近壁策略给出。本模块当前不
伪装成自动壁函数系统。

## 添加模型

新增模型只需要在本目录工作：

1. 新建一个 `.cpp`，从 `Model` 派生；
2. 在构造函数中通过 `Case::scalarField/vectorField/tensorField` 声明模型变量和
   可复用工作场；
3. 在 `correct()` 中只使用公开 `Field`、`math`、`eqn` 和 `solve()` 表达输运；
4. 用 `setEddyViscosity()` 更新 `muEffective`；
5. 在 `model.h` 声明一个局部构造函数，并在 `model.cpp` 的模型选择处增加一个
   分支。

不应修改 Mesh、Field、离散、线性代数、Runtime、MPI，也不应把模型头文件放进
`include/babelsim`。这样新增闭合模型不会扩大普通 Solver 作者的公共概念面。
