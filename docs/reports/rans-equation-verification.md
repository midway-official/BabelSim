# RANS 方程、公开资料与验收对照

对象：2026-09-09 审查 F1–F7 修复后的工作区。模型实现位于 `src/physics/RANS`，
动量耦合位于两个 SIMPLE 模块和 `simple_common.h`。本文区分模型方程核对、离散验证和
工程物理验证；实际命令与日志见 [修复验收记录](f1-f7-implementation.md)。

## 1. 对照范围及来源

| 模型 | 当前采用的版本 | 核对资料 |
|---|---|---|
| SA | 正变量、无 trip 源、保留 ft2 的标准 SA | Allmaras、Johnson、Spalart，ICCFD7-1902（2012）§2.1，式 (1)–(8)，以及 TMR SA 定义 |
| kOmega | Wilcox1988m，高 Re、无 SST 交叉扩散 | Wilcox，AIAA Journal 26(11), 1299–1310（1988），DOI 10.2514/3.10041；TMR 的 Wilcox1988/1988m 节 |
| kEpsilon | 标准高 Re k-ε | Launder、Spalding（1974）的书目信息及标准模型定义；OpenCFD 模型文档、OpenFOAM Foundation v7 的公开实现 |

SA 的逐项对照实际读取了 [2012 公开论文全文](https://www.iccfd.org/iccfd7/assets/pdf/papers/ICCFD7-1902_paper.pdf)
及 [TMR SA 页面](https://tmbwg.github.io/turbmodels/spalart.html)。原始 1994 文献为
Spalart & Allmaras, “A One-Equation Turbulence Model for Aerodynamic Flows”,
Recherche Aérospatiale 1, 5–21；TMR 提供[经作者许可的原论文扫描件](https://tmbwg.github.io/turbmodels/Papers/RechAerosp_1994_SpalartAllmaras.pdf)。
2012 文献明确重述基准模型，不能把其后文 SA-neg 扩展当成本实现已经具备的功能。

Wilcox 的方程和系数以 [TMR 版本定义](https://tmbwg.github.io/turbmodels/wilcox.html) 核对，
与 [1988 原论文](https://doi.org/10.2514/3.10041) 的模型对应。原论文另有公开正文转录，
OCR 公式不清处以 TMR 明确公式为准；不混用 1998、2006、SST 常数。

k-ε 的 [1974 出版社页面](https://www.sciencedirect.com/science/article/pii/0045782574900292)
可确认作者、题名和出处，但本次未取得可逐式核读的该论文全文，不能声称读过其全部公式。
实际方程/常数交叉核对依据 [OpenCFD 标准模型文档](https://doc.openfoam.com/2212/tools/processing/models/turbulence/ras/linear-evm/rtm/kEpsilon/)
和 [OpenFOAM Foundation v7 原始实现](https://cpp.openfoam.org/v7/kEpsilon_8C_source.html)，
后者公开列出生产、耗散和隐式源项。这是标准实现交叉验证，论文全文逐页对照仍是资料覆盖边界。

## 2. 公共物理约定和动量式

ρ 为常数密度，μ 为动力黏度，ν=μ/ρ；phi 为面积积分的体积通量。
令 D 为速度梯度的对称无迹部分，`strainMeasure` 计算 2D:D，模型生产 P=2μt D:D。
当 div(U)=0 时，P 就是标准不可压缩涡黏性生产项；没有额外压缩、浮力、曲率或旋转修正。

修复后，RANS 的动量扩散为：

\[
\nabla\cdot\{(\mu+\mu_t)[\nabla U+(\nabla U)^T-\tfrac23(\nabla\cdot U)I]\}.
\]

`momentumEquation()` 用隐式分量 Laplacian 加显式转置/无迹修正，两个 SIMPLE 共享私有
物理表达。层流分支仍为原先的恒 μ 不可压缩 NS 动量式。全应力的其余各向同性 k 项在
两方程模型中按 p*=平均压力+2ρk/3 的压力定义吸收；p 文件及压力边界要遵守此约定。
本实现没有自动把 p* 输出成平均压力。SA 没有 k 输运量。

与旧 F1 对应的独立制造场：U=(2y,0,0)、μe=1+x。在无散条件下，旧分量 Laplacian 为零，
新增转置贡献为 (0,2,0)。`rans_equations_test.cpp` 在 1/2/4 进程逐 owned 单元检查该值。
该项不是通过重新定义通用离散层的“黏度含义”实现，而是物理层组合张量数学。

## 3. 各输运模型与代码对应

### SA

记 q=nuTilda，涡黏度 μt=ρq fv1。常密度下的代码对应

\[
\partial_t(\rho q)+\nabla\cdot(\rho Uq)
=\nabla\cdot\big[\rho(\nu+q)\nabla q/\sigma\big]
+\rho\{c_{b1}(1-f_{t2})\widetilde S q
-(c_{w1}f_w-c_{b1}f_{t2}/\kappa^2)q^2/d^2
+c_{b2}|\nabla q|^2/\sigma\}.
\]

`updateFunctions()` 对应 2012 基准式 (1)、(4)–(6) 的 χ、fv1/fv2、ft2、S̃、r、fw；
`updateEquationFields()` 对应式 (2)–(3)；`updateViscosity()` 对应式 (1)。
cb1=.1355、cb2=.622、σ=2/3、κ=.41、cw2=.3、cw3=2、cv1=7.1、ct3=1.2、ct4=.5，
cw1=cb1/κ²+(1+cb2)/σ，r 上限为 10。没有 trip 源，不是 SA-noft2。
这些默认值与 [公开基准定义](https://tmbwg.github.io/turbmodels/spalart.html) 对应。

代码把反应项写成 a(q)q，左侧隐式放 ρ max(-a,0)q，右侧显式放 ρ max(a,0)q，
梯度平方项仍显式。两部分在非线性固定点恢复原式。q 下限、d 下限和 S̃≥1e-30 是数值
保护，不属于新增物理项；下限若持续活跃而原方程残差未满足，当前会拒绝收敛。
S̃ 的正值截断是已有稳健性处理，没有实现 2012 后文的平滑替代公式或 SA-neg。

标准无滑移壁面 q=0；标量对称面法向导数为零。q 的内部裁剪保留固定零边界，
派生迹保证壁面 μt=0、μe=μ，扩散系数为 μ/σ。wallDistance 仍须由 case 正确给出最近
真实壁面距离；本轮没有添加壁距求解器。常密度让 2012 式 (9) 中的密度梯度项为零，
不能据此声称支持一般变密度可压缩 SA。

### Wilcox1988m k-ω

`k_omega.cpp` 实现

\[
D_t(\rho k)=\nabla\cdot[(\mu+\sigma_k\mu_t)\nabla k]+P-\beta^*\rho k\omega,
\]
\[
D_t(\rho\omega)=\nabla\cdot[(\mu+\sigma_\omega\mu_t)\nabla\omega]
+\gamma\omega P/k-\beta\rho\omega^2,
\quad \mu_t=\rho k/\omega.
\]

这里 D_t 表示时间加守恒对流项。β*=.09、β=.075、γ=5/9、σk=σω=.5；注意扩散系数是
乘 σ，而 k-ε 是除 σ。左侧 sink 分别为 β*ρω 和 βρω。生产、耗散、涡黏度和扩散均按
[TMR Wilcox1988m 定义](https://tmbwg.github.io/turbmodels/wilcox.html) 核对。模型无交叉扩散、
SST 应力限制器或 2006 修正，壁面 ω 与 k 的具体近壁策略仍由 case 负责。

### 标准高 Re k-ε

`k_epsilon.cpp` 实现

\[
D_t(\rho k)=\nabla\cdot[(\mu+\mu_t/\sigma_k)\nabla k]+P-\rho\epsilon,
\]
\[
D_t(\rho\epsilon)=\nabla\cdot[(\mu+\mu_t/\sigma_\epsilon)\nabla\epsilon]
+C_1\epsilon P/k-C_2\rho\epsilon^2/k,
\quad\mu_t=C_\mu\rho k^2/\epsilon.
\]

Cμ=.09、C1=1.44、C2=1.92、σk=1、σε=1.3。左侧线性 sink 分别为 ρε/k、C2ρε/k；
更新系数后，两项与非线性耗散原式等价。独立测试检查源减 sink×未知量后的总反应，
避免只验证某种拆分写法。默认常数和源项对应[标准模型资料](https://doc.openfoam.com/2212/tools/processing/models/turbulence/ras/linear-evm/rtm/kEpsilon/)。
没有低 Re 阻尼、壁函数、可压缩或浮力项，不把给定零壁面通量的测试称为工程近壁策略。

## 4. 实际验收设计

`make test-rans` 构造独立 case 并保留临时目录及每次运行日志：

1. 三模型分别在 1/2/4 进程检查非零剪切生产、扩散率、总反应源和涡黏度；SA 额外检查
   q=0 边界派生系数，且有 q 梯度贡献。应力制造场单独检查 F1。
2. 完整 transientSimple 在 U=0、均匀变量条件下运行至 t=.1，dt=.01/.005/.0025，
   对照独立 ODE 解，检查时间步细化阶。k-ω：ω=1/(1+.075t)、k=(1+.075t)^(-1.2)；
   k-ε：k=(1+.92t)^(-1/.92)、ε=k^1.92；SA：单独 Python 反应式及 20000 步 RK4。
3. dt=.005 的三模型、两种时间格式还比较 1/2/4 进程字段，绝对差要求 <1e-9。
4. 四种选择 none/SA/kOmega/kEpsilon 的稳态反应扩散边值问题，1/2/4 进程必须终端
   `converged=true`，受测字段绝对差 <1e-8。
5. 故意设 kMin=epsilonMin=1 的衰减算例使裁剪保持变量不变。1/2/4 进程都必须返回 2，
   rTurb 明显非零且不得写成功最终结果，覆盖原 F3 假收敛。

这组试验验证公式、时间积分、非线性停止和并行实现。模型源项测试含非零剪切，
时间阶试验是均匀衰减，稳态是人为边值问题；尚无新增平板/通道/翼型的湍流摩阻、
速度剖面和网格收敛验收。强非定常压力耦合的完整 URANS 时间阶也不由均匀衰减证明。
因此可以接受本轮方程实现修复，不能据此把框架标为工程湍流全面验证完成。
