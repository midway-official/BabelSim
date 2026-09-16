# RANS 模型开发与维护

RANS 模型属于 Physics 模块。每个模型独立实现自己的输运方程、历史、边界、闭合常数和
源项；模型之间不共享一套隐藏的闭合公式。当前模型包括 Spalart–Allmaras、k–omega 和
k–epsilon。

## 与 SIMPLE 的交互

SIMPLE 只看到具有物理意义的对象：

~~~cpp
auto turbulence = rans::load(problem, U, phi);
const auto& muEff = turbulence.effectiveViscosity();

equ::laplacian(momentumEquation, muEff, -1);
equ::source(momentumEquation,
            math::div(turbulence.deviatoricStressRemainder(U)));

const auto transport = turbulence.solveTransport();
~~~

effectiveViscosity 是动量方程的有效黏度。deviatoricStressRemainder 是在隐式黏性项
之外需要显式加入的应力余项；不要把完整应力重复加入。solveTransport 只返回每个模型
自己的线性结果、初始残差和相对变化，SIMPLE 决定如何把它们纳入外层收敛。

## 配置分类

physics 字典放 density、dynamicViscosity、turbulenceModel 和模型常数。solution 字典
放 turbulenceRelaxation、turbulenceTolerance、线性求解器及算法迭代上限。通过
Parameters::word/positive/fraction 等接口读取，模型不访问 entry().tokens。

## 新模型的边界

新模型在 src/physics/RANS/ 中拥有独立源文件和字段。它可以使用公共 Case、Field、math、
equ、time、diagnostics 和 solver API，但不能包含 SIMPLE 私有头、src/internal、MPI、Eigen、
CSR 或原始 Field 存储。若需要新增交互，优先增加具有物理含义的小型值类型或成员函数，
不要把具体模型公式移入通用离散层。

模型必须在至少一个制造或解析案例中验证输运残差、正值限制、边界类型和 1/2/4 rank 一致性。
不要用统一的模型 Manager 把 SA、k–omega 和 k–epsilon 混成一个实现。
