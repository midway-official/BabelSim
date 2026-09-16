# Physics 求解器编程模型

BabelSim Solver 是普通 C++ 函数，不是运行时对象。理想源码可以按下面的顺序阅读：

1. Case 加载文件场；
2. Parameters 读取物性和算法设置；
3. create*Field 创建程序场；
4. createEquation 绑定每个未知量；
5. 时间步或外层迭代开始；
6. Equation reset 后逐项 ddt/div/laplacian/source；
7. solve 写回未知场；
8. 用 math 计算显式量，用 diagnostics 读取残差和变化；
9. Solver 自己决定收敛并把 metric 交给通用 Reporter；
10. 选择字段并 write。

## 两个基本模板

### 稳态

~~~cpp
for (int iter = 0; iter < maxIterations; ++iter) {
    const auto previous = U;
    equation.reset();
    equ::div(equation, phi);
    equ::laplacian(equation, mu, -1);
    const double residual = diagnostics::relativeResidual(equation, U);
    const auto solved = equ::solve(equation, linear);
    // update, diagnostics, report, convergence
}
~~~

### 瞬态

~~~cpp
auto time = time::start(problem);
auto history = time::history(T);
while (time.value() < time.end()) {
    time.advance();
    history.save(T, time.dt());
    equation.reset();
    equ::ddt(equation, capacity, history);
    // other terms, solve, report, write
}
~~~

瞬态 SIMPLE 在自己的 main 中再嵌套一个 SIMPLE 外迭代；稳态 SIMPLE 只有一个主要外循环。
压力 non-orthogonal 修正是可选子循环。两个 SIMPLE 不共享算法文件，Rhie–Chow 只在各自
main 用公开 math 操作组合。

## 层次边界

Case/Field 管数据和生命周期，math 管显式场数学，equ 管离散 Equation，time 管物理时间和
History，diagnostics 管只读度量，monitor 管格式化报告，Runtime/后端管通信和代数。
Physics 不读 token、不碰 MPI/Eigen/CSR、不访问原始数组。RANS 模型各自实现闭合，通过
effectiveViscosity、deviatoricStressRemainder 和 solveTransport 交互。

请用 [dsl-solver-guide.md](dsl-solver-guide.md) 获取完整代码片段和测试清单。
