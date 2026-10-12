# 时间推进与关联历史

[开发指南](README.md) · [公开接口](public-api.md) · [外部求解器示例](../examples/external_heat.cpp)

时间层、物理迭代与步骤接受由求解器明确控制。公共工具处理时间系数、同步历史、预测和可选回滚；`equ::ddt` 的装配公式及后端实现保持原有行为。

## 默认推进：保存关联历史

```cpp
#include "babelsim/time_discretization.h"

time::HistoryGroup histories;
auto& oldT = histories.addHistory(T);
auto& oldC = histories.addHistory(C);
auto clock = time::makeStepper(problem);
problem.finalizeInitialization();
while (!clock.finished()) {
    clock.advance();
    histories.save(clock);
    // 同一时间步的所有非线性/耦合迭代使用这些历史，不再 save。
    equation.reset();
    equ::ddt(equation, capacity, oldT);
    // 装配、求解及接受判断……
}
```

`HistoryGroup::addHistory(field)` 创建并持有历史，返回稳定引用；容器增长不使该引用失效。字段仍由 Case 持有，Case 必须比历史组存活更久。初值修正、投影和模型初始化完成后再调用 addHistory，明确捕获最终初始状态。专家算法也可用 `add(field, existingHistory)` 关联自行持有的历史，但须自行保证它的生存期。一个字段和历史只能登记一次，所有字段属于同一网格；开始推进后关闭登记。

`save(clock)` 检查时钟、网格、连续的步号，以及 `levels/dt/previousDt`。重复保存、跳过物理步、手工推进其中一个历史均会失败。没有 TimeStepper 的高级算法可以使用 `save(dt)`，但绑定时钟后不允许混用两种路径。错误的步长/不同步在保存任何成员之前被拒绝。

默认路径只执行原有历史复制和小规模元数据检查。`trackForRollback()` / `watch()` 登记的回滚状态在此路径不复制；只有显式创建事务或检查点才生成快照。

## 时间系数与显式预测

```cpp
const auto coefficients = time::weights(problem.read.methods().time, oldT);
// current*x_new = previous*x_old + older*x_older （时间项的 Ax=b 形式）
```

Euler：`current=previous=1/dt, older=0`。BDF2 第一次保存只有一层历史时使用 Euler；此后设 `r=dt/previousDt`：

```text
current  = (1+2*r)/((1+r)*dt)
previous = (1+r)/dt
older    = -r*r/((1+r)*dt)
```

这是固定控制体上现有非保守容量时间项的公共描述。测试将公共权重与原 `equ::ddt` 的积分对角和右端逐项比较，底层装配仍调用现有实现。

显式时间预测只在算法要求的位置调用：

```cpp
if (problem.read.methods().time == TimeMethod::BDF2 && oldT.levels() >= 2)
    time::predict(T, oldT, problem.read.methods().time);
```

BDF2 使用 `(1+r)*previous-r*older`，保留目标未知量的物理边界。Euler/启动分支直接恢复保存的 previous，稳态分支不改变目标。内置热传导和输运仅在 BDF2 两层就绪时预测，因此 Euler 没有新增全场遍历。流动求解器继续决定自己的预测/校正顺序，不自动外推全部未知量。

PISO 的时间通量修正保留在物理层，仅通过 `time::requireMatching(U_old, phi_old)` 和公共权重取得旧时间层系数。

## 可选的时间步事务

```cpp
#include "babelsim/time_step.h"

time::HistoryGroup histories;
auto& oldT = histories.addHistory(T);
histories.trackForRollback(pressure, derivedCoefficient); // 没有时间导数，但拒绝步骤时必须恢复
histories.watch(modelState);         // 可复制的模型计数器/状态结构
auto clock = time::makeStepper(problem);
problem.finalizeInitialization();
while (!clock.finished()) {
    auto step = time::beginStep(clock, histories);
    // beginStep 已推进时钟、保存所有关联历史。
    equation.reset();
    equ::ddt(equation, capacity, oldT);
    const auto result = equ::solve(equation);
    if (!result.converged()) {
        step.reject();
        // 调用者可以选择更小的 dt 再试；或返回失败。
        return SolverResult{result.status};
    }
    step.accept();
    write(problem, clock); // 外部输出在接受之后发生
}
```

`Step` 只能移动，不能复制。析构时未接受的步骤自动拒绝。拒绝恢复时钟 `value/step/dt`、完整字段值/边界模式/计算边界迹/边界通量上下文、两层历史与其元数据，以及登记的模型状态。`Field::restoreSnapshot()` 专门提供完整状态恢复，普通赋值仍保持原有仅更新值的语义。

`beginStep(clock, histories, requestedDt)` 支持调用者指定步长，末步截断到 endTime。接受显式短步后，下一次默认步继续从当前时间推进 configured deltaT；拒绝则恢复原推进模式。它提供重试基础设施，尚不包含误差估计、步长控制器或 BDF2 步长比稳定性策略。

所有需回滚的场都必须通过 `addHistory/add/trackForRollback` 登记，模型自己的可复制状态通过 `watch` 登记。Field/History 不能传入 `watch`，避免值赋值误代替完整恢复。自定义状态赋值应可靠且不抛异常；不透明外部资源、日志、文件输出和已执行的外部副作用不属于内存回滚。启动时声明阶段关闭不撤销，必须在开始时间步前完成参数读取和场声明。

待接受步骤期间不能再推进时钟或保存历史，也不能嵌套事务。与历史组无关的其他状态不会自动发现。SIMPLE/PISO/RANS 通过模型 `registerTimeState` 将输运量登记到同一组；`beginStep` 钩子只重置每步模型计数器，不能推进历史。模型原来的 `saveOld` 路径已删除，模型和宿主不再分别保存同一步历史。

## 容量项的适用范围

现有接口离散 `c*d(x)/dt`，容量场取本次装配值。它不等于 `d(c*x)/dt`。公共权重、历史组和事务没有增加保守乘积时间项；时变密度/容量问题需要另外保存系数历史并验证守恒，不能只将常数系数替换为字段。

## 验证入口

```bash
make test
make test-scalar-time
make test-piso-coupling
make test-rans
make test-external
```

`time_tools_test` 覆盖独立时间公式、原装配一致性、预测、错误历史和完整回滚/重试；生产标量时间阶、PISO 耦合和 RANS 方程检查分别验证实际算法迁移。
