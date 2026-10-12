# 03 新增与调用求解器

[开发指南](README.md) · 上一章：[公开接口](public-api.md) · 下一章：[验证与回归](verification.md)

本章完成三个不同层次的操作：调用已提供的求解器；将自写求解器链接成独立程序；将同一算法注册进正式应用。最后用双场耦合示例说明时间推进与二次开发。

## 1. 调用现有求解器

```bash
python3 docs/examples/create_case.py --template heat --output work/callHeat
build/babelsim-solve -case work/callHeat -time firstCall
```

命令行选择 case，case 中 `solver heat` 选择已注册的函数。改变算例中的密度、热容、导热率或源项不需要重编译；新增数学方程或算法步骤需要 C++ 开发。

## 2. 外部热传导求解器

完整源码为 [external_heat.cpp](../examples/external_heat.cpp)，应用入口为 [external_main.cpp](../examples/external_main.cpp)。它们只包含公开头文件。

先编译全部开发示例：

```bash
make -f docs/examples/Makefile -j4
python3 docs/examples/create_case.py --template custom-heat --output work/customHeat
build/examples/custom-heat -case work/customHeat -time external
mpirun -np 4 build/examples/custom-heat -case work/customHeat -time mpi4
```

输出分别位于 `results/external`、`results/mpi4`。示例名为 `customHeat`，由生成器写入 case；如果在外部程序中仍选 `solver heat`，因为未链接正式 heat 注册项，会报未知求解器。

### 装配过程解释

求解函数首先加载 T 并读取类型化物性，在循环外创建 temperature 方程、时间对象和显式历史；随后校验所有已配置参数均被使用。每一步按如下顺序执行：

```cpp
auto step = time::beginStep(clock, histories);
// BDF2 从第二步开始，完整源码中先外推空间修正使用的 T。
equation.reset();
equ::ddt(equation, rho * cp, history);
equ::laplacian(equation, k, -1, "diffusion");
equ::source(equation, Q);
const auto result = equ::solve(equation);
if (!result.converged()) return SolverResult{result.status};
step.accept();
if (clock.step() % interval == 0 || clock.finished()) write(problem, clock);
```

`histories.add(T, history)` 在循环前登记关联历史；未接受的 Step 析构时自动恢复时钟、场和历史。完整时间契约见[时间推进与关联历史](time-integration.md)。

Equation 绑定未知 T，reset 防止上一时间步的贡献叠加；History 在物理步开始保存，不在求解迭代中推进；-1 对应扩散项放左端；写出发生在接受步骤后。完整示例支持与内置热求解器一致的 Euler/BDF2 路径。

注册和入口分别是：

```cpp
const babelsim::SolverRegistration registration("customHeat", runCustomHeat);
```

```cpp
int main(int argc, char** argv) {
    return babelsim::runApplication(argc, argv, reportError);
}
```

错误回调打印消息；应用入口维护 MPI 生命周期与退出码。不要在每个物理求解函数里重复初始化或结束 MPI。

### 外部工程如何链接

示例 Makefile 采用两阶段构建：

```text
g++ -std=c++17 -Iinclude -c <solver-source.cpp> -o solver.o
mpic++ solver.o main.o build/libbabelsim.a <PETSc 与可选 METIS 链接参数> -o mySolver
```

这两行是结构说明，其中尖括号项必须替换。实际、可直接执行且处理当前依赖的命令是 `make -f docs/examples/Makefile`。其中源文件编译只用公开 include；链接参数从同一根 Makefile 获取，避免省略静态框架所需的 MPI、PETSc 和已启用的 METIS。

仓库外工程可复制 `include/` 与 `libbabelsim.a`，并提供匹配的依赖安装。PETSc 参数由该安装的 pkg-config 获得，加入其 lib 动态搜索路径；若框架构建启用了 METIS，还必须链接该库。需要完全避免 METIS 运行依赖时先用 `USE_METIS=0` 构建框架。不同编译器、MPI 或 CPU 优化下的库不是自动兼容的 SDK。

可直接创建并编译仅含公开 SDK 的独立工程：

```bash
make -f docs/examples/Makefile sdk-example
build/examples/sdk-project/custom-heat -case work/customHeat -time sdk
```

生成的工程 Makefile 见[external_project.mk](../examples/external_project.mk)。工程内部有公开头、框架静态库和自己的两个源文件，编译阶段不依赖 src、Eigen 或 MPI/PETSc 头；仍需链接与运行匹配的 MPI/PETSc 库。

## 3. 作为内置求解器注册

热示例源文件不包含 main，可直接作为一个新的内置注册模块：

```bash
mkdir -p src/physics/custom_heat
cp docs/examples/external_heat.cpp src/physics/custom_heat/main.cpp
make -j4
build/babelsim-solve -case work/customHeat -time builtIn
```

只在准备新增该模块时执行这组源码修改命令。根 Makefile 自动收集 `src/physics/*/main.cpp`，直接链接注册对象。不要把注册对象仅放进静态库等待按需抽取；也不要在集中启动器中手工维护求解器名称名单。

手册还提供 `build/examples/registered-heat`，已将同一注册源与正式求解器对象一起链接，可以先验证注册效果：

```bash
build/examples/registered-heat -case work/customHeat -time registered
build/examples/registered-heat -case work/callHeat -time registered
```

这验证一个程序同时分派 customHeat 与正式 heat。示例 Makefile 对按上述路径安装的模块作排除，避免再次加入同一注册；其他目录里的同名注册仍需由开发者消除。

## 4. 二次开发：双场耦合

[coupled_scalar.cpp](../examples/coupled_scalar.cpp)求解：

```text
dT/dt - div(D*grad(T)) + a*T = a*C
dC/dt - div(D*grad(C)) + a*C = a*T
total = T + C
```

热侧 T 的边界输入通过反应交换项影响 C。两个未知量分别使用具名方程 temperature 和 concentration。程序创建派生场 total 并选择输出，用户无需提供 total.field。

```bash
python3 docs/examples/create_case.py --template coupled --output work/coupled
build/examples/coupled-scalar -case work/coupled -time serial
mpirun -np 4 build/examples/coupled-scalar -case work/coupled -time mpi4
build/babelsim-post -case work/coupled -time serial/all -format vtk
```

关键结构是：每个物理步只保存一次旧 T/C；内层重复装配并顺序求解 T、C，用整场变化判断耦合；历史不在内层更新；未收敛返回结构化失败；耦合收敛后更新 total，接受步骤后写出。示例要求 Euler，不能只修改 time bdf2 就宣称得到经过验证的二阶耦合算法。

`couplingIterations` 和 `couplingTolerance` 是该示例自己的顶层 solution 参数。两个方程都需要完整空间选项与线性配置。这样新增模型不需要让通用后端认识新的物理方程。

## 5. 向新的 PDE 扩展

先写清未知量、容量、源项、通量方向和边界条件，再将各项映射到 ddt/div/laplacian/reaction/source。变系数可用 ScalarField，显式右端可用整场 math 运算。非线性系数改变后重新装配；重复求解已冻结的 Equation 不会自动更新它。

新增方程名必须对应新的 methods 与 solution 配置；名称不能绑定另一未知场。先实现小型解析/制造解，再增加成功与失败状态、时间阶、边界、并行一致性及实际输出检查。不要仅用“编译成功”和“某个案例残差下降”作为新算法的验收结论。
