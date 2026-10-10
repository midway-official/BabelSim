# 10 详细操作示例

[用户手册](README.md) · 上一章：[结果与继续计算](results-and-restart.md) · 下一章：[故障排查](troubleshooting.md)

本章各示例由同一个生成器创建，均包含完整输入文件，目录存在时拒绝覆盖。默认使用小网格便于学习软件操作，不代替独立物理验证。热传导的逐文件解释见[首次运行](quick-start.md)，完整生成逻辑见[create_case.py](../examples/create_case.py)。

## 示例 A：给定速度的标量输运

### 问题与输入

128 单元薄层域，初始浓度 C=0、入口 C=1、出口零梯度。给定 U=(1,0,0)，storage=1、diffusivity=0.01、source=0；Euler 推进 5 步到 0.05。求解器只计算 C，不求解 U/p。

```bash
python3 docs/examples/create_case.py --template transport --output work/transport
build/babelsim-solve -case work/transport -time serial
mpirun -np 2 build/babelsim-solve -case work/transport -time mpi2
build/babelsim-post -case work/transport -time serial/all -format vtk
python3 tools/compare_parallel_results.py \
  work/transport/results/serial/0.05 work/transport/results/mpi2/0.05 \
  --atol 1e-8 --rtol 1e-8
```

### 检查与修改

应保存 C 和给定 U，输出 5 个物理时刻，串并行比较通过。观察 C 从入口向域内传播，明确最终时间很短，并不代表整个域达到稳态。将 `equation.transport.convection` 从 upwind 改为 linearUpwind，使用新的运行名对比过渡区域与过冲；保持其余输入不变。

若想求解由 C 反作用于 U 的耦合问题，不能仅修改 solver 的一个参数，应开发新的耦合算法。

## 示例 B：稳态通道流的求解器调用

### 问题与输入

256 个六面体单元，域 4×1×0.1，入口 U=(1,0,0)，上下壁面无滑移，出口 p=0，front/back 使用对称设置；rho=1、mu=0.01、turbulenceModel=none。methods 使用 steady，SIMPLE 上限 1200，速度/压力松弛为 0.5/0.3。

```bash
python3 docs/examples/create_case.py --template simple --output work/channel
build/babelsim-solve -case work/channel -time serial
mpirun -np 4 build/babelsim-solve -case work/channel -time mpi4
build/babelsim-post -case work/channel -time serial -format vtk tecplot
python3 tools/compare_parallel_results.py \
  work/channel/results/serial work/channel/results/mpi4 \
  --atol 1e-7 --rtol 1e-7
```

### 预期与证据范围

迭代日志应出现收敛结束状态、退出码 0，并保存 U/p。稳态结果数值时刻为 startTime=0，并非 control 的 endTime。固定入口和无滑移壁面产生入口发展过程，此小域/粗网格示例不能直接当作完全发展 Poiseuille 解析解精度证明。

更严格的通道和方腔验证入口是 `make validate-poiseuille`、`make validate-cavity`。它们的参考定义与范围见[验证指南](../developer-guide/verification.md)，不要将快速教程与这两个物理基准混同。

## 示例 C：PISO 瞬态通道

```bash
python3 docs/examples/create_case.py --template piso --output work/piso
build/babelsim-solve -case work/piso -time serial
mpirun -np 4 build/babelsim-solve -case work/piso -time mpi4
build/babelsim-post -case work/piso -time serial/all -format vtk
```

几何与层流物性同通道示例，methods 使用 Euler、dt=0.01，5 步到 0.05。关键 solution 设置为：

```text
maxIterations 1
nCorrectors 2
maxCorrectors 4
velocityRelaxation 1
continuityTolerance 1e-8
```

这里外层 maxIterations 与各 `equation.*.maxIterations` 的线性上限不同。示例未设置严格的 couplingTolerance，默认极大阈值允许观察入门运行；检查 `rCoupling`、`CoMax` 和质量指标后，再为研究问题明确配置耦合阈值并做时间步研究。

如果把 `couplingTolerance` 收紧，最大校正次数仍可能不足，此时失败应作为数值诊断依据。增加 maxCorrectors 或减小 dt 后需重新对比结果，不能只屏蔽退出码。

## 示例 D：程序调用与二次开发

```bash
make -f docs/examples/Makefile -j4
python3 docs/examples/create_case.py --template custom-heat --output work/customHeat
build/examples/custom-heat -case work/customHeat -time external
python3 docs/examples/create_case.py --template coupled --output work/coupled
build/examples/coupled-scalar -case work/coupled -time serial
```

热示例说明注册函数如何由 case 调用；耦合示例说明两个方程、显式物理历史、内层耦合与派生 total 场。完整源码逐段解释、外部工程链接和内置模块接入见[新增与调用求解器](../developer-guide/extending-solvers.md)。

## 示例 E：自动验证手册示例

```bash
python3 docs/examples/verify_examples.py --output build/docs-refresh
```

脚本从根目录构建示例，顺序执行六类算例的 1/2/4 进程、串并行比较、热方程 Euler/BDF2 比较、注册分派、VTK/PVD/Tecplot、Euler 场继续计算、单单元网格读取与错误路径。日志、生成输入和 JSON 报告保留在输出目录中。

执行范围是教程运行、公开 API 链接、并行与 I/O；不是所有生产测试，更不是工程物理认证。测试报告记录实际命令、退出码和文件位置，便于交付时复查。
