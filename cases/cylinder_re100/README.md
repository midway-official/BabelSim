# 二维 Re=100 圆柱绕流：20,486 三角柱网格

D=1、U∞=1、rho=1、mu=0.01；层流，PISO，单层挤出厚度 0.1D。计算域 x=[-10,25]D，y=[-10,10]D。圆柱壁无滑移，入口 (1,0,0)，出口 p=0/U 零梯度；上下远场和前后面为 symmetry。

Python 构建 160 周向段、12 层三角形环带（首层 0.015D、增长 1.10），Gmsh 为外部区域生成非结构三角网格，尾流 x=[0,15]D、y=[-2,2]D 目标尺寸 0.11D。Python 合并共形接口、挤出三角柱、写 BABELSIM_MESH 3，并验证体积/内角/非正交/偏斜/拓扑及边界覆盖。脚本依赖 numpy、matplotlib 和 gmsh 可执行程序。

```bash
python3 cases/cylinder_re100/generate_mesh.py
```

该命令会重建网格与初场，已有计算运行期间不要重建。mesh/mesh_quality.json 是生成端检查，mesh/reader_quality.json 是原生读取器检查。check_mesh.cpp 用 libbabelsim.a 和 PETSc 库链接后运行；mesh/mesh.png 是已检查的网格图。

初场基于解析无散流函数，近壁速度平滑衰减到零；仅在圆柱下游有限区域加入一次性小幅流函数扰动，入口没有持续周期强迫。初始压力是启动估计，不是稳态或精确瞬态解。

当前保守配置为 BDF2、linearUpwind、leastSquares、corrected、dt=0.001、alphaU=1，3 次压力校正和 2 次非正交校正。dt=0.005 的原始启动试验失败；同网格 dt=0.001 已通过 1000 步/ t=1 的短算。详见 validation/analysis_and_plan.md，不能把该短算等同于长时间稳定或涡脱落验证。

```bash
python3 cases/cylinder_re100/run.py
```

启动本地后台 worker，mpirun 为 2 ranks/每 rank 1 线程并绑定 CPU core。目标 t=150、每 0.25 保存 U/p，结果写入 results/run。run.py 会拒绝覆盖既有结果。results/run_status.json 记录 worker/mpirun PID 和状态，results/solver.log 保存求解日志。守护检查只拦截明显发散，不能验证物理正确性。求解成功退出后自动运行后处理并生成 GIF。

```bash
python3 cases/cylinder_re100/postprocess.py --dt .001 --gif
```

后处理按 globalCellId 合并两个 rank；使用加权最小二乘梯度计算带符号涡量，输出尾流探针、压力加一层切向黏性牵引估算的升阻力和有限时窗频谱。梯度通过仿射场精确性检查。应结合涡量正负结构、探针周期性以及网格/时间步敏感性评估涡脱落。升阻力为近壁一层梯度估计，尚未做力积分精度验证。

原算例创建阶段的数值失稳分析、对照实验和文献方案位于本地 validation/analysis_and_plan.md；其中结论对应当时的源码与输入。该算例的网格、脚本和配置作为 Re3900 生成依赖随本次提交保存，results 和大规模原始验证输出保留本地。后续通用求解器优化与验收见 [圆柱优化报告](../../docs/reports/cylinder-optimization-2026-10-09.md)。
