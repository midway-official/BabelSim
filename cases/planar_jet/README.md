# 无湍流模型的二维平面射流

本算例使用瞬态不可压层流方程和 PISO，不启用 RANS/LES 湍流闭合，也不求解湍流输运变量。三维网格只有一个展向单元，前后面设为对称边界，以表示二维平面流动。

## 物理与数值设置

- 喷口宽度 D=1，中心速度 Uj=1，密度 rho=1，动力黏度 mu=0.001，因此 Re_D=1000。
- 计算域 x/D=0..20、y/D=-5..5；正交网格 500 x 200 x 1，共 100,000 个六面体，喷口宽度内 20 个单元。
- deltaT=0.01 D/Uj，endTime=60 D/Uj，共 6,000 步，最大入口 Courant 数约为 0.25；BDF2 首步按求解器约定以 Euler 启动。
- 对流项使用 linearUpwind（基于最小二乘梯度的二阶迎风重构），时间项使用 BDF2；PISO 每步做 3 次压力修正。
- 稳态入口速度是平滑顶帽分布。初始速度场包含幅值 0.01 Uj 的横向小扰动，用于触发剪切层不稳定性；没有时间变化的入口表达式。
- 每 100 步（1 D/Uj）写出 U 和 p。运行后可执行 python3 cases/planar_jet/plot_flow.py 生成末时刻涡量图、速度剖面和探针历史。

## 重建与运行

1. 执行 python3 cases/planar_jet/generate_mesh.py
2. 执行 mpirun -n 2 build-petsc/babelsim-solve -case cases/planar_jet
3. 执行 python3 cases/planar_jet/plot_flow.py
4. 执行 python3 cases/planar_jet/make_vorticity_gif.py，生成 60 帧、固定色标的涡量动画。

网格生成器同时写出 mesh/planar_jet.mesh、结构化网格几何报告和初始速度文件 fields/initial/U.initial.dat。二维涡量定义为 omega_z = dUy/dx - dUx/dy。
