# 审查证据

对应报告：[physics-audit-2026-09-26.md](/home/midway/BabelSim/docs/reports/physics-audit-2026-09-26.md)。基线 `005c8181222c451af68b94da277b162a9aea6a35`，PETSc 3.25.5，仓库默认 release 库。`make -q build-petsc/babelsim-solve build-petsc/rans_equations_test` 确认两项生产/模型测试可执行文件与当前依赖一致。没有修改求解器源码。

| 文件 | 内容 |
|---|---|
| operator_probes.cpp / .log | F1：扩散矩阵—面通量恒等式；F2：inletOutlet 已有回流上下文 |
| piso_decay_probe.py / .log / piso_decay.json | F3：生产 PISO，三个模型，松弛 0.7/1，多个 dt，解析衰减比较 |
| piso_failure_probe.py / .log | F4：生产 PISO 接受未收敛的 k 线性方程 |
| piso-failure-performance/rank-0000.json | F4 逐方程 PETSc 状态及真残差 |
| time_order_probe.cpp / time_order.log | F5：与 heat/transport 相同的公共 DSL 装配顺序，固定网格时间阶反例及对照 |
| unit_tests.log | 六项现有测试全部通过 |
| rans_validation.log | 现有三模型验证，包括 1/2/4 rank、一/二阶时间衰减与裁剪拒绝 |

C++ 反例通过公共 DSL 运行；仅在观测实际单元/面数值时使用测试内部读取接口。PISO 两项脚本复用现有 `rans_validation_test.py` 的案例生成函数，但不执行原脚本后半部分测试循环；实际调用生产 `babelsim-solve`。生成的案例存放在 `/tmp/babelsim-rans-validation-*`，位置记录在日志中；报告所需汇总、程序和失败状态已保存于本目录，不依赖临时目录才能阅读。

从仓库根目录复现（根据实际 PETSc 安装修改路径）：

```bash
make -j4 build-petsc/babelsim-solve build-petsc/rans_equations_test
export PETSC_DIR=/home/midway/opt/petsc-3.25.5-openmpi-opt
for probe in operator_probes time_order_probe; do
    mpic++ -std=c++17 -O0 -Iinclude -Isrc -Itests -I/usr/include/eigen3 \
        -I"$PETSC_DIR/include" \
        "docs/reports/physics-audit-2026-09-26-evidence/$probe.cpp" \
        build-petsc/libbabelsim.a -L"$PETSC_DIR/lib" \
        -Wl,-rpath,"$PETSC_DIR/lib" -lpetsc -o "/tmp/babelsim-audit-$probe"
    "/tmp/babelsim-audit-$probe"
done
python3 docs/reports/physics-audit-2026-09-26-evidence/piso_decay_probe.py
python3 docs/reports/physics-audit-2026-09-26-evidence/piso_failure_probe.py
```

现有测试的实际运行命令：

```bash
make -j4 build-petsc/procedural_equation_test build-petsc/numerical_contract_test \
    build-petsc/math_runtime_test build-petsc/operators_test \
    build-petsc/time_history_test build-petsc/field_boundary_test
for t in procedural_equation numerical_contract math_runtime operators time_history field_boundary; do
    "build-petsc/${t}_test" || exit 1
done
python3 tests/rans_validation_test.py --solver build-petsc/babelsim-solve \
    --equations build-petsc/rans_equations_test
```

F1/F2/F5 反例未在 MPI 模式运行。F5 是时间分辨率自收敛试验，不以空间连续解析解作比较；无需把正交方案在剪切网格上的空间偏差当成时间误差。F3 的 k/omega/epsilon 为解析参考，SA 使用独立小步长 RK4 参考。程序打印的问题是待修复行为，不是将错误结果定义为正确回归期望。
