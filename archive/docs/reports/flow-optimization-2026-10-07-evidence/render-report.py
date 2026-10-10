#!/usr/bin/env python3
"""Render the completed report from the measured summary and field checks."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NAME = ROOT.name


def main():
    summary = json.loads((ROOT/'final-summary.json').read_text())
    checks = json.loads((ROOT/'field-comparisons.json').read_text())
    assert checks and all(check['passed'] for check in checks)
    rows = summary['scaling']
    table = '| 求解器 | MPI 进程 | 样本数 | 纯计算中位数 / 秒 | 加速比 | 标准并行效率 | 峰值 RSS / GiB | 峰值 PSS / GiB |\n'
    table += '|---|---:|---:|---:|---:|---:|---:|---:|\n'
    for r in rows:
        table += f"| {r['case'].upper()} | {r['ranks']} | {r['samples']} | {r['computeMedian']:.2f} | {r['speedup']:.2f}× | {100*r['efficiency']:.1f}% | {r['peakRssGiB']:.2f} | {r['peakPssGiB']:.2f} |\n"
    comparison = '| 求解器 / 进程 | 内存优化前 / 秒 | 最终 / 秒 | 耗时变化 | RSS 变化 | PSS 变化 | 原 / 新进程交换页峰值 / MiB |\n'
    comparison += '|---|---:|---:|---:|---:|---:|---:|\n'
    for r in summary['comparisons']:
        comparison += f"| {r['case'].upper()} / {r['ranks']} | {r['referenceSeconds']:.2f} | {r['finalSeconds']:.2f} | {100*r['computeChange']:+.1f}% | {100*r['rssChange']:+.1f}% | {100*r['pssChange']:+.1f}% | {r['referenceSwapMiB']:.1f} / {r['finalSwapMiB']:.1f} |\n"
    field_table = '| 对照 | 字段 | 全局单元数 | 最大绝对差 |\n|---|---|---:|---:|\n'
    for check in checks:
        for field, result in check['fields'].items():
            field_table += f"| {check['label']} | {field} | {result['cells']:,} | {result['maxAbsDifference']:.3g} |\n"
    conclusions = []
    for case in ('simple','piso'):
        selected = [r for r in rows if r['case']==case]
        twelve = next(r for r in selected if r['ranks']==12)
        best = min(selected, key=lambda r:r['computeMedian'])
        conclusions.append(f"{case.upper()} 的 12 进程加速比为 {twelve['speedup']:.2f} 倍、效率 {100*twelve['efficiency']:.1f}%；本次测得最短时间出现在 {best['ranks']} 进程（{best['computeMedian']:.2f} 秒）。")
    text = f'''# BabelSim 通用速度、并行与内存优化：2026-10-07

## 当前测量结果

{' '.join(conclusions)}这些是同一最终二进制、相同输入及固定工作量的强扩展结果。
网格为 500,094 个扭曲三角棱柱，1,266,111 个面、262,144 个顶点，3 层 halo。
这轮实现与验证完成；结果适用于下面所列硬件与算例，不能直接外推到多节点网络或其他 AMG 难度。

{table}
12 进程数据取独立重测样本的中位数；其他进程数只有一个持续样本。
原始最小/最大值、样本数见 [最终汇总 JSON]({NAME}/final-summary.json)。
单样本点不能给出统计置信区间，也不能证明小于几百分点的差异稳定存在。

![同版强扩展]({NAME}/final-scaling.png)

## 测量契约

- AMD Ryzen AI MAX+ 392，12 个物理核、24 个逻辑线程，单 NUMA；WSL2 分配约 7.6 GiB RAM、2 GiB swap。
- 每个 MPI rank 单线程，OpenMPI 4.1.2 绑定物理核；PETSc 3.25.5、METIS 5.1.0、GCC 11。
- 所有测试逐项执行。性能窗口中没有并行启动另一个测试、编译、Python 模型或绘图任务。
- SIMPLE 固定 128 次外迭代，计量后 127 次间隔。退出 `maxIterations` 是吞吐实验预期状态，
  **这不是达到稳态物理收敛所需的时间**；原有收敛 SIMPLE 回归另外检查。
- PISO 使用同一隐式 Euler 输入、`dt=0.01`，接受 6 步至 `t=0.06`，计量后 5 个间隔。
  保持原耦合判据与线性容差，要求日志报告接受、收敛及线性检查正常。
- `T_P=max_rank(iterationComputeSeconds)`。首次迭代前的网格读取、Mesh/场/方程创建、初次稀疏模式
  和预条件器建立，以及完整首次迭代均排除；扣除窗口内结果写出。每个正式样本的该窗口都超过 60 秒。
- 强扩展保持工作量一致，所以单进程样本超过 3 分钟。`S_P=T_1/T_P`，`E_P=S_P/P`；
  4→12 的增量效率是 `T_4/(3T_12)`，与标准效率不同。

总 RSS 每 0.1 秒采样，PSS 每秒采样。RSS 累计共享页时可能重复计数，PSS 按比例分摊。
表中内存峰值覆盖完整进程生命周期，包括初始化与输出；多样本取这些峰值的中位数。
它们不是纯计算窗口内的内存峰值。进程 VmSwap 与系统换页记录单独保存。

## 内存与划分相对于前一速度优化版本的取舍

这里的参考版本已包含预条件器复用和 RHS 工作量保护；它不是最初 Git HEAD。
输入哈希与最终版本相同。4 进程参考各一个样本；12 进程参考和最终版本各重测并取中位数。

{comparison}
![完整运行内存峰值]({NAME}/final-memory.png)

12 进程参考版本触发了本机有限 RAM 下的进程交换页，新版的对应值由表列出。
初始化已从计算计时排除，但初始化造成的换页可能延续到计算阶段。
因此 12 进程对照反映这台资源受限机器上的实际表现；4 进程对照另用于判断无换页时的代价。
用户允许少量速度损失以降低内存占用，报告保留速度回退项，不以省内存掩盖它。
先前游戏干扰的数据保留在 [排除记录]({NAME}/engineering-screen/excluded.json)，不进入上述汇总。

## 通用实现与抽象

Physics 源码没有改动，物理程序继续通过既有 `equ`、`math`、场和 Case 接口组织计算。
方程离散、压力/速度松弛、物理边界、线性真实残差目标和步接受规则保持各测试的原设置。
改动集中在 Backend、Parallel 和运行观测层：

1. 精确比较已装配矩阵数值；任一 rank 有变化就共同更新，完全不变时只更新 RHS。
   同一方程保留 Mat/KSP/PC 和工作数组，冷启动利用 `r0=b` 省去多余 SpMV，warm start 仍计算 `b-Ax0`。
2. 预条件器默认复用预算为 8 个矩阵版本及 8 个有实际迭代工作的 RHS。
   完全不变矩阵的 RHS 不单独触发周期重建；迭代增长请求重建，旧 PC 失败时以原始初值和同一残差目标重试一次。
   Block Jacobi 局部因子准备计入 setup 阶段，计数保留恢复尝试。
3. COO 行列仅在建立矩阵模式时临时存在，长期计划只保存 upper/lower 面索引，owned 顺序复用 Mesh。
   数值快照交换缓冲避免整段复制；三个向量分量顺序共用一对 RHS/解 Vec，分别保留各分量收敛记录。
4. rank 0 读取完整网格并计算一次所有权，逐 rank 分发 owned 与所需 halo 的面环、点坐标和全局元数据。
   工作 rank 只构造局部 Mesh；root 与工作 rank 在局部几何构造前释放完整网格或传输缓冲。
   设置阶段错误在进入下一阶段前集体传播，保持物理 patch、面方向、全局 ID 和 stencil 契约。
5. 可选 METIS 采用通信量目标、单元数/面关联数双约束、连通输入的连续分区，随后恢复精确配额。
   没有 METIS 时保留确定性图生长实现。公开 SDK 不引入 METIS 类型。

预条件器选项的运行行为见 [数值配置](../numerical-configuration.md)，
后端生命周期和计数边界见 [性能方法](../performance/README.md)。

## 先建模型，再接入框架

先研究 13 类符合面连接网格契约的图，包含不同长宽比、旋转、障碍、分支、局部加密、
三角棱柱、四面体、变面数 Voronoi、不连通域和打乱编号。
最终候选比较为 52 组问题、572 次顺序实验。
所选固定种子 0 的 `metis-volume-contig` 三层 ghost 几何平均减少 23.6%，最差一组增加约 0.30%。
这些是图通信/存储指标，不能换算为 PDE 加速比。

策略同时约束精确 owned 配额和面关联工作量，使用割边代价优先队列与最多展开 64 节点的
局部连通证明完成配额修正。预算只限制一次性证明开销，不改变 halo 宽度。
同一 512,000 单元 Voronoi 问题中，修正从 23.694 秒降到 1.491 秒，所有权哈希不变。
完整修正过的 Voronoi 工程基底由真实面环和框架几何独立验证，旧失真工程图已排除。

实际 500,094 单元案例的 C++ CSR、权重、初始/最终所有权与模型逐单元对齐。
工程规模审计发现快速浮点倒数变换曾改变一次性 METIS 比例的末位；最终版本仅对归一化使用普通除法，
求解循环编译优化保持原设置。修正前二进制 `44f8bb…` 的记录保留，未混入最终强扩展。
模型方法、选择依据、逐图结果和官方参考资料见 [划分研究报告](partition-model-2026-10-07.md)。
实际局部面负载、halo 与邻接 rank 指标见 [最终 CSV]({NAME}/final-scaling.csv)。

## 数值、兼容性与边界检查

最终版本 `make -j4 test-full` 全部通过；编译可并行，测试组和组内求解顺序运行。
包括 78 项源码/头文件分层检查、标量/向量后端、矩阵局部变化的集体更新、RHS 重复求解、
复用增长与失败恢复、warm start、typed 配置、框架外 SDK 链接、时间推进、RANS 阶数、
收敛 SIMPLE 全字段一致性，以及 MPI 网格/数学/非结构网格/通道/腔流/输运/热传导/Poiseuille。
root 读取与工作 rank 无效网格路径、root 错误的集体传播、不连通凹多边形棱柱及完整面环均有检查。
不启用 METIS 的划分模块单独编译并通过 4 rank 回退测试；没有将它称为完整无 METIS 构建矩阵。
[完整日志]({NAME}/tests/test-full-normalized-partition.log)保留原断言和容差。

工程 PISO 最终输出按全局 ID 比较 U、p，要求单元覆盖唯一且完全相同、字段有限，
采用现有 MPI 一致性检查的 `atol=rtol=5e-6`。报告实际最大差，不以该检查改变求解容差。

{field_table}
SIMPLE 吞吐达到外迭代上限时正式求解器按原语义不输出最终场，因此没有宣称进行了该吞吐状态的
工程字段对照；收敛 SIMPLE 的字段与守恒验证来自完整回归。PISO 的六步实验也不是长期物理模型验证。

阶段计数覆盖完整求解回调，包含首次迭代，与纯计算窗口不是相同范围。
显式 halo/归约计数不包括 PETSc/Hypre 内部通信，SpMV/PC 应用计数尚不覆盖它们的内部调用。
所以不能用这些计数为零来判断线性求解没有 SpMV，也不能由割面下降断言 AMG 扩展性必然改善。
压力求解与预条件器的阶段时间及进程数曲线用于界定下一处成本，AMG 粗层与网络/带宽贡献尚未单独分离。

## 可复现证据

最终二进制 SHA256：`{summary['finalBinarySha256']}`。
内存优化前参考 SHA256：`{summary['referenceBinarySha256']}`。
输入、绑定、构建参数、完整计数、stdout、RSS/PSS/交换页 trace 见
[同版强扩展]({NAME}/final-strong-scaling/)、[12 进程重测]({NAME}/final-repeat-12/)。
源码、原 HEAD、构建开关与 Physics 未改动检查保存在 [源码清单]({NAME}/final-source/manifest.json)；
tracked 补丁与新增源码分别保存，临时二进制路径记录在各 study.json 中。

```bash
make -j4 all
make -j4 test-full
python3 tools/benchmark_solver_optimization.py \\
  --variant final=/absolute/path/babelsim-memory-partition-normalized \\
  --case simple=/absolute/path/simple128 --throughput simple=128 \\
  --case piso=/absolute/path/piso6 --accepted-steps piso=6 \\
  --ranks 1,4,8,12 --repeat 1 --warmup 0 \\
  --output /tmp/new-evidence --work-root /tmp/new-runs
python3 docs/reports/{NAME}/analyze-results.py
python3 docs/reports/{NAME}/render-report.py
```

正式输入小文件已逐项保存，网格 SHA256 记录在各 study.json 的 manifest 中。
原始生成器及参数 `gridN=63, warpFraction=0.15` 可重建网格；输入快照保留正式吞吐/时间推进配置。
CSV 与 PNG/SVG/PDF 图表由分析脚本生成。早期 `clean-pilot` 的 SIMPLE 48 次外迭代曾由 130.867 秒降到
66.013 秒；那是复用策略筛选记录，不是最终 128 次外迭代的速度比。
没有收益的 transpose-cache 候选（PISO 148.691/148.738 秒）未保留到当前实现。
'''
    destination = ROOT.parent/'flow-optimization-2026-10-07.md'
    destination.write_text(text)
    print(destination)


if __name__ == '__main__':
    main()
