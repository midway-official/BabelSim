# 05 网格与边界条件

[用户手册](README.md) · 上一章：[算例管理](case-management.md) · 下一章：[选择物理模型](models.md)

## 网格输入

正式读取器使用 `BABELSIM_MESH 3`：文件包含顶点坐标、每个面的顶点连接、owner/neighbour 和边界 patch。几何为三维；薄层二维算例仍需要正的厚度与体积。运行时不直接读取其他网格软件的文件，转换程序必须生成合法的原生格式。

生成一个 32 × 8 × 1 的直通道网格：

```bash
python3 tools/generate_hexa_mesh.py work/channel.mesh \
  --cells 32 8 1 --bounds 0 0 0 4 1 0.1 \
  --patch minus_x inlet inlet --patch plus_x outlet outlet \
  --patch minus_y lower_wall wall --patch plus_y upper_wall wall \
  --patch minus_z front symmetry --patch plus_z back symmetry
```

`--bounds` 依次为 `x0 y0 z0 x1 y1 z1`；三个方向的上界均须大于下界。`--patch` 依次指定几何侧面、patch 名和 patch 类型。可用 `--tanh-clustering` 调整生成网格的坐标聚集程度，但网格质量仍需独立检查。

完整格式、编号与面朝向规则见[文件格式参考](../reference/file-formats.md)。手册提供的[单单元网格](../examples/one-cell.mesh)可直接由正式读取器使用。

## 网格验收

至少核对：所有单元体积为正；内部面只有一对 owner/neighbour；边界面属于有效 patch；面顶点朝向与 owner 的外法向一致；封闭性、尺寸、非正交程度和边界位置符合计算要求。读取器会拒绝多类非法输入，但成功读入不等于网格适合工程求解。

所有物理 patch 都必须在每个输入场的 `boundary` 中出现一次。MPI 分区产生的 processor 边界由框架管理，不需要写入用户的场文件。

## 初始场和边界文件

热传导完整标量场示例：

```text
field T
{
    type scalar
    location cell
    internal uniform (0)
    boundary
    {
        hot   { type fixedValue value (1) }
        cold  { type fixedValue value (0) }
        lower { type symmetry }
        upper { type symmetry }
        front { type symmetry }
        back  { type symmetry }
    }
}
```

`field T`、文件名 `T.field` 及求解器加载的 `T` 必须一致。`scalar` 值有 1 个分量；`vector` 有 3 个，如 `(1 0 0)`；`tensor` 有 9 个，按行优先排列。文件加载的场目前使用 `location cell`；面场由 C++ 接口显式创建。

## 条件的物理含义

| 类型 | 配置 | 含义与用途 |
| --- | --- | --- |
| `fixedValue` | `value (...)` | 指定边界场值；定温、指定入口速度、无滑移速度等 |
| `fixedGradient` | `value (...)` | 指定沿外法向的导数；不是已经乘扩散系数的热流/通量 |
| `zeroGradient` | 无 `value` | 外法向梯度为零，常用于出口或绝热标量边界 |
| `symmetry` | 无 `value` | 标量零法向梯度，向量使用法向投影约束 |
| `inletOutlet` | `value (...)` | 面通量流入时采用指定值，流出时零梯度 |

对于导热，外向热流密度为 `q_out = -conductivity * dT/dn`。若要给定热流，应先将热流换成 `fixedGradient` 的梯度值并核对符号。`inletOutlet` 需要求解或运算具有面通量上下文，不能当作与流向无关的普通定值边界。

网格的 `wall` 类型只标记几何边界，不自动设置无滑移；速度场中仍需写 `type fixedValue value (0 0 0)`。压力出口常采用 `fixedValue (0)`，入口与无滑移壁面常采用 `zeroGradient`。实际边界组合应与物理问题一致。

## 非均匀初值

将 `internal uniform (...)` 改为：

```text
internal file ../../initial-values/T.csv
```

路径相对于当前 `.field` 文件；数据使用全局单元编号，必须覆盖网格所有单元且不能重复。也可指向已有结果时刻目录，读取其中的 `rank-*/T.csv`。每个场仍保留自己的边界设置。数据格式及继续计算步骤见[结果处理与继续计算](results-and-restart.md)。

## 需要专门核对的边界

当前向量 `symmetry` 的隐式扩散使用共享标量系数与边界投影，切向响应需要按具体方向和网格独立验证；不能仅凭 `symmetry` 名称认定已完整验证任意方向的自由滑移。标量绝热条件和薄层标量教程不涉及这一向量问题。

纯 Neumann 扩散问题还需满足体源与边界通量相容。参考压力/标量值只是消除常数不确定性的一部分，不能修复不相容的输入。当前字典没有提供 Robin、周期边界或任意时间函数边界；如需这些行为，应开发并验证相应接口。
