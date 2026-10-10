# 算例字典参考

[配置参考](README.md) · [算例管理](../user-guide/case-management.md)

## 通用语法

每行一个 `key value...`，空白分词，`#` 开始注释。普通参数字典不使用花括号、引号、分号或表达式。字段名与枚举大小写敏感。必需值缺失、重复项、未知项、错误数字和未被求解器消费的配置会报错。

数字需满足各项有限性与范围要求。参数接口的布尔值允许 `true/yes/on/1` 和 `false/no/off/0`，布尔词不区分大小写；具体读取为整数的参数仍必须使用其指定整数形式。

## case.bs

| 键 | 必需/默认 | 含义与约束 |
| --- | --- | --- |
| `solver` | 必需 | 已链接的求解器注册名 |
| `mesh` | 必需 | 原生 mesh 3 文件路径 |
| `fields` | 必需 | 直接存放 `<name>.field` 的目录 |
| `physics` | 必需 | 物理参数字典路径 |
| `methods` | 必需 | 数值离散字典路径 |
| `solution` | 必需 | 算法与线性配置字典路径 |
| `control` | 必需 | 时间字典路径 |
| `output` | 必需 | 写出字典路径，可采用默认内容 |
| `ghostLayers` | `3` | 幽灵单元层数，整数且至少 3 |

输入路径按算例根目录解析，可以指定绝对输入路径；交付推荐使用相对路径以便搬迁。路径格式不支持带空格的引号字符串。`fields` 不会自动追加 `/initial`。

```text
solver piso
mesh mesh/domain.mesh
fields fields/initial
physics physics/model.bs
methods numerics/methods.bs
solution numerics/solution.bs
control control.bs
output output.bs
ghostLayers 3
```

## control.bs

| 键 | 必需/默认 | 约束 |
| --- | --- | --- |
| `startTime` | 必需 | 有限数值 |
| `endTime` | 必需 | 有限，且不小于 startTime |
| `deltaT` | 必需 | 有限正数，步数须在可表示范围内 |

Euler 时间驱动可截短最后一步。BDF2 的当前配置驱动要求区间是固定步长的整数倍，否则报 `BDF2 requires an integral number of uniform time steps`。公开低层时间算子能够使用上一时间步长度，这不改变当前配置驱动的限制。

`time steady` 在 methods 中选择，不写在 control。稳态 SIMPLE 不通过物理步循环推进；热传导/输运的 steady 装配仍由非空 control 时间区间触发。

## output.bs

| 键 | 默认 | 含义 |
| --- | --- | --- |
| `directory` | `results` | 相对算例根目录的输出位置 |
| `timeName` | `final` | 单个非数字运行名；命令行 -time 可覆盖 |
| `writeInterval` | `1` | 按物理时间步编号写出间隔，正整数 |
| `writeFields` | 空 | 空时默认选择从文件加载的单元场；否则限定字段列表 |
| `excludeFields` | 空 | 取消所列场的默认输出选择 |

接受兼容别名 `time_name`、`fields`、`exclude`，分别对应 `timeName`、`writeFields`、`excludeFields`。新文档统一使用表中的正式写法，同一参数不要通过两个名称重复定义。

输出目录必须为相对路径且不包含 `.` 或 `..` 组件。`timeName` 与求解 `-time` 的运行名必须为单组件，不能是纯数值、`all` 或 `latest`。空文件可使用输出默认值，case 的 `output` 文件路径本身仍为必需。

```text
directory results
timeName baseline
writeInterval 5
writeFields U p mut
excludeFields wallDistance
```

这里 `mut` 必须确实由启用的模型创建；层流算例中不能保留这个不存在的输出场。程序创建的单元场可由 `Case::output(field)` 显式选择；面场当前不支持写入此结果格式。

## physics、methods、solution

物理键依赖实际 solver 和 turbulenceModel，见[物理参数](physical-models.md)。methods 使用完整的逐方程配置，solution 使用算法控制和逐方程线性配置，见[数值配置](numerical-controls.md)。没有“所有求解器都生效”的物理或算法参数全集文件；新增 C++ 求解器自行声明并消费相应参数。

## 错误定位

错误通常包含文件、键或配置名称。先检查 `case.bs` 路径，再检查该文件里的准确拼写与值。多个实例写入同一结果位置、旧字典项未删除或新场未配置全部 patch，均会破坏复现或触发校验。
