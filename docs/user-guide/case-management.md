# 04 算例管理

[用户手册](README.md) · 上一章：[首次运行](quick-start.md) · 下一章：[网格与边界条件](mesh-and-boundaries.md)

## 推荐目录

```text
myCase/
  case.bs
  control.bs
  output.bs
  mesh/domain.mesh
  fields/initial/T.field       或 U.field、p.field 等
  physics/model.bs
  numerics/methods.bs
  numerics/solution.bs
  results/                    计算产生
  post/                       后处理产生
  logs/                       使用者保存的启动与诊断日志
```

软件按 `case.bs` 中的路径读取输入，子目录名称可以更改。`fields` 条目必须直接指向存放 `<场名>.field` 的目录；不能写 `fields fields` 后期待软件自动追加 `initial`。所有相对路径均以算例根目录为基准，只有 `internal file` 的相对路径以对应 `.field` 文件所在目录为基准。

## 修改文件的顺序

| 需要改变 | 应修改 |
| --- | --- |
| 求解器、网格或输入文件位置 | `case.bs` |
| 密度、导热率、黏度、模型及源项 | `physics` 指向的字典 |
| 初始值、入口值、壁面与出口条件 | `fields` 指向的场文件 |
| 时间离散、插值、梯度、对流、扩散 | `methods` 指向的字典 |
| 外层耦合、松弛、线性容差与预条件器 | `solution` 指向的字典 |
| 起始时间、结束时间、时间步 | `control` 指向的字典 |
| 保存目录、运行名称、场列表与间隔 | `output` 指向的字典 |

更换求解器时，应同步修改字段、物性和方程配置，删除旧求解器专用且不再使用的键。严格配置校验会拒绝未使用或拼错的项，防止设置看似存在却没有影响计算。

## 文本格式

普通字典每行一个键及其值，以空白分隔，`#` 开始注释。键和值大小写敏感，重复键、未知键、缺少必需项和错误数值会报错。不要添加 JSON 的引号、逗号或其他软件字典的分号。

```text
# 时间设置
startTime 0
endTime 1
deltaT 0.01
```

场文件使用花括号组织 `boundary` 和 patch 条件。参数文件并不是 C++ 或 Python 脚本，不支持表达式求值、单位自动换算或通用 `include` 指令。路径含空格的输入不适合当前字典语法，交付时使用简单路径名称。

## 可复现的运行记录

给每组物理与数值设置建立独立输入目录。保存源码版本或交付包校验值、构建选项、MPI 进程数、运行命令和退出码。运行名推荐使用 `baseline`、`dt001`、`meshFine` 等单个非数字目录名；`all`、`latest` 和纯数字名保留给后处理与物理时刻。

同一输入可以用不同运行名比较 MPI：

```bash
build/babelsim-solve -case myCase -time np1
mpirun -np 4 build/babelsim-solve -case myCase -time np4
```

同名运行写入同一结果位置，不能作为历史版本留存机制。修改时间范围或写出间隔后重用旧结果目录，还可能保留先前的时刻文件；正式比较使用新运行名和完整输入快照。

## 不同类型的目录

`cases/` 中提供产品示例与研究算例；各算例 README 给出具体生成器、规模和已有证据。`docs/examples` 提供手册的小型可执行示例。研究报告与历史开发记录分别保存在 `archive/docs/reports`、`archive/docs/history`，已移出交付 docs。历史报告中的旧参数或构建命令只说明当时的实验，不应代替当前配置参考。
