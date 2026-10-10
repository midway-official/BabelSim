# 命令行参考

[配置参考](README.md) · [运行手册](../user-guide/running-and-parallel.md)

## babelsim-solve

```text
build/babelsim-solve -case <case-directory> [-time <run-name>] [-performance <directory>]
```

| 参数 | 必需 | 含义 |
| --- | --- | --- |
| `-case` | 是 | 包含 `case.bs` 的目录 |
| `-time` | 否 | 覆盖 output 的结果运行名；不是起始物理时刻或重启选择 |
| `-performance` | 否 | 写出每个进程的性能 JSON 到指定目录 |

解析器拒绝未知参数和缺少参数值的调用。求解器注册名从 `case.bs` 读取；正式程序包含 `heat`、`transport`、`simple`、`transientSimple`、`piso`。注册名称区分大小写。没有 `--help`、`--version`、`--list-solvers`、`--dry-run` 或命令行字典覆盖选项。

外部示例程序使用相同的 `runApplication` 参数格式，但只包含链接进该程序的注册项。不能期待只链接 `libbabelsim.a` 的自定义程序自动拥有全部正式求解器注册项。

应用成功返回 0；配置/启动错误返回 1；求解未收敛或数值失败返回 2。MPI 错误可能由启动器额外包装。生产算法必须显式报告状态与写出数据，入口不会替算法自动写出成功结果。

## babelsim-post

```text
build/babelsim-post -case <case-directory> [-time <selection>] -format vtk [tecplot]
```

`-format` 必须提供至少一个当前支持的格式：`vtk`、`tecplot`；可以只选其中任意一个，也可同时选择。后处理本身使用单进程启动，不必再次用 mpirun 分发原始 rank 数据。

| selection | 选择规则 | 典型产物 |
| --- | --- | --- |
| `runA` | `results/runA` 的快照 | `post/runA.vtu` |
| `runA/0.05` | 指定命名运行时刻 | `post/runA/0.05.vtu` |
| `runA/latest` | 命名运行最后数值时刻 | `post/runA/<实际时刻>.vtu` |
| `runA/all` | 命名运行全部数值时刻 | `post/runA/series.pvd` 与逐时刻 VTU |
| `latest`、`all` | 结果根目录未命名序列 | `post/<时刻>.vtu` 或 `post/series.pvd` |
| 省略 | 使用 output 的 `timeName` | 对应快照文件 |

Tecplot 对应扩展名 `.dat`。VTK 序列才生成 PVD；单快照不创建序列索引。路径组件不能借助 `..` 等形式越出结果目录。后处理失败返回非零，检查错误日志中的缺失字段、rank、元数据和网格差异。

## 手册示例生成器

```bash
python3 docs/examples/create_case.py --template heat --output work/heat
```

`--template` 必需，可选 `heat`、`transport`、`simple`、`piso`、`custom-heat`、`coupled`；`--output` 必需，指定新的算例根目录。存在的目录或符号链接不会被覆盖。前四种使用正式求解器，后两种使用[示例程序](../examples/README.md)。

```bash
make -f docs/examples/Makefile -j4
python3 docs/examples/verify_examples.py --output build/docs-refresh
```

示例 Makefile 必须从项目根目录执行，复用根 Makefile 的库与链接依赖。验证工具使用默认 `build/` 中的程序，创建独立目录并保留日志、输入和报告；不适用于用其他 `BUILD` 的程序验证。

## 常用维护命令

```bash
./dependencies/install.sh --help
python3 tools/check_docs.py --json build/docs-refresh/docs-check.json
python3 tools/check_docs.py --delivery-only --json build/docs-refresh/delivery-check.json
python3 tools/compare_parallel_results.py --help
```

第二个文档检查入口只检查产品手册、根 README 和依赖说明，不要求 Git 或可选 archive。安装参数见[安装手册](../user-guide/installation.md)。不要把 Python 工具或安装脚本的 `--help` 推广为正式求解程序也支持该选项。
