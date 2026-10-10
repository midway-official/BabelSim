# 12 交付、维护与验收

[用户手册](README.md) · 上一章：[故障排查](troubleshooting.md) · [开发验证指南](../developer-guide/verification.md)

## 交付内容

BabelSim 当前采用源码交付与本地构建方式。交付包应当使接收方在准备系统依赖后，从同一目录完成 `make`、运行教程、查询全部配置并编写外部求解器。

| 内容 | 交付作用 |
| --- | --- |
| README、Makefile、LICENSE | 产品入口、构建入口和项目许可 |
| include、src | 公开开发接口与实现 |
| dependencies | 精简 PETSc/Eigen/Hypre 源码、版本/校验清单和安装脚本 |
| docs/user-guide | 完整使用手册 |
| docs/reference | 精确参数与文件格式 |
| docs/developer-guide、docs/examples | 二次开发、验证和可运行代码 |
| tools、tests、必要 cases | 网格生成、比较工具、回归与产品算例 |
| docs/developer-guide 的性能章节 | 高级测量方法 |
| 验收清单与报告 | 本次交付版本、环境、实际命令及结果 |

`build/` 是可重建产物；通常不需要作为源码依赖交付。用户计算结果、长期研究日志、历史冻结源码和开发归档可单独提供，产品操作文档不以这些大型数据为前提。

## 版本与来源

记录本次源码提交或交付包 SHA-256，保留 dependencies/sources.json 与 pruning.json。当前依赖清单固定 PETSc 3.25.5、Eigen 3.4.0 和所需 Hypre 修复版本。文档参考当前源码接口，日期用于记录核对时间，不虚构发布版本或服务承诺。

```bash
git rev-parse HEAD
git status --short
sha256sum dependencies/distfiles/hypre-*.tar.gz
```

对于带本地修改的工作目录，提交号不能完整描述交付内容，应同时保留修改文件或最终打包校验值。离线源码构建仍需要预先安装的编译器、MPI 和 BLAS/LAPACK；缺失系统包的安装可能需要网络和系统权限。

## 重建与首次验收

在新解压位置执行：

```bash
./dependencies/install.sh --check
make -j4 PETSC_JOBS=4
make petsc-check
python3 tools/check_docs.py --delivery-only --json build/docs-check.json
python3 docs/examples/verify_examples.py --output build/tutorial-validation
```

检查依赖来源为本包，库生成在当前 build/petsc，正式应用与示例都能运行。交付检查不要求 .git 或可选 archive 目录，适用于解压后的源码包。按报告核对串行/MPI 数据、导出文件及错误路径。安装脚本另可直接 `./dependencies/install.sh` 检查并安装缺少系统依赖后构建。

需要验收完整框架行为时，按[测试矩阵](../developer-guide/verification.md)顺序运行 test-full 和额外 PISO 耦合组。教程通过、完整回归通过与物理基准通过分别记录，不合并为单个“全部正确”的声明。

## 搬迁与运行库

本地构建包含对应 PETSc lib 的动态搜索路径，复制已编译程序到另一绝对目录后不一定能直接运行。源码交付推荐在新位置重新构建；如交付二进制，应核对 `ldd`、MPI/BLAS ABI、CPU 指令集和运行时路径，并在接收环境实际运行。

默认本机优化不保证跨 CPU；可移植构建需要同时修改框架与 PETSc 的编译选项，见[安装手册](installation.md)。精简依赖只支持当前选择的配置，扩展 GPU、外部包或 Eigen 模块应更新来源和验收范围。

## 升级与维护

升级前保存独立输入、结果及旧版本证据；新版本使用新构建目录或清理对应旧产物。检查配置名、场名、数值格式与公开 API 是否变化，先运行小型教程，再比较相同输入和时刻，最后重复项目所需物理验证。

旧平铺文档已迁移到用户/参考/开发三类目录，迁移表见[文档维护规范](../developer-guide/documentation.md)。archive/docs 中的旧 API、命令和性能数据是追溯资料，不用于当前安装和操作，也不是交付 docs 的组成部分。

## 清理与许可

`make clean` 删除当前构建目录，不删除用户算例结果或卸载 apt 包。卸载本地源码部署时先保存需要的算例与数据，再移除对应源码/构建目录；没有单独 make uninstall 目标。

BabelSim 项目采用 [MIT License](../../LICENSE)。PETSc、Eigen、Hypre 等依赖保留各自许可和版权文件，见[依赖说明](../../dependencies/README.md)。再次分发时应保留源码来源、依赖清单和相应许可文件；商业交付说明应依据实际交付范围，不添加未验证的性能指标或物理认证声明。
