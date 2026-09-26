# MPI 审查证据

对应 [审查报告](../mpi-audit-2026-09-26.md)。基线 commit 为 `deecd86`，工作树改动见报告。

- `halo-before.log`：新不对称更新回归在原后端上失败，mpirun 返回 1。
- `mpi-tests.log`、`unit-tests.log`：最终 MPI / 串行测试输出。
- `backend-extra.log`：PETSc 1/3/8 进程补充验证。
- `partition.log`：2/3/4/8 进程，三种网格的单元覆盖、均衡、halo 正确性及规模指标。
- `consistency-3layers.json`、`consistency-extra.json`：共 45 次流场运行和 110 次状态比较，包括二进制 SHA256、基线 commit、容差与误差。JSON 中 `/tmp` 为本次运行原始输出路径，不是永久交付位置。
- `partition_benchmark.cpp`：基线与优化函数的直接 owners 对照，15 组；`partition-benchmark.csv` 每列为网格边长、分区数、基线秒、新版本秒、归属一致性。单次采样，无统计意义上的加速保证。

复现：

```sh
make -j4 test
make -j4 test-mpi
python3 tests/simple_parallel_consistency_test.py --solver build-petsc/babelsim-solve
TMPDIR=/tmp mpirun --oversubscribe -np 8 build-petsc/parallel_partition_test
TMPDIR=/tmp mpirun --oversubscribe -np 8 build-petsc/petsc_backend_test
mpic++ -O3 -std=c++17 -DOMPI_SKIP_MPICXX=1 -Iinclude -Isrc -Itests \
  -I/usr/include/eigen3 docs/reports/mpi-audit-2026-09-26-evidence/partition_benchmark.cpp \
  build-petsc/libbabelsim.a -o /tmp/babelsim-partition-benchmark
/tmp/babelsim-partition-benchmark
```

命令从仓库根目录执行。oversubscribe 仅用于正确性覆盖，不用于性能评估。
