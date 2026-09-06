#!/usr/bin/env python3
"""汇总 BabelSim 日志中的完成时间、固定外迭代吞吐和后端计数。"""

import argparse
import re
from pathlib import Path


PERFORMANCE = re.compile(r"BabelSim performance (?P<body>.*)")
TIME_TO_SOLUTION = re.compile(r"BabelSim timeToSolution=([^\s]+)")
RESULT = re.compile(r"BabelSim result .*?\bsteps=(\d+)\b")
ITERATION = re.compile(r"(?:Transient )?SIMPLE\s+(\d+)\b")


def read_log(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")
    performance = None
    for match in PERFORMANCE.finditer(text):
        values = {}
        for token in match.group("body").split():
            if "=" not in token:
                continue
            name, value = token.split("=", 1)
            if "/" in value:
                count, seconds = value.split("/", 1)
                values[name + "Count"] = float(count)
                values[name + "Seconds"] = float(seconds)
            else:
                values[name] = float(value)
        performance = values
    if performance is None:
        raise ValueError(f"{path}: 未找到 BabelSim performance 记录")
    wall = [float(value) for value in TIME_TO_SOLUTION.findall(text)]
    if not wall:
        raise ValueError(f"{path}: 未找到 BabelSim timeToSolution 记录")
    performance["timeToSolution"] = wall[-1]
    iterations = [int(value) for value in ITERATION.findall(text)]
    steps = [int(value) for value in RESULT.findall(text)]
    outer = iterations[-1] if iterations else (steps[-1] if steps else 0)
    return performance, outer


def main():
    parser = argparse.ArgumentParser(
        description="同时报告 time-to-solution 与固定外迭代的单位工作量性能")
    parser.add_argument("logs", nargs="+", type=Path)
    parser.add_argument(
        "--require-outer", type=int,
        help="要求每个日志恰好完成该外迭代数，用于固定工作量比较")
    arguments = parser.parse_args()

    print("| 日志 | time-to-solution/s | 运行域/s | 外迭代 | s/外迭代 | Krylov/外迭代 | SpMV | halo | Allreduce | 装配/s | 预条件构建/s | 预条件应用/s | 线性求解/s |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for path in arguments.logs:
        values, outer = read_log(path)
        if arguments.require_outer is not None and outer != arguments.require_outer:
            raise SystemExit(
                f"{path}: 外迭代为 {outer}，期望 {arguments.require_outer}；"
                "该日志不能用于固定工作量比较")
        elapsed = values["runtimeElapsed"]
        time_to_solution = values["timeToSolution"]
        divisor = max(outer, 1)
        print(
            f"| {path.name} | {time_to_solution:.6g} | {elapsed:.6g} | {outer} | "
            f"{elapsed/divisor:.6g} | "
            f"{values.get('krylovIterations', 0)/divisor:.6g} | "
            f"{int(values.get('spmv', 0))} | {int(values.get('halo', 0))} | "
            f"{int(values.get('allreduce', 0))} | "
            f"{values.get('assemblySeconds', 0):.6g} | "
            f"{values.get('preconditionerSetupSeconds', 0):.6g} | "
            f"{values.get('preconditionerApplySeconds', 0):.6g} | "
            f"{values.get('linearSeconds', 0):.6g} |")


if __name__ == "__main__":
    main()
