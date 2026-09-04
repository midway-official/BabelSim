#!/usr/bin/env python3
"""将 BabelSim 方腔中心线速度与 Ghia 等（1982）的数据比较。"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


GHIA_Y = np.array(
    [1.0, .9766, .9688, .9609, .9531, .8516, .7344, .6172,
     .5, .4531, .2813, .1719, .1016, .0703, .0625, .0547, 0.0]
)
GHIA_U = {
    100: np.array(
        [1.0, .84123, .78871, .73722, .68717, .23151, .00332, -.13641,
         -.20581, -.21090, -.15662, -.10150, -.06434, -.04775, -.04192,
         -.03717, 0.0]),
    400: np.array(
        [1.0, .75837, .68439, .61756, .55892, .29093, .16256, .02135,
         -.11477, -.17119, -.32726, -.24299, -.14612, -.10338, -.09266,
         -.08186, 0.0]),
    1000: np.array(
        [1.0, .65928, .57492, .51117, .46604, .33304, .18719, .05702,
         -.06080, -.10648, -.27805, -.38289, -.29730, -.22220, -.20196,
         -.18109, 0.0]),
}
GHIA_X = np.array(
    [1.0, .9688, .9609, .9531, .9453, .9063, .8594, .8047,
     .5, .2344, .2266, .1563, .0938, .0781, .0703, .0625, 0.0]
)
GHIA_V = {
    100: np.array(
        [0.0, -.05906, -.07391, -.08864, -.10313, -.16914, -.22445,
         -.24533, .05454, .17527, .17507, .16077, .12317, .10890, .10091,
         .09233, 0.0]),
    400: np.array(
        [0.0, -.12146, -.15663, -.19254, -.22847, -.23827, -.44993,
         -.38598, .05186, .30174, .30203, .28124, .22965, .20920, .19713,
         .18360, 0.0]),
    1000: np.array(
        [0.0, -.21388, -.27669, -.33714, -.39188, -.51550, -.42665,
         -.31966, .02526, .32235, .33075, .37095, .32627, .30353, .29012,
         .27485, 0.0]),
}


def _read_rows(source: Path) -> np.ndarray:
    """读取串行 CSV，或按 global_id 合并并行 rank 目录。"""
    files = [source]
    if source.is_dir():
        files = sorted(source.glob("rank-*/U.csv"))
        if not files:
            raise ValueError(f"{source} 中没有 rank-*/U.csv")
    rows = [np.atleast_1d(np.genfromtxt(path, delimiter=",", names=True))
            for path in files]
    names = rows[0].dtype.names
    if names is None or any(row.dtype.names != names for row in rows):
        raise ValueError(f"{source} 中的 CSV 表头不一致")
    data = np.concatenate(rows)
    if source.is_dir():
        if "global_id" not in names:
            raise ValueError(f"{source} 的并行结果缺少 global_id")
        identifiers = data["global_id"].astype(np.int64)
        if np.unique(identifiers).size != identifiers.size:
            raise ValueError(f"{source} 的并行结果包含重复 global_id")
        data = data[np.argsort(identifiers)]
    return data


def _grid(source: Path) -> tuple[np.ndarray, ...]:
    data = _read_rows(source)
    required = {"x", "y"}
    if data.dtype.names is None or not required.issubset(data.dtype.names):
        raise ValueError(f"{source} 必须包含列 {sorted(required)}")
    # 通用输出使用 value0/value1；旧的单文件回归格式使用 u/v，二者均可读取。
    if {"u", "v"}.issubset(data.dtype.names):
        u, v = data["u"], data["v"]
    elif {"value0", "value1"}.issubset(data.dtype.names):
        u, v = data["value0"], data["value1"]
    else:
        raise ValueError(f"{source} 必须包含 u/v 或 value0/value1")
    if not all(np.isfinite(data[name]).all() for name in required) or \
       not np.isfinite(u).all() or not np.isfinite(v).all():
        raise ValueError("方腔 CSV 包含非有限值")
    # 多面体质心可能相差几个 ulp，拓扑判断不能依赖坐标逐位相等。
    rounded_x = np.round(data["x"], 12)
    rounded_y = np.round(data["y"], 12)
    xs = np.unique(rounded_x)
    ys = np.unique(rounded_y)
    if xs.size * ys.size != data.size:
        raise ValueError("方腔 CSV 不是单个结构化 z 层")
    order = np.lexsort((rounded_x, rounded_y))
    shape = (ys.size, xs.size)
    return (
        data["x"][order].reshape(shape),
        data["y"][order].reshape(shape),
        u[order].reshape(shape),
        v[order].reshape(shape),
    )


def centreline(source: Path) -> tuple[np.ndarray, ...]:
    """返回含壁面端点的 y/u 与 x/v 中心线。"""
    x, y, u, v = _grid(source)
    u_line = np.array([np.interp(0.5, x[row], u[row])
                       for row in range(y.shape[0])])
    v_line = np.array([np.interp(0.5, y[:, column], v[:, column])
                       for column in range(x.shape[1])])
    return (
        np.r_[0.0, y[:, 0], 1.0], np.r_[0.0, u_line, 1.0],
        np.r_[0.0, x[0], 1.0], np.r_[0.0, v_line, 0.0],
    )


def ghia_errors(source: Path, reynolds: int) -> tuple[np.ndarray, ...]:
    """在 Ghia 表格坐标处计算误差；Re=400 排除已知异常点。"""
    y, u, x, v = centreline(source)
    calculated_u = np.interp(GHIA_Y, y, u)
    calculated_v = np.interp(GHIA_X, x, v)
    u_error = calculated_u - GHIA_U[reynolds]
    v_error = calculated_v - GHIA_V[reynolds]
    v_mask = np.ones(v_error.shape, dtype=bool)
    if reynolds == 400:
        v_mask[np.flatnonzero(np.isclose(GHIA_X, .9063))[0]] = False
    return calculated_u, calculated_v, u_error, v_error, v_mask


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path,
                        help="串行 U.csv，或包含 rank-*/U.csv 的并行结果目录")
    parser.add_argument("--re", type=int, choices=sorted(GHIA_U), default=100)
    parser.add_argument("--u-l2-limit", type=float, default=0.01)
    parser.add_argument("--v-l2-limit", type=float, default=0.01)
    args = parser.parse_args()
    if GHIA_Y.shape != GHIA_U[args.re].shape or \
       GHIA_X.shape != GHIA_V[args.re].shape:
        raise RuntimeError("Ghia 坐标与参考数组长度不一致")
    calculated_u, calculated_v, u_error, v_error, v_mask = ghia_errors(
        args.csv, args.re)
    # 原论文 Re=400、x=0.9063 的 v=-0.23827 与相邻点和论文曲线不连续，后续文献
    # 普遍将其视为排印错误。保留原值供追溯，但不让这一点污染误差范数。
    u_l2 = float(np.sqrt(np.mean(u_error * u_error)))
    v_l2 = float(np.sqrt(np.mean(v_error[v_mask] * v_error[v_mask])))
    x, y, _, _ = _grid(args.csv)
    print(
        f"Re={args.re} grid={x.shape[1]}x{y.shape[0]} "
        f"u_center={calculated_u[8]:.8f} v_center={calculated_v[8]:.8f} "
        f"u_Linf={np.max(np.abs(u_error)):.8g} u_L2={u_l2:.8g} "
        f"v_Linf={np.max(np.abs(v_error[v_mask])):.8g} v_L2={v_l2:.8g}"
    )
    if u_l2 > args.u_l2_limit or v_l2 > args.v_l2_limit:
        raise SystemExit("Ghia 中心线误差超过配置门限")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
