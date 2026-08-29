#!/usr/bin/env python3
"""Merge BabelSim rank-owned field CSV files into VTK and Tecplot files.

The rank writer stores cell-centre samples, not cell connectivity.  Therefore
the VTK output is a point cloud (one vertex per owned cell), which is directly
readable by ParaView and preserves global ids without duplicating halo cells.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


FIELDS = ("x", "y", "z", "u", "v", "w", "p")


def _rank_files(directory: Path, prefix: str, ranks: int | None) -> list[Path]:
    if not directory.is_dir():
        raise ValueError(f"output directory does not exist: {directory}")
    if not prefix or Path(prefix).name != prefix:
        raise ValueError("prefix must be a non-empty file prefix")

    exact = sorted(directory.glob(f"{prefix}_*.csv"))
    if exact:
        files = exact
    elif ranks is not None:
        files = sorted(directory.glob(f"{prefix}{ranks}_*.csv"))
    else:
        candidates = sorted(directory.glob(f"{prefix}[0-9]*_*.csv"))
        groups: dict[str, list[Path]] = {}
        for path in candidates:
            match = re.fullmatch(
                rf"{re.escape(prefix)}([0-9]+)_([0-9]+)\.csv", path.name)
            if match:
                groups.setdefault(match.group(1), []).append(path)
        if len(groups) != 1:
            raise ValueError(
                f"no unique rank-file set for prefix {prefix!r}; pass --ranks")
        files = next(iter(groups.values()))

    if ranks is not None:
        expected = {
            directory / f"{prefix}_{rank}.csv" for rank in range(ranks)}
        if not exact:
            expected = {
                directory / f"{prefix}{ranks}_{rank}.csv"
                for rank in range(ranks)}
        if set(files) != expected:
            missing = sorted(str(path) for path in expected - set(files))
            extra = sorted(str(path) for path in set(files) - expected)
            details = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if extra:
                details.append("unexpected " + ", ".join(extra))
            raise ValueError("; ".join(details))
    if not files:
        raise ValueError(f"no rank CSV files found in {directory}")
    return files


def _load(files: list[Path]) -> list[tuple[int, list[float]]]:
    rows: dict[int, list[float]] = {}
    for path in files:
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames is None:
                raise ValueError(f"{path} has no CSV header")
            required = {"global_id", *FIELDS}
            missing = required.difference(reader.fieldnames)
            if missing:
                raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
            for line, row in enumerate(reader, start=2):
                try:
                    identifier = int(row["global_id"])
                    values = [float(row[name]) for name in FIELDS]
                except (TypeError, ValueError) as error:
                    raise ValueError(f"invalid values in {path}:{line}") from error
                if identifier < 0 or not all(math.isfinite(value) for value in values):
                    raise ValueError(f"non-finite or negative id in {path}:{line}")
                if identifier in rows:
                    raise ValueError(f"duplicate global id {identifier}")
                rows[identifier] = values
    return [(identifier, rows[identifier]) for identifier in sorted(rows)]


def _write_vtk(path: Path, rows: list[tuple[int, list[float]]]) -> None:
    with path.open("w") as stream:
        stream.write("# vtk DataFile Version 3.0\n")
        stream.write("BabelSim merged owned-cell fields\nASCII\n")
        stream.write("DATASET POLYDATA\n")
        stream.write(f"POINTS {len(rows)} double\n")
        for _, values in rows:
            stream.write(f"{values[0]:.17g} {values[1]:.17g} {values[2]:.17g}\n")
        stream.write(f"VERTICES {len(rows)} {2 * len(rows)}\n")
        for point in range(len(rows)):
            stream.write(f"1 {point}\n")
        stream.write(f"POINT_DATA {len(rows)}\n")
        stream.write("SCALARS global_id int 1\nLOOKUP_TABLE default\n")
        stream.write("".join(f"{identifier}\n" for identifier, _ in rows))
        stream.write("VECTORS velocity double\n")
        for _, values in rows:
            stream.write(f"{values[3]:.17g} {values[4]:.17g} {values[5]:.17g}\n")
        stream.write("SCALARS pressure double 1\nLOOKUP_TABLE default\n")
        stream.write("".join(f"{values[6]:.17g}\n" for _, values in rows))


def _write_tecplot(path: Path, title: str, rows: list[tuple[int, list[float]]]) -> None:
    with path.open("w") as stream:
        stream.write(
            'TITLE="BabelSim merged owned-cell fields"\n'
            'VARIABLES="X","Y","Z","U","V","W","P","GlobalID"\n')
        stream.write(f'ZONE T="{title}", I={len(rows)}, F=POINT\n')
        for identifier, values in rows:
            stream.write(
                " ".join(f"{value:.17g}" for value in values) +
                f" {identifier}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--prefix", default="fields")
    parser.add_argument(
        "--ranks", type=int,
        help="rank count when files use the size-qualified prefix")
    parser.add_argument(
        "--output", type=Path,
        help="output basename; defaults to <directory>/<prefix>")
    parser.add_argument(
        "--format", choices=("both", "vtk", "tecplot"), default="both")
    args = parser.parse_args()
    if args.ranks is not None and args.ranks <= 0:
        parser.error("--ranks must be positive")

    files = _rank_files(args.directory, args.prefix, args.ranks)
    rows = _load(files)
    if not rows:
        raise ValueError("rank CSV files contain no owned cells")
    basename = args.output or args.directory / args.prefix
    basename.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    if args.format in ("both", "vtk"):
        vtk = Path(str(basename) + ".vtk")
        _write_vtk(vtk, rows)
        written.append(vtk)
    if args.format in ("both", "tecplot"):
        tecplot = Path(str(basename) + ".dat")
        _write_tecplot(tecplot, args.prefix, rows)
        written.append(tecplot)
    print(
        f"ranks={len(files)} cells={len(rows)} " +
        " ".join(f"{path.suffix[1:]}={path}" for path in written))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        raise SystemExit(f"merge_parallel_csv.py: {error}")
