#!/usr/bin/env python3
"""Build and execute the documented tutorials; write logs and a JSON evidence report."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "build/docs-refresh")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="tutorials-", dir=output))
    records: list[dict[str, object]] = []
    env = dict(os.environ, TMPDIR="/tmp", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

    def run(label: str, command: list[object], expected: int = 0) -> str:
        command = [str(value) for value in command]
        result = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                capture_output=True, timeout=240)
        log = work / f"{label}.log"
        log.write_text(result.stdout + result.stderr)
        records.append(dict(check=label, command=command, exitCode=result.returncode,
                            expectedExitCode=expected, log=str(log)))
        if result.returncode != expected:
            raise RuntimeError(f"{label}: exit {result.returncode}, expected {expected}; see {log}")
        print(f"PASS {label}", flush=True)
        return result.stdout

    def load(directory: Path, name: str, count: int) -> dict[int, list[float]]:
        values: dict[int, list[float]] = {}
        for path in directory.glob(f"rank-*/{name}.csv"):
            with path.open() as stream:
                for row in csv.DictReader(stream):
                    identifier = int(row["global_id"])
                    assert identifier not in values, (path, identifier)
                    components = [float(row[key]) for key in row if key.startswith("value")]
                    assert components and all(math.isfinite(value) for value in components), path
                    values[identifier] = components
        assert set(values) == set(range(count)), (directory, name, len(values))
        return values

    report: dict[str, object] = dict(createdUTC=datetime.now(timezone.utc).isoformat(),
                                   workDirectory=str(work), scope="tutorial execution, API linkage and I/O; no physical validation")
    try:
        run("build-examples", ["make", "-f", "docs/examples/Makefile", "-j4"])
        run("build-sdk-project", ["make", "-f", "docs/examples/Makefile", "sdk-example"])
        solve = ROOT / "build/babelsim-solve"
        cases: dict[str, Path] = {}
        for template in ("heat", "transport", "simple", "piso", "custom-heat", "coupled"):
            case = work / template
            cases[template] = case
            run(f"create-{template}", [sys.executable, ROOT / "docs/examples/create_case.py",
                                      "--template", template, "--output", case])
            program = (ROOT / "build/examples/custom-heat" if template == "custom-heat" else
                       ROOT / "build/examples/coupled-scalar" if template == "coupled" else solve)
            for ranks in (1, 2, 4):
                command: list[object] = [program, "-case", case, "-time", f"np{ranks}"]
                if ranks > 1:
                    command = ["mpirun", "-np", ranks] + command
                run(f"run-{template}-np{ranks}", command)
                snapshot = case / "results" / f"np{ranks}"
                fields = ("T", "C", "total") if template == "coupled" else (
                    ("U", "p") if template in {"simple", "piso"} else
                    ("C", "U") if template == "transport" else ("T",))
                for name in fields:
                    load(snapshot, name, 256 if template in {"simple", "piso"} else 128)
                if template != "simple":
                    assert sorted(p.name for p in snapshot.iterdir() if p.is_dir() and re.match(r"^0\.", p.name)) == [
                        "0.01", "0.02", "0.03", "0.04", "0.05"], snapshot
                if ranks > 1:
                    tolerance = "1e-7" if template in {"simple", "piso"} else "1e-8"
                    run(f"compare-{template}-np{ranks}", [sys.executable, ROOT / "tools/compare_parallel_results.py",
                        case / "results/np1", snapshot, "--atol", tolerance, "--rtol", tolerance])
        heat = cases["heat"]
        run("sdk-custom-heat", [ROOT / "build/examples/sdk-project/custom-heat", "-case", cases["custom-heat"], "-time", "sdk"])
        run("compare-sdk-heat", [sys.executable, ROOT / "tools/compare_parallel_results.py",
            cases["custom-heat"] / "results/np1", cases["custom-heat"] / "results/sdk", "--atol", "1e-12", "--rtol", "1e-12"])
        for ranks in (2, 4):
            run(f"sdk-custom-heat-np{ranks}", ["mpirun", "-np", ranks, ROOT / "build/examples/sdk-project/custom-heat",
                "-case", cases["custom-heat"], "-time", f"sdk{ranks}"])
            run(f"compare-sdk-heat-np{ranks}", [sys.executable, ROOT / "tools/compare_parallel_results.py",
                cases["custom-heat"] / "results/sdk", cases["custom-heat"] / f"results/sdk{ranks}", "--atol", "1e-8", "--rtol", "1e-8"])
        run("registered-custom-heat", [ROOT / "build/examples/registered-heat", "-case", cases["custom-heat"], "-time", "registered"])
        run("registered-built-in-heat", [ROOT / "build/examples/registered-heat", "-case", heat, "-time", "registered"])
        run("compare-public-heat", [sys.executable, ROOT / "tools/compare_parallel_results.py",
            heat / "results/np1", cases["custom-heat"] / "results/np1", "--atol", "1e-12", "--rtol", "1e-12"])
        for template in ("heat", "custom-heat"):
            methods = cases[template] / "numerics/methods.bs"
            methods.write_text(methods.read_text().replace("time euler", "time bdf2"))
            program = solve if template == "heat" else ROOT / "build/examples/custom-heat"
            run(f"bdf2-{template}", [program, "-case", cases[template], "-time", "bdf2"])
            methods.write_text(methods.read_text().replace("time bdf2", "time euler"))
        run("compare-bdf2-public-heat", [sys.executable, ROOT / "tools/compare_parallel_results.py",
            heat / "results/bdf2", cases["custom-heat"] / "results/bdf2", "--atol", "1e-12", "--rtol", "1e-12"])
        run("post-series", [ROOT / "build/babelsim-post", "-case", heat, "-time", "np1/all", "-format", "vtk"])
        pvd = ET.parse(heat / "post/np1/series.pvd")
        assert [float(node.attrib["timestep"]) for node in pvd.findall(".//DataSet")] == [.01, .02, .03, .04, .05]
        for node in pvd.findall(".//DataSet"):
            tree = ET.parse(heat / "post/np1" / node.attrib["file"])
            assert tree.find(".//Piece").attrib["NumberOfCells"] == "128"
        run("post-latest", [ROOT / "build/babelsim-post", "-case", heat, "-time", "np1/latest", "-format", "vtk", "tecplot"])
        assert (heat / "post/np1/0.05.dat").is_file()
        restart = work / "restart"
        run("create-restart", [sys.executable, ROOT / "docs/examples/create_case.py", "--template", "heat", "--output", restart])
        initial = restart / "fields/initial/T.field"
        relative = os.path.relpath(heat / "results/np1/0.03", initial.parent)
        initial.write_text(initial.read_text().replace("internal uniform (0)", f"internal file {relative}"))
        (restart / "control.bs").write_text("startTime 0.03\nendTime 0.05\ndeltaT 0.01\n")
        run("restart-np2", ["mpirun", "-np", 2, solve, "-case", restart, "-time", "continued"])
        run("compare-restart", [sys.executable, ROOT / "tools/compare_parallel_results.py",
            heat / "results/np1/0.05", restart / "results/continued/0.05", "--atol", "1e-8", "--rtol", "1e-8"])
        single = work / "one-cell"
        run("create-one-cell", [sys.executable, ROOT / "docs/examples/create_case.py", "--template", "heat", "--output", single])
        shutil.copyfile(ROOT / "docs/examples/one-cell.mesh", single / "mesh/domain.mesh")
        run("read-one-cell-mesh", [solve, "-case", single, "-time", "meshFormat"])
        load(single / "results/meshFormat", "T", 1)
        inputs = sorted(path for path in heat.rglob("*") if path.is_file() and "results" not in path.parts and "post" not in path.parts)
        before = {path: path.read_bytes() for path in inputs}
        run("refuse-overwrite", [sys.executable, ROOT / "docs/examples/create_case.py", "--template", "heat", "--output", heat], expected=1)
        assert all(path.read_bytes() == value for path, value in before.items())
        dangling = work / "dangling-case"
        missing = work / "must-remain-missing"
        dangling.symlink_to(missing, target_is_directory=True)
        run("refuse-dangling-symlink", [sys.executable, ROOT / "docs/examples/create_case.py", "--template", "heat", "--output", dangling], expected=1)
        assert dangling.is_symlink() and not missing.exists()
        (single / "case.bs").write_text((single / "case.bs").read_text().replace("solver heat", "solver unknownTutorial"))
        run("unknown-solver", [solve, "-case", single], expected=1)
        T = load(cases["coupled"] / "results/np1", "T", 128)
        C = load(cases["coupled"] / "results/np1", "C", 128)
        total = load(cases["coupled"] / "results/np1", "total", 128)
        assert max(abs(total[i][0] - T[i][0] - C[i][0]) for i in total) < 1e-12
        report.update(status="passed", checks=records, tutorialRuns=18, mpiRanks=[1, 2, 4],
                      assertions=["finite unique owned cells", "five saved timesteps", "serial/MPI comparison",
                                  "public/built-in Euler and BDF2 agreement", "VTK geometry and PVD times",
                                  "Euler field restart across rank counts", "one-cell native reader",
                                  "no-overwrite policy", "unknown solver error", "derived-field arithmetic"])
    except (AssertionError, RuntimeError, OSError, subprocess.TimeoutExpired) as error:
        report.update(status="failed", error=str(error), checks=records)
    (output / "validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(f"{report['status']}: {output / 'validation.json'}", flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
