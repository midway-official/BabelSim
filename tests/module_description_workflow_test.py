#!/usr/bin/env python3
"""Check public metadata CLI and execute its generated scalar tutorials."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile


def run(command: list[str], ok: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if ok and result.returncode != 0:
        raise AssertionError(f"{command} failed ({result.returncode})\n{result.stdout}\n{result.stderr}")
    if not ok and result.returncode == 0:
        raise AssertionError(f"{command} accepted an invalid request")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", required=True, type=Path)
    parser.add_argument("--external", type=Path)
    args = parser.parse_args()
    solver = str(args.solver.resolve())
    listed = run([solver, "-list"]).stdout
    for name in ("heat", "transport", "simple", "transientSimple", "piso", "SA", "kOmega", "kEpsilon"):
        assert f"{name} (" in listed, (name, listed)
        schema = json.loads(run([solver, "-describe", name, "-format", "json"]).stdout)
        assert schema["title"] == name
        assert schema["x-boundaryTypes"], "static initialization lost declared boundaries"
        assert schema["x-equations"], name
        assert schema["x-timeMethods"], name
        if name in ("heat", "transport"):
            assert schema["required"] == ["physics"]
        markdown = run([solver, "-describe", name]).stdout
        assert f"# {name}" in markdown
        for dictionary in schema["properties"].values():
            for parameter in dictionary["properties"]:
                assert parameter in markdown
    run([solver, "-describe", "unknown"], ok=False)
    run([solver, "-describe", "heat", "-case", "cases/heat"], ok=False)
    with tempfile.TemporaryDirectory(prefix="babelsim-description-workflow-") as directory:
        work = Path(directory)
        for name, field in (("heat", "T"), ("transport", "C")):
            case = work / name
            run([solver, "-template", name, "-output", str(case)])
            assert (case / "README.md").is_file()
            schema = json.loads((case / "module.schema.json").read_text())
            assert schema["title"] == name
            original = (case / "physics/model.bs").read_bytes()
            run([solver, "-template", name, "-output", str(case)], ok=False)
            assert (case / "physics/model.bs").read_bytes() == original
            run([solver, "-case", str(case), "-time", "generated"])
            with (case / "results/generated/rank-0000" / f"{field}.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            assert len(rows) == 1
            value = float(rows[0]["value0"])
            assert abs(value - 0.05) < 1e-10, (name, value)
            methods = case / "numerics/methods.bs"
            methods.write_text(methods.read_text().replace("time euler", "time bdf2"))
            run([solver, "-case", str(case), "-time", "generatedBdf2"])
            with (case / "results/generatedBdf2/rank-0000" / f"{field}.csv").open() as stream:
                value = float(next(csv.DictReader(stream))["value0"])
            assert abs(value - 0.05) < 1e-10, (name, "bdf2", value)
        run([solver, "-template", "piso", "-output", str(work / "unsupported")], ok=False)
        assert not (work / "unsupported").exists()
        if args.external:
            external = str(args.external.resolve())
            schema = json.loads(run([external, "-describe", "customHeat", "-format", "json"]).stdout)
            assert schema["title"] == "customHeat"
            case = work / "external"
            run([external, "-template", "customHeat", "-output", str(case)])
            run([external, "-case", str(case), "-time", "generated"])
    print("module_description_workflow_test: discovery, schema/docs, no-overwrite templates, Euler/BDF2 tutorials passed")


if __name__ == "__main__":
    main()
