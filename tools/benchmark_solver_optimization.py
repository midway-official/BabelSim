#!/usr/bin/env python3
"""Sequential binary A/B and strong-scaling measurements of production solvers.

Use identical input dictionaries for every binary and rank count. Fixed SIMPLE
iterations are throughput evidence; accepted PISO steps are transient evidence.
Only max-rank iterationComputeSeconds determines speedup. This interval starts
after the first completed iteration, excluding setup and its first solve. Output
within the interval is subtracted per rank. A machine-wide advisory
lock prevents concurrent invocations of this driver. Repetitions are interleaved;
all runs, including warmups and failures, retain their logs and rank-local JSON.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import random
import re
import shutil
import statistics

import benchmark_backend as bench


def named(items: list[str]) -> dict[str, Path]:
    result = {}
    for item in items:
        name, path = item.split("=", 1)
        if not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name in result:
            raise ValueError(f"invalid or duplicate label: {name}")
        result[name] = Path(path).resolve()
    return result


def summarize(records: list[dict], output: Path) -> None:
    groups = {}
    for run in records:
        if not run["warmup"] and run["validWork"] and run.get("measurementValid", True):
            groups.setdefault((run["caseLabel"], run["variant"], run["ranks"]), []).append(run)
    rows = []
    for (case, variant, ranks), runs in sorted(groups.items()):
        metric = lambda key: statistics.median(max(p["local"][key] for p in r["performance"]) for r in runs)
        values = [r["computeSeconds"] for r in runs]
        measured = statistics.median(values)
        serial = groups.get((case, variant, 1))
        speedup = statistics.median(r["computeSeconds"] for r in serial) / measured if serial else None
        rows.append(dict(case=case, variant=variant, ranks=ranks, samples=len(runs),
                         sustainedSamples=sum(r["sustainedSample"] for r in runs),
                         computeMedian=measured, computeMin=min(values), computeMax=max(values),
                         speedup=speedup, efficiency=speedup / ranks if speedup is not None else None,
                         workMedian=statistics.median(r["workUnits"] for r in runs),
                         setups=metric("preconditionerSetups"), iterations=metric("krylovIterations"),
                         setupSeconds=metric("preconditionerSeconds"), solveSeconds=metric("linearSolveSeconds"),
                         haloSeconds=metric("haloSeconds"), reductionSeconds=metric("globalReductionSeconds")))
    (output / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    if rows:
        with (output / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def archive_fields(run: dict, target: Path, index: Path) -> dict:
    """Retain exported owned-cell fields after the timed process has exited.

    Copies are independent of the work directory, including any future rerun.
    An unavailable old temporary result remains explicitly unavailable.
    """
    if index.exists():
        previous = json.loads(index.read_text())
        for entry in previous["files"]:
            path = Path(previous["root"]) / entry["path"]
            if not path.is_file() or bench.sha256(path) != entry["sha256"]:
                raise ValueError(f"archived field changed: {path}")
        return previous
    fields = [entry for entry in run["results"]
              if any(part.startswith("rank-") for part in Path(entry["path"]).parts)
              and (entry["path"].endswith(".csv") or Path(entry["path"]).name == "metadata.bs")]
    available = all((Path(run["case"]) / entry["path"]).is_file() for entry in fields)
    if not available:
        return dict(root=str(target), files=[], status="original temporary fields unavailable")
    for entry in fields:
        source = Path(run["case"]) / entry["path"]
        destination = target / entry["path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        if bench.sha256(destination) != entry["sha256"]:
            raise ValueError(f"exported field differs from run manifest: {source}")
    result = dict(root=str(target), files=fields,
                  status="retained" if fields else "no final field exported")
    index.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", action="append", required=True, help="label=/absolute/binary")
    parser.add_argument("--case", action="append", required=True, help="label=/absolute/case")
    parser.add_argument("--throughput", action="append", default=[], help="label=required-SIMPLE-outer-count")
    parser.add_argument("--accepted-steps", action="append", default=[], help="label=required-accepted-PISO-step-count")
    parser.add_argument("--ranks", type=bench.parse_ranks, default=[1, 2, 4, 8, 12])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--warmup", type=int, choices=(0, 1), default=1)
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--minimum-compute", type=float, default=60,
                        help="mark shorter samples; retain them without making a sustained-run claim")
    parser.add_argument("--seed", type=int, default=20261007)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True, help="isolated copies of cases and their results")
    parser.add_argument("--field-archive-root", type=Path,
                        help="optional persistent directory for checked exported CSV fields; copy after timing")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.repeat < 1 or args.minimum_compute < 0:
        parser.error("repeat must be positive and minimum-compute must be nonnegative")
    variants, cases = named(args.variant), named(args.case)
    throughput = {name: int(count) for name, count in (value.split("=", 1) for value in args.throughput)}
    steps = {name: int(count) for name, count in (value.split("=", 1) for value in args.accepted_steps)}
    if any(name not in cases or count < 1 for name, count in throughput.items()):
        parser.error("throughput must name a provided case and positive outer count")
    if any(name not in cases or name in throughput or count < 1 for name, count in steps.items()):
        parser.error("accepted-steps must name a non-throughput case and positive step count")
    args.output = args.output.resolve()
    args.work_root = args.work_root.resolve()
    if args.field_archive_root is not None:
        args.field_archive_root = args.field_archive_root.resolve()
        if args.field_archive_root.is_relative_to(args.work_root):
            parser.error("field archive must be outside the work directory")
    with open("/tmp/babelsim-solver-benchmark.lock", "w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            parser.error("another solver benchmark is running; wait for it to finish")
        args.output.mkdir(parents=True, exist_ok=True)
        args.work_root.mkdir(parents=True, exist_ok=True)
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ[key] = "1"
        if os.environ.get("PETSC_OPTIONS"):
            parser.error("clear PETSC_OPTIONS so profiling/options cannot change the A/B comparison")
        metadata = dict(host=bench.host_metadata(), build=bench.build_metadata(), git=bench.git_metadata(),
                        variants={n: {"binary": str(p), "sha256": bench.sha256(p)} for n, p in variants.items()},
                        cases={n: {"path": str(p), "manifest": bench.case_manifest(p)} for n, p in cases.items()},
                        ranks=args.ranks, repeat=args.repeat, warmup=args.warmup, throughput=throughput,
                        acceptedSteps=steps,
                        seed=args.seed, minimumCompute=args.minimum_compute,
                        metric="max_rank(iterationComputeSeconds)", workRoot=str(args.work_root))
        meta_path = args.output / "study.json"
        if meta_path.exists():
            old = json.loads(meta_path.read_text())
            keys = ("variants", "cases", "ranks", "repeat", "warmup", "throughput", "acceptedSteps", "seed", "minimumCompute", "workRoot")
            if not args.resume or any(old[key] != metadata[key] for key in keys):
                parser.error("evidence exists or inputs changed; choose a new output/work-root")
        else:
            meta_path.write_text(json.dumps(metadata, indent=2) + "\n")
        copies = {}
        for variant in variants:
            for name, source in cases.items():
                target = args.work_root / variant / name
                manifest = metadata["cases"][name]["manifest"]
                if not target.exists():
                    target.mkdir(parents=True)
                    for entry in manifest["files"]:
                        relative = Path(entry["path"])
                        dest = target / relative
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(source / relative, dest)
                        if relative.suffix in (".bs", ".field", ".json"):
                            snapshot = args.output / "inputs" / name / relative
                            snapshot.parent.mkdir(parents=True, exist_ok=True)
                            shutil.copy2(source / relative, snapshot)
                if bench.case_manifest(target)["sha256"] != manifest["sha256"]:
                    parser.error(f"case copy differs from authoritative inputs: {target}")
                copies[variant, name] = target
        jobs = [(variant, name, rank) for variant in variants for name in cases for rank in args.ranks]
        rng = random.Random(args.seed)
        records = []
        for repeat in range(0 if args.warmup else 1, args.repeat + 1):
            order = jobs.copy()
            rng.shuffle(order)
            for variant, name, rank in order:
                print(f"RUN variant={variant} case={name} ranks={rank} repeat={repeat}", flush=True)
                run_args = argparse.Namespace(case=copies[variant, name], binary=variants[variant],
                    mode="throughput" if name in throughput else "complete", resume=args.resume,
                    timeout=args.timeout, live_log=True, memory_trace=True)
                run_output = args.output / variant / name
                run = bench.run_one(run_args, run_output, rank, repeat, warmup=repeat == 0)
                if args.field_archive_root is not None:
                    run_label = f"warmup-{rank:02d}" if repeat == 0 else f"repeat-{repeat:02d}"
                    relative = Path(variant) / name / f"ranks-{rank}" / run_label
                    run["fieldArchive"] = archive_fields(run, args.field_archive_root / relative,
                                                         args.output / relative / "fields.json")
                log = (run_output / f"ranks-{rank}" / (f"warmup-{rank:02d}" if repeat == 0 else f"repeat-{repeat:02d}") / "stdout.log").read_text()
                outer = [int(n) for n in re.findall(r"^SIMPLE ([0-9]+) ", log, re.MULTILINE)]
                accepted = len(re.findall(r"^PISO [0-9]+ .*accepted=true.*converged=true", log, re.MULTILINE))
                iterations = re.findall(r"^(?:SIMPLE|PISO) [0-9]+ .*", log, re.MULTILINE)
                work = max(outer, default=0) if name in throughput else accepted
                valid = (len(run["performance"]) == rank and
                         bool(iterations) and all("linear=ok" in line for line in iterations) and
                         ((name in throughput and run["status"] == "maxIterations" and work == throughput[name]) or
                          (name not in throughput and run["status"] == "converged" and
                           (name not in steps or accepted == steps[name]))))
                if valid and work > 0:
                    valid = all(p["reportedIterationIntervals"] == work - 1 for p in run["performance"])
                run.update(variant=variant, caseLabel=name, workUnits=work, validWork=valid,
                           measurementValid=True,
                           computeSeconds=max((p["iterationComputeSeconds"] for p in run["performance"]), default=0.0))
                run["sustainedSample"] = valid and run["computeSeconds"] >= args.minimum_compute
                records.append(run)
                (args.output / "records.json").write_text(json.dumps(records, indent=2) + "\n")
                summarize(records, args.output)
                print(f"DONE status={run['status']} compute={run['computeSeconds']:.3f}s work={work} sustained={run['sustainedSample']}", flush=True)
                if not valid:
                    print("Stopping after invalid work; inspect raw evidence before continuing.", flush=True)
                    return 1
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
