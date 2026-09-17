#!/usr/bin/env python3
"""Run reproducible BabelSim backend experiments and collect structured evidence.

The driver deliberately treats the solver as a black box.  It never changes a
case's numerical dictionaries and never decides physical convergence.  The
application reports the numerical status; this script records it, classifies
process failures, and computes wall-clock statistics from independent runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import statistics
import subprocess
import sys
import threading
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def cpu_topology() -> dict[str, Any]:
    """Return CPUs visible to this process and their physical-core mapping."""
    allowed = set(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    records: list[tuple[int, int, int]] = []
    topology = run_text(["lscpu", "-p=CPU,CORE,SOCKET"])
    for line in topology.splitlines():
        if not line or line.startswith("#"):
            continue
        try:
            cpu, core, socket = (int(value) for value in line.split(",")[:3])
        except (TypeError, ValueError):
            continue
        if allowed is None or cpu in allowed:
            records.append((cpu, core, socket))
    if not records:
        logical = len(allowed) if allowed is not None else (os.cpu_count() or 1)
        return {
            "allowedCpus": sorted(allowed) if allowed is not None else [],
            "logicalCpus": logical,
            "physicalCores": logical,
        }
    return {
        "allowedCpus": sorted(cpu for cpu, _, _ in records),
        "logicalCpus": len(records),
        "physicalCores": len({(core, socket) for _, core, socket in records}),
    }


def run_text(command: list[str], *, timeout: float = 30.0) -> str:
    try:
        result = subprocess.run(
            command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired) as error:
        return f"<unavailable: {error}>"
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def case_manifest(case: Path) -> dict[str, Any]:
    """Hash every input file so a result cannot be detached from its case."""
    output_roots = {"results"}
    output_control = case / "output.bs"
    if output_control.exists():
        for line in output_control.read_text().splitlines():
            tokens = line.split()
            if len(tokens) == 2 and tokens[0] == "directory":
                output_roots.add(tokens[1].strip("/"))
    files = []
    aggregate = hashlib.sha256()
    for path in sorted(p for p in case.rglob("*") if p.is_file()):
        relative = path.relative_to(case).as_posix()
        if any(relative == root or relative.startswith(root + "/") for root in output_roots):
            continue
        digest = sha256(path)
        size = path.stat().st_size
        files.append({"path": relative, "bytes": size, "sha256": digest})
        aggregate.update(relative.encode("utf-8"))
        aggregate.update(b"\0")
        aggregate.update(digest.encode("ascii"))
        aggregate.update(b"\n")
    return {"sha256": aggregate.hexdigest(), "files": files}


def result_manifest(case: Path, run_name: str) -> list[dict[str, Any]]:
    """Return checksums for files written below the run-specific output tree."""
    results = []
    for path in sorted(p for p in case.rglob("*") if p.is_file() and run_name in p.parts):
        results.append({
            "path": path.relative_to(case).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    return results


def git_metadata() -> dict[str, Any]:
    return {
        "commit": run_text(["git", "rev-parse", "HEAD"]),
        "branch": run_text(["git", "branch", "--show-current"]),
        "status": run_text(["git", "status", "--short"]),
        "diff": run_text(["git", "diff", "--binary"]),
    }


def host_metadata() -> dict[str, Any]:
    cpuinfo = Path("/proc/cpuinfo")
    topology = cpu_topology()
    return {
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "cpu": run_text(["lscpu"]),
        "memory": run_text(["free", "-h"]),
        "mpi": run_text(["mpirun", "--version"]),
        "compiler": run_text(["g++", "--version"]),
        "mpi_compiler": run_text(["mpic++", "--showme:command"]),
        "eigen": run_text([
            "bash", "-lc",
            "awk '/EIGEN_WORLD_VERSION|EIGEN_MAJOR_VERSION|EIGEN_MINOR_VERSION/ {print}' "
            "/usr/include/eigen3/Eigen/src/Core/util/Macros.h"
        ]),
        "cpu_flags": (
            next((line.split(":", 1)[1].strip() for line in cpuinfo.read_text().splitlines()
                  if line.lower().startswith("flags")), "<unavailable>")
            if cpuinfo.exists() else "<unavailable>"
        ),
        "cpuTopology": topology,
        "environment": {
            name: os.environ.get(name, "")
            for name in ("CXXFLAGS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                         "MKL_NUM_THREADS", "GOMP_CPU_AFFINITY", "OMP_PROC_BIND")
        },
    }


def build_metadata() -> dict[str, Any]:
    """Record the repository build recipe alongside the binary hash."""
    makefile = ROOT / "Makefile"
    try:
        text = makefile.read_text()
    except OSError as error:
        return {"makefile": f"<unavailable: {error}>"}

    def assignment(name: str) -> str:
        lines = text.splitlines()
        for index, line in enumerate(lines):
            if not line.startswith(name) or "?=" not in line:
                continue
            parts = [line.split("?=", 1)[1].strip().rstrip("\\").strip()]
            while parts and index + len(parts) < len(lines):
                continuation = lines[index + len(parts)]
                if not continuation.lstrip().startswith("-") and not continuation.startswith(" "):
                    break
                previous = lines[index + len(parts) - 1]
                if not previous.rstrip().endswith("\\"):
                    break
                parts.append(continuation.strip().rstrip("\\").strip())
            return " ".join(parts)
        return ""

    return {
        "makefileSha256": sha256(makefile),
        "optimizationFlags": assignment("OPTFLAGS"),
        "compilerFlags": assignment("CXXFLAGS"),
        "backendDefaults": {
            name: assignment(name)
            for name in ("ASYNC_HALO", "CSR_SPMV", "SERIAL_CSR_SPMV")
        },
    }


def binding_options(rank: int, topology: dict[str, Any]) -> tuple[list[str], str]:
    """Select a non-oversubscribed binding policy for physical cores or SMT."""
    if rank <= 0 or rank > int(topology["logicalCpus"]):
        raise ValueError(
            f"rank count {rank} exceeds visible logical CPUs {topology['logicalCpus']}")
    if rank <= int(topology["physicalCores"]):
        return ["--bind-to", "core", "--map-by", "core"], "physical-core"
    return ["--bind-to", "hwthread", "--map-by", "hwthread"], "smt-hwthread"


def parse_ranks(value: str) -> list[int]:
    ranks: list[int] = []
    for token in value.split(","):
        try:
            rank = int(token)
        except ValueError as error:
            raise argparse.ArgumentTypeError(f"invalid rank count: {token}") from error
        if rank <= 0:
            raise argparse.ArgumentTypeError("rank counts must be positive")
        ranks.append(rank)
    if not ranks:
        raise argparse.ArgumentTypeError("at least one rank count is required")
    return list(dict.fromkeys(ranks))


def terminate_group(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=5)
    except (ProcessLookupError, subprocess.TimeoutExpired):
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=5)


class ProcessTreeSampler:
    """Sample concurrent RSS for mpirun and its solver descendants.

    Linux exposes the process tree and VmRSS without requiring ptrace.  The
    sampler is intentionally best-effort: missing /proc entries are expected
    while MPI ranks exit.  It never participates in solver decisions.
    """

    def __init__(self, root_pid: int, solver_name: str, interval: float = 0.1) -> None:
        self.root_pid = root_pid
        self.solver_name = solver_name
        self.interval = interval
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._lock = threading.Lock()
        self.peak_tree_kib = 0
        self.peak_solver_kib = 0

    @staticmethod
    def _children(pid: int) -> list[int]:
        try:
            text = Path(f"/proc/{pid}/task/{pid}/children").read_text()
        except (OSError, ValueError):
            return []
        values: list[int] = []
        for token in text.split():
            try:
                values.append(int(token))
            except ValueError:
                continue
        return values

    @staticmethod
    def _rss_kib(pid: int) -> int:
        try:
            for line in Path(f"/proc/{pid}/status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
        except (OSError, ValueError, IndexError):
            pass
        return 0

    @staticmethod
    def _command_line(pid: int) -> str:
        try:
            return Path(f"/proc/{pid}/cmdline").read_bytes().replace(b"\0", b" ").decode(
                "utf-8", errors="replace")
        except OSError:
            return ""

    def _sample_once(self) -> None:
        pending = [self.root_pid]
        seen: set[int] = set()
        tree_rss = 0
        solver_rss = 0
        while pending:
            pid = pending.pop()
            if pid in seen:
                continue
            seen.add(pid)
            rss = self._rss_kib(pid)
            tree_rss += rss
            if self.solver_name in self._command_line(pid):
                solver_rss = max(solver_rss, rss)
            pending.extend(self._children(pid))
        with self._lock:
            self.peak_tree_kib = max(self.peak_tree_kib, tree_rss)
            self.peak_solver_kib = max(self.peak_solver_kib, solver_rss)

    def _run(self) -> None:
        while not self._stop.is_set():
            self._sample_once()
            self._stop.wait(self.interval)
        self._sample_once()

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> dict[str, int]:
        self._stop.set()
        self._thread.join(timeout=2.0)
        with self._lock:
            return {
                "peakTreeRssKiB": self.peak_tree_kib,
                "peakSolverRssKiB": self.peak_solver_kib,
            }


def status_from_process(returncode: int | None, timed_out: bool, reports: list[dict[str, Any]]) -> str:
    if timed_out:
        return "timeout"
    if returncode is None or returncode != 0:
        if not reports:
            return "infrastructureFailure"
        statuses = {report.get("status") for report in reports}
        if statuses and statuses <= {"maxIterations"}:
            return "maxIterations"
        return "numericalFailure"
    statuses = {report.get("status") for report in reports}
    if "converged" in statuses and statuses <= {"converged"}:
        return "converged"
    if "maxIterations" in statuses:
        return "maxIterations"
    if "numericalFailure" in statuses:
        return "numericalFailure"
    return "missingPerformanceData"


def read_reports(directory: Path) -> list[dict[str, Any]]:
    reports = []
    for path in sorted(directory.glob("rank-*.json")):
        try:
            value = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            raise RuntimeError(f"invalid performance report {path}: {error}") from error
        value["path"] = str(path)
        reports.append(value)
    return reports


def run_one(
    args: argparse.Namespace,
    output: Path,
    rank: int,
    repeat: int,
    *,
    warmup: bool = False,
) -> dict[str, Any]:
    run_label = f"warmup-{rank:02d}" if warmup else f"repeat-{repeat:02d}"
    run_dir = output / f"ranks-{rank}" / run_label
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest = run_dir / "run.json"
    if args.resume and manifest.exists():
        try:
            previous = json.loads(manifest.read_text())
            if previous.get("status") in {"converged", "maxIterations", "numericalFailure", "timeout"}:
                previous["resumed"] = True
                return previous
        except json.JSONDecodeError:
            pass

    performance_dir = run_dir / "performance"
    run_name = f"benchmark-r{rank}-{'warmup' if warmup else f'n{repeat:02d}'}"
    topology = cpu_topology()
    binding, binding_name = binding_options(rank, topology)
    command = [
        "mpirun", *binding, "--report-bindings", "-np", str(rank), str(args.binary),
        "-case", str(args.case), "-time", run_name,
        "-performance", str(performance_dir),
    ]
    if args.mode == "throughput":
        # The production solver still owns its stopping rule.  This label
        # only changes the evidence classification; it never fakes convergence.
        mode_note = "fixed-production-run"
    else:
        mode_note = "production-convergence"
    environment = dict(os.environ)
    environment.setdefault("TMPDIR", "/tmp")
    tracked_environment = {
        name: environment.get(name, "")
        for name in (
            "TMPDIR", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "GOMP_CPU_AFFINITY", "OMP_PROC_BIND", "BABELSIM_ASYNC_KRYLOV_HALO",
            "BABELSIM_CSR_SPMV", "BABELSIM_SERIAL_CSR_SPMV",
        )
    }
    start = time.monotonic()
    timed_out = False
    sampler: ProcessTreeSampler | None = None
    try:
        process = subprocess.Popen(
            command, cwd=ROOT, env=environment, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            start_new_session=True)
        sampler = ProcessTreeSampler(process.pid, args.binary.name)
        sampler.start()
        try:
            stdout, _ = process.communicate(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            terminate_group(process)
            stdout, _ = process.communicate()
    except OSError as error:
        process = None
        stdout = f"cannot start benchmark process: {error}\n"
    memory = sampler.stop() if sampler is not None else {
        "peakTreeRssKiB": 0,
        "peakSolverRssKiB": 0,
    }
    elapsed = time.monotonic() - start
    (run_dir / "stdout.log").write_text(stdout)
    reports = read_reports(performance_dir) if performance_dir.exists() else []
    if not reports and not timed_out and process is not None and process.returncode == 0:
        raise RuntimeError(f"solver exited successfully without performance data: {run_dir}")
    status = status_from_process(
        None if process is None else process.returncode, timed_out, reports)
    result = {
        "case": str(args.case),
        "binary": str(args.binary),
        "ranks": rank,
        "repeat": repeat,
        "warmup": warmup,
        "binding": binding_name,
        "cpuTopology": topology,
        "mode": args.mode,
        "modeNote": mode_note,
        "command": command,
        "environment": tracked_environment,
        "elapsedSeconds": elapsed,
        "returncode": None if process is None else process.returncode,
        "status": status,
        "memory": memory,
        "performance": reports,
        "results": result_manifest(args.case, run_name),
    }
    manifest.write_text(json.dumps(result, indent=2) + "\n")
    return result


def write_summary(
    args: argparse.Namespace,
    output: Path,
    runs: list[dict[str, Any]],
    input_manifest: dict[str, Any],
) -> None:
    measured_runs = [run for run in runs if not run.get("warmup", False)]
    by_rank: dict[str, list[float]] = {}
    all_by_rank: dict[str, dict[str, list[float]]] = {}
    for run in measured_runs:
        all_by_rank.setdefault(str(run["ranks"]), {}).setdefault(
            run["status"], []).append(run["elapsedSeconds"])
        if run["status"] == "converged":
            by_rank.setdefault(str(run["ranks"]), []).append(run["elapsedSeconds"])
    statistics_by_rank = {}
    for rank, values in by_rank.items():
        statistics_by_rank[rank] = {
            "samples": len(values),
            "medianSeconds": statistics.median(values),
            "minSeconds": min(values),
            "maxSeconds": max(values),
            "coefficientOfVariation": (
                statistics.pstdev(values) / statistics.mean(values)
                if len(values) > 1 and statistics.mean(values) else 0.0
            ),
        }
    counter_names = (
        "linearSolves", "krylovIterations", "sparseMatvecs", "haloExchanges", "haloBytes",
        "globalReductions", "equationAssemblies", "preconditionerSetups",
        "preconditionerApplications", "outputWrites", "elapsedSeconds", "assemblySeconds",
        "preconditionerSeconds", "preconditionerApplySeconds", "linearSolveSeconds",
        "sparseMatvecSeconds", "haloSeconds", "globalReductionSeconds", "outputSeconds",
            )

    wall_clock_by_rank = {}
    for rank, statuses in all_by_rank.items():
        wall_clock_by_rank[rank] = {}
        for status, values in statuses.items():
            mean = statistics.mean(values)
            wall_clock_by_rank[rank][status] = {
                "samples": len(values),
                "medianSeconds": statistics.median(values),
                "minSeconds": min(values),
                "maxSeconds": max(values),
                "coefficientOfVariation": (
                    statistics.pstdev(values) / mean
                    if len(values) > 1 and mean else 0.0
                ),
            }
    # Keep the requested MPI size separate from the rank-local index.  A
    # single "rank-0" bucket across -np 1,2,4,... would mix unrelated runs.
    counters_by_rank: dict[str, dict[str, dict[str, dict[str, list[float]]]]] = {}
    partition_by_rank: dict[str, dict[str, dict[str, dict[str, list[float]]]]] = {}
    phase_names = (
        "caseSetupSeconds", "solverSeconds", "solverComputeSeconds",
        "applicationSeconds",
    )
    phase_times_by_rank: dict[str, dict[str, dict[str, dict[str, list[float]]]]] = {}
    critical_phase_samples: dict[str, dict[str, list[float]]] = {}
    for run in measured_runs:
        size_key = str(run["ranks"])
        run_reports = run.get("performance", [])
        critical = critical_phase_samples.setdefault(size_key, {})
        for name in phase_names:
            values = [report[name] for report in run_reports if name in report]
            if values:
                critical.setdefault(name, []).append(max(values))
        for report in run_reports:
            local = report.get("local", report)
            rank_key = str(report.get("rank", 0))
            aggregate = counters_by_rank.setdefault(size_key, {}).setdefault(rank_key, {})
            for name in counter_names:
                if name in local:
                    aggregate.setdefault(name, {}).setdefault("samples", []).append(local[name])
            phases = phase_times_by_rank.setdefault(size_key, {}).setdefault(rank_key, {})
            for name in phase_names:
                if name in report:
                    phases.setdefault(name, {}).setdefault("samples", []).append(report[name])
            partition = report.get("partition")
            if partition:
                metrics = partition_by_rank.setdefault(size_key, {}).setdefault(rank_key, {})
                for name, value in partition.items():
                    metrics.setdefault(name, {}).setdefault("samples", []).append(value)
    for rank_groups in counters_by_rank.values():
        for aggregate in rank_groups.values():
            for name, value in list(aggregate.items()):
                samples = value.pop("samples")
                value.update({"min": min(samples), "mean": statistics.mean(samples), "max": max(samples)})
    for rank_groups in phase_times_by_rank.values():
        for phases in rank_groups.values():
            for name, value in list(phases.items()):
                samples = value.pop("samples")
                value.update({"min": min(samples), "mean": statistics.mean(samples), "max": max(samples)})
    for rank_groups in partition_by_rank.values():
        for metrics in rank_groups.values():
            for name, value in list(metrics.items()):
                samples = value.pop("samples")
                value.update({"min": min(samples), "mean": statistics.mean(samples), "max": max(samples)})
    critical_path_by_rank = {}
    for size_key, phases in critical_phase_samples.items():
        critical_path_by_rank[size_key] = {}
        for name, samples in phases.items():
            critical_path_by_rank[size_key][name] = {
                "samples": len(samples),
                "min": min(samples),
                "mean": statistics.mean(samples),
                "max": max(samples),
            }
    summary = {
        "case": str(args.case),
        "caseManifest": input_manifest,
        "binary": str(args.binary),
        "binarySha256": sha256(args.binary),
        "mode": args.mode,
        "timeoutSeconds": args.timeout,
        "repeat": args.repeat,
        "warmup": args.warmup,
        "git": git_metadata(),
        "host": host_metadata(),
        "build": build_metadata(),
        "runs": runs,
        "wallClockByRank": wall_clock_by_rank,
        "convergedWallClock": statistics_by_rank,
        "criticalPathByRank": critical_path_by_rank,
        "phaseTimesByRank": phase_times_by_rank,
        "partitionByRank": partition_by_rank,
        "localCountersByRank": counters_by_rank,
    }
    (output / "metadata.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (output / "results.csv").open("w") as csv:
        csv.write("case,ranks,repeat,mode,status,elapsed_seconds,returncode\n")
        for run in runs:
            csv.write(",".join([
                str(args.case), str(run["ranks"]), str(run["repeat"]),
                run["mode"], run["status"], f'{run["elapsedSeconds"]:.9g}',
                str(run["returncode"]),
            ]) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--ranks", type=parse_ranks, default=[1])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument(
        "--warmup", type=int, default=1,
        help="number of uncounted warmup runs per rank (default: 1)")
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--mode", choices=("complete", "throughput"), default="complete")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--binary", type=Path, default=ROOT / "build/babelsim-solve")
    args = parser.parse_args()
    args.case = args.case.resolve()
    args.output = args.output.resolve()
    args.binary = args.binary.resolve()
    if not args.case.is_dir():
        parser.error(f"case directory does not exist: {args.case}")
    if not args.binary.is_file():
        parser.error(f"solver binary does not exist: {args.binary}")
    if args.repeat <= 0 or args.warmup < 0 or args.timeout <= 0:
        parser.error("repeat and timeout must be positive; warmup cannot be negative")
    args.output.mkdir(parents=True, exist_ok=True)
    input_manifest = case_manifest(args.case)
    runs = []
    for rank in args.ranks:
        for warmup in range(args.warmup):
            print(f"benchmark warmup case={args.case} ranks={rank} run={warmup + 1}", flush=True)
            runs.append(run_one(args, args.output, rank, 0, warmup=True))
            write_summary(args, args.output, runs, input_manifest)
        for repeat in range(1, args.repeat + 1):
            print(f"benchmark case={args.case} ranks={rank} repeat={repeat}", flush=True)
            runs.append(run_one(args, args.output, rank, repeat))
            write_summary(args, args.output, runs, input_manifest)
    print(f"wrote {args.output / 'metadata.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
