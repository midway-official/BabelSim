#!/usr/bin/env python3
"""Build the bundled PETSc/Hypre without putting object files in the source tree."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], source: Path, env: dict[str, str], log: Path) -> None:
    print(f"PETSc: {command[0]} {command[1:]}\n  log: {log}", flush=True)
    with log.open("w") as stream:
        result = subprocess.run(command, cwd=source, env=env,
                                stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        print("\n".join(log.read_text(errors="replace").splitlines()[-60:]), file=sys.stderr)
        raise RuntimeError(f"PETSc command failed ({result.returncode}); see {log}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--arch", required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--cc", default="mpicc")
    parser.add_argument("--cxx", default="mpicxx")
    parser.add_argument("--fc", default="mpifort")
    parser.add_argument("--optflags", default="-O3 -march=native -mtune=native")
    parser.add_argument("--clean", action="store_true")
    args = parser.parse_args()
    if args.build.is_symlink():
        raise ValueError("BUILD must be a real directory, not a symlink")
    source, build = args.source.resolve(), args.build.resolve()
    if not re.fullmatch(r"arch-babelsim-[A-Za-z0-9_-]+", args.arch):
        raise ValueError("The bundled build requires an arch-babelsim-* architecture name")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    if ROOT == build or ROOT.is_relative_to(build) or source.is_relative_to(build):
        raise ValueError("BUILD must not contain the project or PETSc source directory")
    work = build / "petsc-build"
    prefix = build / "petsc"
    arch = source / args.arch
    if arch.is_symlink():
        if arch.resolve() != work:
            raise ValueError(f"{arch} points to another build directory; choose a different PETSC_ARCH")
    elif arch.exists():
        raise ValueError(f"Refusing to replace the existing source directory {arch}")
    if args.clean:
        if arch.is_symlink():
            arch.unlink()
        if build.exists():
            shutil.rmtree(build)
        print(f"Removed build products: {build}")
        return

    if not (source / "configure").is_file():
        raise FileNotFoundError(f"Missing bundled PETSc source: {source}")
    env = dict(os.environ, PETSC_DIR=str(source), PETSC_ARCH=args.arch)
    # PETSc controls its own package parallelism. Do not propagate top-level
    # command-line assignments into its Makefiles or reset the jobserver.
    env.pop("MAKEFLAGS", None)
    env.pop("MFLAGS", None)
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    logs = build / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    make = shutil.which("make")
    if not make:
        raise FileNotFoundError("GNU Make is required")
    make_args = [make, f"PETSC_DIR={source}", f"PETSC_ARCH={args.arch}"]
    metadata = json.loads((ROOT / "dependencies/sources.json").read_text())
    package = ROOT / metadata["hypre"]["archive"]
    if hashlib.sha256(package.read_bytes()).hexdigest() != metadata["hypre"]["sha256"]:
        raise ValueError(f"Bundled Hypre checksum mismatch: {package}")
    work.mkdir(parents=True, exist_ok=True)
    if not arch.is_symlink():
        arch.symlink_to(os.path.relpath(work, source), target_is_directory=True)
    configure = [sys.executable, str(source / "configure"), f"PETSC_ARCH={args.arch}",
                 f"--prefix={prefix}", "--with-debugging=0", "--with-shared-libraries=1",
                 "--with-clanguage=c", "--with-scalar-type=real", "--with-precision=double",
                 "--with-fortran-bindings=0", "--with-petsc4py=0",
                 f"--with-cc={args.cc}", f"--with-cxx={args.cxx}", f"--with-fc={args.fc}",
                 f"--with-make-np={args.jobs}", f"--download-hypre={package}",
                 f"COPTFLAGS={args.optflags}", f"CXXOPTFLAGS={args.optflags}",
                 f"FOPTFLAGS={args.optflags}"]
    run(configure, source, env, logs / "petsc-configure.log")
    run(make_args + [f"-j{args.jobs}", "all"], source, env, logs / "petsc-build.log")
    run(make_args + ["install-lib"], source, env, logs / "petsc-install.log")
    conf = (prefix / "include/petscconf.h").read_text()
    if "#define PETSC_HAVE_HYPRE 1" not in conf:
        raise RuntimeError("The built PETSc installation is missing Hypre")
    if not (prefix / "lib/pkgconfig/PETSc.pc").is_file():
        raise RuntimeError("PETSc install did not produce PETSc.pc")
    manifest = dict(sources=metadata, configure=configure, source=str(source),
                    work=str(work), prefix=str(prefix), arch=args.arch)
    (prefix / ".babelsim-built").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"PETSc/Hypre installed: {prefix}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as error:
        sys.exit(str(error))
