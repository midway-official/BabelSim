#!/usr/bin/env python3
"""Reproduce review observations; success of this runner does not mean defects pass.

Prerequisite: make -j4 (same native build used in the review).
All generated executables and case outputs go into a fresh temporary directory.
"""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
WORK = Path(tempfile.mkdtemp(prefix="babelsim-review-reproduce-"))
ENV = dict(os.environ, TMPDIR="/tmp")


def run(command, name):
    result = subprocess.run(command, cwd=ROOT, env=ENV, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (WORK / (name + ".log")).write_text(result.stdout)
    print(name, "exit=", result.returncode)
    print(result.stdout[-1800:])
    result.check_returncode()


for name in ["rans", "sa", "post"]:
    shutil.copytree(HERE / (name + "-case"), WORK / (name + "-case"))

for name in ["scale", "rans", "sa-boundary", "gg", "gg3", "ls",
             "inletoutlet", "transient"]:
    source = (HERE / (name + ".cpp")).read_text()
    for case in ["rans", "sa"]:
        source = source.replace("/tmp/babelsim-review-" + case + "-case",
                                str(WORK / (case + "-case")))
    target = WORK / (name + ".cpp")
    target.write_text(source)
    exe = WORK / name
    run(["mpic++", "-std=c++17", "-O2", "-march=native", "-fno-lto",
         "-Iinclude", "-Isrc", "-Itests", "-I/usr/include/eigen3", str(target),
         "build/libbabelsim.a", "-o", str(exe)], name + "-build")
    command = (["mpirun", "-np", "4", str(exe)]
               if name in ["gg", "gg3", "ls"] else [str(exe)])
    run(command, name)

run([str(ROOT / "build/babelsim-post"), "-case", str(WORK / "post-case"),
     "-time", "original", "-format", "vtk"], "post")
print("Evidence:", WORK)
