"""从新 Solver 作者的角度检查：公共 API、真实时间输出和串并行工作流。"""
import csv
import os
import re
from pathlib import Path
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
ENV = dict(os.environ, TMPDIR="/tmp", OMPI_ALLOW_RUN_AS_ROOT="1",
           OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")

# 普通 Solver 的头文件闭包不得重新带入执行/代数/并行实现。
seen = set()
def public_header(name):
    if name in seen:
        return
    seen.add(name)
    assert name not in {"runtime.h", "parallel.h", "linear_solver.h", "assembly.h",
                        "distributed_solver.h", "simple_control.h"}, name
    text = (ROOT / "include/babelsim" / name).read_text()
    assert not re.search(r"#include\s*[<\"](?:mpi|Eigen)", text), name
    for child in re.findall(r'#include "babelsim/([^"]+)"', text):
        public_header(child)

for name in ("case.h", "solver.h", "simple.h"):
    public_header(name)
for path in list((ROOT / "src/physics").glob("*/main.cpp")) + [
        ROOT / "tests/examples/coupled_scalar.cpp"]:
    assert not re.search(r"MPI_|ParallelContext|HaloExchange|mutableData|\.data\(|"
                         r"std::vector|unique_ptr|shared_ptr|RunTime|SparseAssembly", path.read_text()), path


def run(*args, success=True):
    process = subprocess.run([str(a) for a in args], cwd=ROOT, env=ENV,
                             text=True, capture_output=True, timeout=90)
    if (process.returncode == 0) != success:
        raise AssertionError(f"{args}: exit={process.returncode}\n{process.stdout}\n{process.stderr}")
    return process


def clone_case(base, name, source):
    target = base / name
    shutil.copytree(ROOT / "cases" / source, target,
                    ignore=shutil.ignore_patterns("results", "post"))
    return target


def values(directory, field):
    result = {}
    for file in sorted(directory.glob(f"rank-*/{field}.csv")):
        with file.open() as stream:
            for row in csv.DictReader(stream):
                index = int(row["global_id"])
                assert index not in result
                result[index] = tuple(float(row[key]) for key in row if key.startswith("value"))
    assert result, directory
    return result


def compare(left, right, field, tolerance=1e-8):
    a, b = values(left, field), values(right, field)
    assert a.keys() == b.keys()
    difference = max(abs(x-y) for cell in a for x, y in zip(a[cell], b[cell]))
    assert difference < tolerance, (left, right, field, difference)
    return difference


def sequence(case, name, expected):
    directory = case / "results" / name
    actual = sorted(float(p.name) for p in directory.iterdir()
                    if p.is_dir() and not p.name.startswith("rank-"))
    assert actual == expected, (actual, expected)
    for time in expected:
        text = (directory / f"{time:g}" / "rank-0000" / "metadata.bs").read_text()
        assert f"time {time:g}\n" in text
    return directory


def check_pvd(case, name, expected):
    run(ROOT / "build/babelsim-post", "-case", case, "-time", f"{name}/all",
        "-format", "vtk", "tecplot")
    path = case / "post" / name / "series.pvd"
    datasets = ET.parse(path).findall("./Collection/DataSet")
    assert [float(item.attrib["timestep"]) for item in datasets] == expected
    for item in datasets:
        file = path.parent / item.attrib["file"]
        assert file.suffix == ".vtu" and file.exists()
        piece = ET.parse(file).find("./UnstructuredGrid/Piece")
        count = int(piece.attrib["NumberOfCells"])
        for array in piece.findall("./CellData/DataArray"):
            assert len(array.text.split()) == count * int(array.attrib["NumberOfComponents"])
    return path


with tempfile.TemporaryDirectory(prefix="babelsim-workflow-") as temporary:
    base = Path(temporary)
    heat = clone_case(base, "heat", "heat")
    (heat / "physics/thermal.bs").write_text("density 1\nheatCapacity 1\nconductivity 0\nsource 2\n")
    (heat / "output.bs").write_text("directory results\ntimeName final\nwriteInterval 2\n")
    for count in (1, 2, 4):
        name = f"np{count}"
        run("mpirun", "-np", count, ROOT / "build/babelsim-solve", "-case", heat, "-time", name)
        directory = sequence(heat, name, [0.02, 0.04, 0.05])
        for time in (0.02, 0.04, 0.05):
            assert all(abs(v[0] - 2*time) < 1e-11 for v in values(directory / f"{time:g}", "T").values())
            if count > 1:
                compare(heat / "results/np1" / f"{time:g}", directory / f"{time:g}", "T")
    check_pvd(heat, "np4", [0.02, 0.04, 0.05])

    transport = clone_case(base, "transport", "transport")
    for count in (1, 2, 4):
        run("mpirun", "-np", count, ROOT / "build/babelsim-solve",
            "-case", transport, "-time", f"np{count}")
        sequence(transport, f"np{count}", [0.01, 0.02, 0.03, 0.04, 0.05])
        if count > 1:
            for time in (0.01, 0.03, 0.05):
                compare(transport / "results/np1" / f"{time:g}",
                        transport / "results" / f"np{count}" / f"{time:g}", "C")

    coupled = clone_case(base, "coupled", "heat")
    (coupled / "control.bs").write_text("startTime 0\nendTime 0.25\ndeltaT 0.1\n")
    (coupled / "physics/thermal.bs").write_text("diffusivity 0.1\ncoupling 1\n")
    (coupled / "numerics/solution.bs").write_text(
        "scalarSolver bicgstab ilut 1e-14 1e-12 1000\n"
        "couplingIterations 100\ncouplingTolerance 1e-12\n")
    field = (coupled / "fields/initial/T.field").read_text()
    field = field.replace("type fixedValue value (1)", "type zeroGradient")
    field = field.replace("type fixedValue value (0)", "type zeroGradient")
    (coupled / "fields/initial/T.field").write_text(field.replace("uniform (0)", "uniform (1)"))
    (coupled / "fields/initial/C.field").write_text(field.replace("field T", "field C"))
    tensor = field.replace("field T", "field stress").replace("type scalar", "type tensor")
    tensor = tensor.replace("uniform (0)", "uniform (1 2 3 4 5 6 7 8 9)")
    (coupled / "fields/initial/stress.field").write_text(tensor)
    T, C = 1.0, 0.0
    expected = {}
    for time, dt in ((0.1, 0.1), (0.2, 0.1), (0.25, 0.05)):
        T, C = (T + dt*C)/(1-dt*dt), (C + dt*T)/(1-dt*dt)
        expected[time] = (T, C)
    for count in (1, 2, 4):
        run("mpirun", "-np", count, ROOT / "build/case_programming_test", coupled, f"np{count}")
        directory = sequence(coupled, f"np{count}", [0.1, 0.2, 0.25])
        assert not list(directory.rglob("previous.csv")), "temporary mathematical field was written"
        assert not list(directory.rglob("scratch*.csv")), "temporary vector field was written"
        assert all(value == (1, 2, 3, 4, 5, 6, 7, 8, 9)
                   for value in values(directory / "0.25", "stress").values())
        for time, exact in expected.items():
            for field_name, answer in zip(("T", "C"), exact):
                error = max(abs(v[0]-answer)
                            for v in values(directory / f"{time:g}", field_name).values())
                assert error < 1e-9, (count, time, field_name, answer, error)
    check_pvd(coupled, "np4", [0.1, 0.2, 0.25])
    conflict = run("mpirun", "-np", 2, ROOT / "build/case_programming_test",
                   coupled, "np4", success=False)
    assert "incompatible partitions" in conflict.stderr

    # 数值健康/内迭代失败不能写出一个被标成完成的最终步。
    path = coupled / "numerics/solution.bs"
    path.write_text(path.read_text().replace("couplingIterations 100", "couplingIterations 1"))
    failed = run("mpirun", "-np", 2, ROOT / "build/case_programming_test",
                 coupled, "failed", success=False)
    assert not (coupled / "results/failed").exists()

    # 数值时间排序不是字典序排序；ParaView 读取的必须是物理时间，不是序号。
    (heat / "control.bs").write_text("startTime 0\nendTime 10\ndeltaT 2\n")
    (heat / "output.bs").write_text("directory results\ntimeName final\nwriteInterval 1\n")
    run(ROOT / "build/babelsim-solve", "-case", heat, "-time", "ordered")
    pvd = check_pvd(heat, "ordered", [2.0, 4.0, 6.0, 8.0, 10.0])
    run(ROOT / "build/babelsim-post", "-case", heat, "-time", "ordered/latest", "-format", "vtk")
    assert (pvd.parent / "10.vtu").exists()

    # 有 ParaView 时检查真实读取器，不以 XML 语法正确冒充可视化软件可用。
    if shutil.which("pvpython"):
        check = base / "read_pvd.py"
        check.write_text(
            "import site\nsite.addsitedir('/usr/lib/python3/dist-packages')\n"
            "from paraview.simple import PVDReader\n"
            f"reader = PVDReader(FileName={str(pvd)!r})\n"
            "assert list(reader.TimestepValues) == [2, 4, 6, 8, 10]\n"
            "reader.UpdatePipeline(time=10)\n"
            "assert reader.GetDataInformation().GetNumberOfCells() == 32\n"
            "assert reader.CellData.GetArray('T').GetRange()[0] > 19.999999\n")
        run("pvpython", check)
        print("ParaView: PVD time steps and VTU cell values passed")

    (heat / "output.bs").write_text("directory results\ntimeName final\nwriteInterval 0\n")
    run("mpirun", "-np", 2, ROOT / "build/babelsim-solve", "-case", heat, success=False)
    (heat / "output.bs").write_text("directory results\ntimeName final\n")
    with (heat / "physics/thermal.bs").open("a") as stream:
        stream.write("misspelledProperty 1\n")
    error = run("mpirun", "-np", 2, ROOT / "build/babelsim-solve", "-case", heat, success=False)
    assert "unused or unknown entry" in error.stderr

print("solver_workflow_test: heat/transport/coupled, 1/2/4 ranks, physical times, failure paths passed")
