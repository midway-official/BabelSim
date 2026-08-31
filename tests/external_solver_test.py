"""真正仓库外构建：只用公开头和预编译库，同时验收禁止访问的实现边界。"""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ENV = dict(os.environ, TMPDIR="/tmp", OMPI_ALLOW_RUN_AS_ROOT="1",
           OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")


def run(*args, cwd, success=True):
    result = subprocess.run([str(arg) for arg in args], cwd=cwd, env=ENV,
                            text=True, capture_output=True, timeout=120)
    assert (result.returncode == 0) == success, (args, result.stdout, result.stderr)
    return result


def rows(directory, name):
    values = {}
    for path in directory.glob(f"rank-*/{name}.csv"):
        with path.open() as stream:
            for row in csv.DictReader(stream):
                cell = int(row["global_id"])
                assert cell not in values
                values[cell] = row
    assert len(values) == 32, (directory, name, len(values))
    return values


with tempfile.TemporaryDirectory(prefix="babelsim-external-") as temporary:
    work = Path(temporary)
    # staged SDK 不含 src/、Eigen 或 MPI 头；Solver 编译阶段只需标准 C++ 编译器。
    shutil.copytree(ROOT / "include", work / "include")
    shutil.copy(ROOT / "build/libbabelsim.a", work / "libbabelsim.a")
    shutil.copy(ROOT / "tests/external/solver.cpp", work / "solver.cpp")
    run("g++", "-std=c++17", "-Wall", "-Wextra", "-Werror", "-Iinclude",
        "-c", "solver.cpp", "-o", "solver.o", cwd=work)
    run("mpic++", "solver.o", "libbabelsim.a", "-o", "external-solver", cwd=work)
    (work / "single.cpp").write_text(
        '#include "babelsim/application.h"\nint solveCase(babelsim::Case&){return 0;}\n'
        'int main(int argc,char** argv){ return babelsim::runApplication(argc,argv,{"single",solveCase}); }\n')
    run("g++", "-std=c++17", "-Iinclude", "-fsyntax-only", "single.cpp", cwd=work)

    negative = {
        "field_data": 'auto pointer = field.data(); (void)pointer;',
        "field_index": 'field[0] = 1.0;',
        "field_resize": 'field.mutableData().resize(10);',
        "mesh_arrays": 'mesh.cell_volumes[0] = 0.0;',
        "mesh_partition": 'auto cells = mesh.ownedCellCount(); (void)cells;',
        "mesh_mutation": 'mesh.setOwnership({1,1,1}, 0, 0, 1, 0);',
        "mesh_replacement": 'mesh = Mesh::cartesian({2,1,1}, {}, {1,1,1});',
        "equation_storage": 'auto e = fvm::ddt(field) == 0.0; e.discrete();',
        "face_kernel": 'fvc::integratedNormalGradient(field, field, 0);',
    }
    expected_diagnostics = {
        "field_data": "data", "field_index": "operator[]", "field_resize": "mutableData",
        "mesh_arrays": "cell_volumes", "mesh_partition": "ownedCellCount",
        "mesh_mutation": "setOwnership", "mesh_replacement": "operator=",
        "equation_storage": "discrete", "face_kernel": "integratedNormalGradient",
    }
    for name, body in negative.items():
        source = work / f"{name}.cpp"
        source.write_text('#include "babelsim/solver.h"\nusing namespace babelsim;\n'
                          'int main(){ Mesh mesh = Mesh::cartesian({1,1,1}, {}, {1,1,1});\n'
                          'ScalarField field(mesh, FieldLocation::Cell);\n' + body + '\n}\n')
        failure = run("g++", "-std=c++17", "-Iinclude", "-fsyntax-only", source, cwd=work, success=False)
        assert expected_diagnostics[name] in failure.stderr, (name, failure.stderr)
    (work / "private.cpp").write_text('#include "internal/field_access.h"\n')
    run("g++", "-Iinclude", "-fsyntax-only", "private.cpp", cwd=work, success=False)

    # 后处理读取器以普通 g++ 链接且不带 libmpi，防止头文件/静态库传递依赖回归。
    (work / "reader.cpp").write_text(
        '#include "babelsim/result_reader.h"\n'
        'int main(int argc,char** argv){ if(argc!=2)return 2;\n'
        'return babelsim::readParallelResults(argv[1], 32).fields.empty(); }\n')
    run("g++", "-std=c++17", "-Iinclude", "reader.cpp", "libbabelsim.a", "-o", "reader", cwd=work)

    for solver, properties in (
        ("transport_extension", "diffusivity 0.1\nsource 2\n"),
        ("coupled_extension", "diffusivity 0.1\ncoupling 1\n"),
        ("vector_extension", "strength 2\n"),
    ):
        case = work / solver
        shutil.copytree(ROOT / "cases/heat", case, ignore=shutil.ignore_patterns("results", "post"))
        path = case / "case.bs"
        path.write_text(path.read_text().replace("solver heat", f"solver {solver}"))
        (case / "physics/thermal.bs").write_text(properties)
        field = (case / "fields/initial/T.field").read_text().replace("field T", "field C")
        field = field.replace("type fixedValue value (1)", "type zeroGradient")
        field = field.replace("type fixedValue value (0)", "type zeroGradient")
        (case / "fields/initial/C.field").write_text(field)
        for count in (1, 2, 4):
            run("mpirun", "-np", count, work / "external-solver", "-case", case,
                "-time", f"np{count}", cwd=work)
            result = case / "results" / f"np{count}" / "0.05"
            run(work / "reader", result, cwd=work)
            if solver == "transport_extension":
                assert all(abs(float(row["value0"]) - 0.1) < 1e-10 for row in rows(result, "C").values())
            elif solver == "coupled_extension":
                T, C, dt = 1.0, 0.0, 0.01
                for _ in range(5):
                    T, C = (T + dt*C)/(1-dt*dt), (C + dt*T)/(1-dt*dt)
                for name, exact in (("T", T), ("C", C)):
                    assert all(abs(float(row["value0"]) - exact) < 1e-9 for row in rows(result, name).values())
                assert not list(result.rglob("previous*.csv"))
            else:
                for row in rows(result, "U").values():
                    for component, coordinate in enumerate(("x", "y", "z")):
                        exact = 0.1 * (component + 1 + float(row[coordinate]))
                        assert abs(float(row[f"value{component}"]) - exact) < 1e-9
                for name, exact in (("p", 3.0), ("rAU", 0.01)):
                    assert all(abs(float(row["value0"]) - exact) < 1e-9 for row in rows(result, name).values())
                for row in rows(result, "energy").values():
                    exact = 0.005*sum((i + 1 + float(row[axis]))**2 for i, axis in enumerate(("x", "y", "z")))
                    assert abs(float(row["value0"]) - exact) < 1e-9
                for row in rows(result, "stress").values():
                    assert [float(row[f"value{i}"]) for i in range(9)] == list(range(1, 10))
                assert not list(result.rglob("force.csv"))

    # 这是启动器的负向夹具，故意模拟不一致返回码；不是普通 Solver 编程示例。
    (work / "failure.cpp").write_text(
        '#include "babelsim/application.h"\n#include "babelsim/case.h"\n#include <cstdlib>\n'
        'int failure(babelsim::Case& problem){\n'
        ' (void)problem.properties().number("strength");\n'
        ' auto& field=problem.scalarField("failed",0.0); problem.output(field);\n'
        ' const char* rank=std::getenv("OMPI_COMM_WORLD_RANK");\n'
        ' return rank && rank[0]==\'0\' ? -1 : 0; }\n'
        'int main(int argc,char** argv){\n'
        '#if defined(TEST_DUPLICATE)\n'
        ' const babelsim::SolverEntry entries[]={{"vector_extension",failure},{"vector_extension",failure}};\n'
        ' return babelsim::runApplication(argc,argv,entries);\n'
        '#elif defined(TEST_NULL)\n'
        ' return babelsim::runApplication(argc,argv,nullptr,1);\n'
        '#else\n'
        ' return babelsim::runApplication(argc,argv,{"vector_extension",failure});\n'
        '#endif\n}\n')
    for name, definition, expected in (("negative", "TEST_NEGATIVE", None),
                                       ("duplicate", "TEST_DUPLICATE", "duplicate solver entry"),
                                       ("null", "TEST_NULL", "empty solver table")):
        run("g++", "-std=c++17", "-Iinclude", f"-D{definition}", "-c", "failure.cpp", "-o", "failure.o", cwd=work)
        run("mpic++", "failure.o", "libbabelsim.a", "-o", name, cwd=work)
        failure = run("mpirun", "-np", 2, work / name, "-case", case, "-time", name, cwd=work, success=False)
        if expected:
            assert expected in failure.stderr
        assert not (case / "results" / name).exists(), "failed application wrote successful output"

print("external_solver_test: out-of-tree equation/coupled/vector solvers, 1/2/4 ranks, 10 negative API checks, 3 application failures, MPI-free reader passed")
