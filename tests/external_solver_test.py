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
    # 两个新头文件分别可用；math 单独求值不需要包含方程 API。
    (work / "math_api.cpp").write_text(
        '#include "babelsim/math.h"\nusing namespace babelsim;\n'
        'void gradient(const ScalarField& p, VectorField& result){ math::evaluate(math::grad(p),result); }\n')
    (work / "eqn_api.cpp").write_text(
        '#include "babelsim/eqn.h"\nusing namespace babelsim;\n'
        'auto heat(ScalarField& T){ return eqn::ddt(T) == eqn::laplacian(1.0,T) + eqn::source(2.0); }\n'
        'auto momentum(ScalarField& phi,VectorField& U,ScalarField& p){\n'
        ' return eqn::div(phi,U) == -math::grad(p) + eqn::laplacian(0.1,U); }\n')
    for name in ("math_api.cpp", "eqn_api.cpp"):
        run("g++", "-std=c++17", "-Wall", "-Wextra", "-Werror", "-Iinclude",
            "-fsyntax-only", name, cwd=work)
    # physics() 与 case.bs 的条目同名；旧接口不再作为并存别名暴露。
    for method, success in (("physics", True), ("properties", False)):
        (work / "case_api.cpp").write_text(
            '#include "babelsim/case.h"\n'
            f'double density(babelsim::Case& problem){{ return problem.{method}().positive("density"); }}\n')
        result = run("g++", "-std=c++17", "-Wall", "-Wextra", "-Werror", "-Iinclude",
                     "-fsyntax-only", "case_api.cpp", cwd=work, success=success)
        if not success:
            assert method in result.stderr
    for name in ("fvm", "fvc"):
        assert not (work / "include/babelsim" / (name + ".h")).exists()
        (work / "retired_api.cpp").write_text(
            '#include "babelsim/eqn.h"\nusing namespace babelsim;\n'
            f'void retired(ScalarField& T){{ {name}::laplacian(1.0,T); }}\n')
        failure = run("g++", "-std=c++17", "-Iinclude", "-fsyntax-only",
                      "retired_api.cpp", cwd=work, success=False)
        assert name in failure.stderr
    shutil.copy(ROOT / "tests/external/solver.cpp", work / "solver.cpp")
    run("g++", "-std=c++17", "-Wall", "-Wextra", "-Werror", "-Iinclude",
        "-c", "solver.cpp", "-o", "solver.o", cwd=work)
    run("mpic++", "solver.o", "libbabelsim.a", "-o", "external-solver", cwd=work)

    # 把整个 SIMPLE 模块（不是只有调用 SimpleSolver 的短 main）当成外部源码构建。
    # 不提供 src/、MPI/Eigen 头或框架私有接口；算法自己的 state.h 可正常使用。
    shutil.copytree(ROOT / "src/physics/simple", work / "simple")
    simple_objects = []
    for source in sorted((work / "simple").glob("*.cpp")):
        output = source.with_suffix(".o")
        run("g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", "-Iinclude",
            "-c", source, "-o", output, cwd=work)
        simple_objects.append(output)
    (work / "simple_entry.cpp").write_text(
        '#include "babelsim/application.h"\n'
        'int main(int argc,char** argv){ return babelsim::runApplication(argc,argv); }\n')
    run("g++", "-std=c++17", "-Iinclude", "-c", "simple_entry.cpp", "-o", "simple_entry.o", cwd=work)
    run("mpic++", "simple_entry.o", *simple_objects, "libbabelsim.a", "-o", "external-simple", cwd=work)
    simple_case = work / "simple_case"
    shutil.copytree(ROOT / "cases/poiseuille", simple_case, ignore=shutil.ignore_patterns("results", "post"))
    for count in (1, 2, 4):
        result = run("mpirun", "-np", count, work / "external-simple", "-case", simple_case,
                     "-time", f"np{count}", cwd=work)
        assert "converged=true" in result.stdout, result.stdout
        # 一次诊断由框架统一输出，不能让每个进程都打印一份相同日志。
        assert sum("converged=true" in line for line in result.stdout.splitlines()) == 1
    for count in (2, 4):
        run("python3", ROOT / "tools/compare_parallel_results.py", simple_case / "results/np1",
            simple_case / f"results/np{count}", "--atol", "5e-6", "--rtol", "5e-6", cwd=work)
    print("external_solver_test: complete SIMPLE module rebuilt with public headers; 1/2/4 ranks passed")
    (work / "single.cpp").write_text(
        '#include "babelsim/application.h"\nint solveCase(babelsim::Case&){return 0;}\n'
        'const babelsim::SolverRegistration single("single",solveCase);\n'
        'int main(int argc,char** argv){ return babelsim::runApplication(argc,argv); }\n')
    run("g++", "-std=c++17", "-Iinclude", "-fsyntax-only", "single.cpp", cwd=work)

    negative = {
        "field_data": 'auto pointer = field.data(); (void)pointer;',
        "field_index": 'field[0] = 1.0;',
        "field_resize": 'field.mutableData().resize(10);',
        "mesh_arrays": 'mesh.cell_volumes[0] = 0.0;',
        "mesh_partition": 'auto cells = mesh.ownedCellCount(); (void)cells;',
        "mesh_mutation": 'mesh.setOwnership({1,1,1}, 0, 0, 1, 0);',
        "mesh_replacement": 'mesh = Mesh::cartesian({2,1,1}, {}, {1,1,1});',
        "equation_storage": 'auto e = eqn::ddt(field) == 0.0; e.discrete();',
        "face_kernel": 'math::integratedNormalGradient(field, field, 0);',
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

    for solver, physics in (
        ("transport_extension", "diffusivity 0.1\nsource 2\n"),
        ("coupled_extension", "diffusivity 0.1\ncoupling 1\n"),
        ("vector_extension", "strength 2\n"),
    ):
        case = work / solver
        shutil.copytree(ROOT / "cases/heat", case, ignore=shutil.ignore_patterns("results", "post"))
        path = case / "case.bs"
        path.write_text(path.read_text().replace("solver heat", f"solver {solver}"))
        (case / "physics/thermal.bs").write_text(physics)
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
        ' (void)problem.physics().number("strength");\n'
        ' auto& field=problem.scalarField("failed",0.0); problem.output(field);\n'
        ' const char* rank=std::getenv("OMPI_COMM_WORLD_RANK");\n'
        ' return rank && rank[0]==\'0\' ? -1 : 0; }\n'
        '#if defined(TEST_EMPTY)\n'
        '#elif defined(TEST_NULL)\n'
        ' const babelsim::SolverRegistration entry(nullptr,failure);\n'
        '#elif defined(TEST_NULL_FUNCTION)\n'
        ' const babelsim::SolverRegistration entry("vector_extension",nullptr);\n'
        '#elif defined(TEST_EMPTY_NAME)\n'
        ' const babelsim::SolverRegistration entry("",failure);\n'
        '#elif defined(TEST_UNKNOWN)\n'
        ' const babelsim::SolverRegistration entry("different_solver",failure);\n'
        '#else\n'
        ' const babelsim::SolverRegistration entry("vector_extension",failure);\n'
        '#endif\n'
        'int main(int argc,char** argv){ return babelsim::runApplication(argc,argv); }\n')
    # 重名来自另一个源文件，确保错误检查不依赖翻译单元的静态初始化顺序。
    (work / "duplicate.cpp").write_text(
        '#include "babelsim/application.h"\nint failure(babelsim::Case&);\n'
        'const babelsim::SolverRegistration duplicate("vector_extension",failure);\n')
    run("g++", "-std=c++17", "-Iinclude", "-c", "duplicate.cpp", "-o", "duplicate.o", cwd=work)
    for name, definition, expected in (("negative", "TEST_NEGATIVE", None),
                                       ("duplicate", "TEST_DUPLICATE", "duplicate solver registration"),
                                       ("empty", "TEST_EMPTY", "no registered solvers"),
                                       ("null", "TEST_NULL", "invalid solver registration"),
                                       ("null_function", "TEST_NULL_FUNCTION", "invalid solver registration"),
                                       ("empty_name", "TEST_EMPTY_NAME", "empty solver name"),
                                       ("unknown", "TEST_UNKNOWN", "unknown BabelSim solver")):
        run("g++", "-std=c++17", "-Iinclude", f"-D{definition}", "-c", "failure.cpp", "-o", "failure.o", cwd=work)
        objects = ["failure.o", "duplicate.o"] if name == "duplicate" else ["failure.o"]
        run("mpic++", *objects, "libbabelsim.a", "-o", name, cwd=work)
        failure = run("mpirun", "-np", 2, work / name, "-case", case, "-time", name, cwd=work, success=False)
        if expected:
            assert expected in failure.stderr
        assert not (case / "results" / name).exists(), "failed application wrote successful output"

print("external_solver_test: out-of-tree equation/coupled/vector solvers, 1/2/4 ranks, 13 negative API checks, 7 application failures, MPI-free reader passed")
