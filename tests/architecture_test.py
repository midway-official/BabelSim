"""维护者约束：按真实头文件依赖检查层级，而不把目录名称当成隔离证明。"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
files = sorted(path for folder in (ROOT / "include", ROOT / "src")
               for path in folder.rglob("*") if path.suffix in (".h", ".cpp"))
texts = {path: path.read_text() for path in files}
edges = {}
for path, text in texts.items():
    dependencies = set()
    for delimiter, name in re.findall(r'#\s*include\s*([<"])([^">]+)[">]', text):
        candidates = [directory / name for directory in (path.parent, ROOT / "include", ROOT / "src")]
        target = next((candidate.resolve() for candidate in candidates if candidate.is_file()), None)
        if target is None and delimiter == '<' and not name.startswith(("babelsim/", "internal/")):
            continue  # 标准库/系统头不属于项目依赖图。
        assert target in texts, (path, "unresolved project include", name)
        dependencies.add(target)
    edges[path] = dependencies


def closure(path, stack=()):
    assert path not in stack, "include cycle: " + " -> ".join(str(p.relative_to(ROOT)) for p in (*stack, path))
    result = {path}
    for dependency in edges[path]:
        result.update(closure(dependency, (*stack, path)))
    return result


closures = {path: closure(path) for path in files}
# 公开语义只使用 eqn/math；保留内部 FvmExecution 这个有限体积后端名称。
for retired in ("include/babelsim/fvm.h", "include/babelsim/fvc.h",
                "src/discretization/fvm_expression.cpp"):
    assert not (ROOT / retired).exists(), retired
for path, text in texts.items():
    assert not re.search(r'\b(?:fvm|fvc)\s*::|namespace\s+(?:fvm|fvc)\b|'
                         r'babelsim/(?:fvm|fvc)\.h|\b(?:FvmTermKind|ScalarFvmTerm|VectorFvmTerm)\b', text), path
# 串行代数只认识 A/b/x 和求解配置；不能经装配头间接依赖 Mesh/Field/FVM。
for path in closures[ROOT / "include/babelsim/linear_solver.h"]:
    assert path.name in {"linear_solver.h", "solver_control.h"}, path
for name in ("mesh.h", "field.h"):
    assert all(path.name in {"mesh.h", "field.h", "vector.h"}
               for path in closures[ROOT / "include/babelsim" / name]), name

implementation_headers = {
    "runtime.h", "parallel.h", "mpi_support.h", "linear_solver.h", "assembly.h",
    "distributed_solver.h", "discrete_equation.h", "operators.h",
}
for name in ("case.h", "solver.h", "math.h", "eqn.h", "application.h", "postprocess.h",
             "result_reader.h", "simple.h", "simple_control.h"):
    for path in closures[ROOT / "include/babelsim" / name]:
        assert path.name not in implementation_headers, (name, path)
        assert not re.search(r'#include\s*[<"](?:mpi|Eigen)', texts[path]), path

# include/ 自身提供全部项目头依赖，不能偷偷要求使用者再添加 -Isrc。
for path in (ROOT / "include").rglob("*.h"):
    assert all(dependency.is_relative_to(ROOT / "include")
               for dependency in closures[path]), path

# 应用只声明 Solver 分派表；格式实现归 IO，MPI 生命周期归运行基础设施。
for path in (ROOT / "src/apps").glob("*.cpp"):
    for dependency in closures[path]:
        assert dependency.name not in implementation_headers, (path, dependency)
        assert not re.search(r'#include\s*[<"](?:mpi|Eigen)', texts[dependency]), dependency
    assert not re.search(r'MPI_|ParallelContext|HaloExchange|\.data\(|std::ofstream|'
                         r'face_vertices|owned_cells', texts[path]), path

# RunTime 负责运行生命周期与时间推进；FVM 的装配/工作区不能重新回到这个文件。
runtime = texts[ROOT / "src/runtime/runtime.cpp"]
assert not re.search(r'Eigen::|SparseAssembly|DiscreteEquation|EquationTerm|gradient_workspace', runtime)
assert "integratedNormalGradient" not in texts[ROOT / "include/babelsim/math.h"]
assert all(path.name != "runtime.h" for path in
           closures[ROOT / "src/discretization/fvm_execution.cpp"])

# 底层不能导入 Solver 或它的私有状态；不存在 SIMPLE 专用的跨层桥接。
for path in files:
    if path.is_relative_to(ROOT / "src") and path.relative_to(ROOT / "src").parts[0] in {
            "core", "discretization", "algebra", "runtime", "parallel", "io"}:
        for dependency in closures[path]:
            assert not dependency.is_relative_to(ROOT / "src/physics"), (path, dependency)
            assert not dependency.is_relative_to(ROOT / "src/apps"), (path, dependency)
            assert dependency.name not in {"simple.h", "simple_control.h"}, (path, dependency)
assert not (ROOT / "src/internal/simple_discretization.h").exists()
assert not (ROOT / "src/discretization/simple_discretization.cpp").exists()


def check_solver(path, text, dependencies):
    # 本算法内的私有状态合法，但所有跨模块能力必须来自公开 Solver API。
    module = ROOT / "src/physics" / path.relative_to(ROOT / "src/physics").parts[0]
    for dependency in dependencies:
        assert dependency.is_relative_to(ROOT / "include") or dependency.is_relative_to(module), (path, dependency)
        assert dependency.name not in implementation_headers, (path, dependency)
    code = re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"', '', text, flags=re.S)
    assert not re.search(r'\bdetail\s*::|MPI_|ParallelContext|HaloExchange|mutableData|'
                         r'\.(?:data|values|internal)\s*\(|std::vector|SparseAssembly|'
                         r'\b(?:CSR|LDU|Eigen|owned_cells|owned_faces|ghost|communicator)\b', code), path


for path in files:
    if path.is_relative_to(ROOT / "src/physics"):
        check_solver(path, texts[path], closures[path])

# 负向验收：检查器本身必须拒绝曾经漏过的私有接口、存储访问和运行后端包含链。
momentum = ROOT / "src/physics/simple/momentum.cpp"
for body, dependency in (
    ('detail::fieldData(m_U)[0] = Vec3{};', None),
    ('detail::meshData(m_U.mesh());', None),
    ('m_U.mutableData();', None),
    ('', ROOT / "src/internal/field_access.h"),
    ('', ROOT / "include/babelsim/runtime.h"),
    ('', ROOT / "include/babelsim/operators.h"),
):
    try:
        check_solver(momentum, texts[momentum] + body,
                     closures[momentum] | ({dependency} if dependency else set()))
    except AssertionError:
        continue
    raise AssertionError("architecture guard accepted an internal dependency: " + body + str(dependency))

# 主 Solver 不操作底层对象；数值算法的修正循环不属于存储循环。
for path in list((ROOT / "src/physics").glob("*/main.cpp")) + [
        ROOT / "tests/examples/coupled_scalar.cpp", ROOT / "tests/external/solver.cpp"]:
    assert not re.search(r'MPI_|ParallelContext|HaloExchange|mutableData|\.data\(|'
                         r'std::vector|unique_ptr|shared_ptr|RunTime|SparseAssembly', path.read_text()), path
control = texts[ROOT / "include/babelsim/simple_control.h"]
assert "LinearSolverConfig" not in control and "simpleRunTimeControl" not in control
for name in ("thermal.h", "transport.h", "equation.h", "solvers.h"):
    assert not (ROOT / "include/babelsim" / name).exists(), name

print(f"architecture_test: {len(files)} sources/headers, acyclic includes and layer boundaries passed")
