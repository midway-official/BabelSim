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
    for name in re.findall(r'#include\s*"([^"]+)"', text):
        candidates = [directory / name for directory in (path.parent, ROOT / "include", ROOT / "src")]
        target = next((candidate.resolve() for candidate in candidates if candidate.is_file()), None)
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
for name in ("case.h", "solver.h", "fvc.h", "fvm.h", "application.h", "postprocess.h",
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
assert not re.search(r'Eigen::|SparseAssembly|DiscreteEquation|FvmTerm|gradient_workspace', runtime)
assert "integratedNormalGradient" not in texts[ROOT / "include/babelsim/fvc.h"]
assert all(path.name != "runtime.h" for path in
           closures[ROOT / "src/discretization/fvm_execution.cpp"])

# 专用数值桥接可调用通用执行层，但底层不能导入 Solver 或它的私有状态。
for path in files:
    if path.is_relative_to(ROOT / "src") and path.relative_to(ROOT / "src").parts[0] in {
            "core", "discretization", "algebra", "runtime", "parallel", "io"}:
        for dependency in closures[path]:
            assert not dependency.is_relative_to(ROOT / "src/physics"), (path, dependency)
            assert not dependency.is_relative_to(ROOT / "src/apps"), (path, dependency)
            assert dependency.name not in {"simple.h", "simple_control.h"}, (path, dependency)
for path in (ROOT / "src/discretization/operators.cpp", ROOT / "include/babelsim/operators.h"):
    assert all(dependency.name != "simple_discretization.h" for dependency in closures[path]), path
    assert "applyMomentumInterpolation" not in texts[path], path

# 算法方程与主步骤不需要执行器定义；仅初始化和日志实现使用执行上下文。
for name in ("momentum.cpp", "pressure.cpp", "simple_solver.cpp", "state.h"):
    for dependency in closures[ROOT / "src/physics/simple" / name]:
        assert dependency.name not in implementation_headers, (name, dependency)
for path in files:
    if path.is_relative_to(ROOT / "src/physics"):
        assert not re.search(r'MPI_|ParallelContext|HaloExchange|mutableData|\.data\(|'
                             r'std::vector|SparseAssembly|Eigen::', texts[path]), path

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
