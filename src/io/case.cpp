#include "babelsim/case.h"

#include "babelsim/field_io.h"
#include "babelsim/mpi_support.h"
#include "babelsim/numerics_io.h"
#include "babelsim/parallel_writer.h"
#include "babelsim/runtime.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace babelsim {
namespace {

ParallelContext activeParallel() {
    int initialized = 0;
    detail::checkMpi(MPI_Initialized(&initialized), "MPI_Initialized");
    if (!initialized) return {};
    int finalized = 0;
    detail::checkMpi(MPI_Finalized(&finalized), "MPI_Finalized");
    if (finalized) throw std::logic_error("cannot create Case after MPI_Finalize");
    return ParallelContext::world();
}

RuntimeControl runtimeControl(const CaseDefinition& definition, const Parameters& solution) {
    RuntimeControl result;
    result.methods = readMethodsFile(definition.methods_file);
    result.time = readTimeControlFile(definition.control_file);
    // Case 必须显式选择两类线性求解器；缺项由 Parameters 给出文件路径和键名。
    readLinearSolverLine(definition.solution_file, solution.entry("scalarSolver"), result.scalar_solver);
    readLinearSolverLine(definition.solution_file, solution.entry("vectorSolver"), result.vector_solver);
    result.validate();
    return result;
}

void requireRelativePath(const std::filesystem::path& path) {
    if (path.empty() || path.is_absolute()) throw std::invalid_argument("output path must be relative");
    for (const auto& part : path) {
        if (part == ".." || part == ".") throw std::invalid_argument("output path cannot contain . or ..");
    }
}

void requireName(const std::string& name) {
    requireRelativePath(name);
    if (std::filesystem::path(name).has_parent_path())
        throw std::invalid_argument("expected a single field/run name: " + name);
}

void requireRunName(const std::string& name) {
    requireName(name);
    std::size_t consumed = 0;
    try { std::stod(name, &consumed); } catch (const std::exception&) {}
    if (consumed == name.size() || name == "all" || name == "latest")
        throw std::invalid_argument("run/final label must not be a time or reserved selection: " + name);
}

std::string timeName(double time) {
    std::ostringstream text;
    text << std::setprecision(15) << time;
    return text.str();
}

}  // 匿名命名空间

struct Case::Implementation {
    Implementation(const std::filesystem::path& directory, const std::string& run_name)
        : definition(readCase(directory)), output(readOutputControl(definition)),
          physics(definition.physics_file), solution(definition.solution_file),
          parallel(activeParallel()), mesh(readDistributedMesh(definition.mesh_file, parallel)),
          run_time(RunTime::forMesh(mesh, runtimeControl(definition, solution)))
    {
        requireRelativePath(output.directory);
        requireRunName(output.time_name);
        if (!run_name.empty()) requireRunName(run_name);
        const auto base = definition.root / output.directory;
        // 命名运行隔离时间序列，避免不同进程数/参数的实验互相覆盖。
        series_directory = run_name.empty() ? base : base / run_name;
        final_directory = base / (run_name.empty() ? output.time_name : run_name);
    }

    template <typename T>
    Field<T>& field(std::vector<std::unique_ptr<Field<T>>>& fields,
                    const std::string& name, FieldLocation location,
                    bool read_file = true, T initial = T{}) {
        requireName(name);
        for (const auto& value : fields) {
            if (value->name() == name && value->location() == location) return *value;
        }
        for (const auto& info : field_names) {
            if (info == name) throw std::invalid_argument("field name reused with another type/location: " + name);
        }
        if (started) throw std::logic_error("create fields before starting the time/algorithm loop");
        auto value = std::make_unique<Field<T>>(mesh, location, name, initial);
        if (location == FieldLocation::Cell && read_file) {
            readFieldFile(definition.fields_directory / (name + ".field"), *value);
            output_names.push_back(name);
        }
        fields.push_back(std::move(value));
        field_names.push_back(name);
        return *fields.back();
    }

    template <typename T>
    void writeFields(const std::filesystem::path& directory,
                     const std::vector<std::unique_ptr<Field<T>>>& fields,
                     const char* type, std::vector<FieldOutputInfo>& info) {
        for (const auto& field : fields) {
            if (std::find(output_names.begin(), output_names.end(), field->name()) ==
                output_names.end()) continue;
            writeOwnedFieldCsv(directory, *field, parallel);
            info.push_back({field->name(), type, field->location()});
        }
    }

    void write(const std::filesystem::path& directory) {
        // 不删除旧实验，也不把不同分区数量写进同一结果集。预检结果必须全局一致，
        // 否则某个进程抛异常、其他进程进入 writer 的 collective 会造成死锁。
        int incompatible = 0;
        if (parallel.rank == 0) {
            try {
                int ranks = 0;
                if (std::filesystem::exists(directory)) {
                    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
                        if (entry.is_directory() &&
                            entry.path().filename().string().rfind("rank-", 0) == 0) ++ranks;
                    }
                }
                if (ranks != 0 && ranks != parallel.size) incompatible = 1;
            } catch (const std::exception&) { incompatible = 1; }
        }
        if (parallel.maximum(incompatible))
            throw std::runtime_error("existing result has incompatible partitions; use a new run label: " +
                                     directory.string());
        std::vector<FieldOutputInfo> info;
        writeFields(directory, scalars, "scalar", info);
        writeFields(directory, vectors, "vector", info);
        writeFields(directory, tensors, "tensor", info);
        writeOwnedResultMetadata(directory, mesh, parallel, timeName(run_time.time()), info);
    }

    void writeStep(bool force) {
        if (run_time.step() == last_written_step) return;
        if (!force && run_time.step() % output.write_interval != 0) return;
        write(series_directory / timeName(run_time.time()));
        last_written_step = run_time.step();
    }

    CaseDefinition definition;
    OutputControl output;
    Parameters physics;
    Parameters solution;
    ParallelContext parallel;
    Mesh mesh;
    std::vector<std::unique_ptr<ScalarField>> scalars;
    std::vector<std::unique_ptr<VectorField>> vectors;
    std::vector<std::unique_ptr<TensorField>> tensors;
    std::vector<std::string> field_names;
    std::vector<std::string> output_names;
    std::filesystem::path series_directory;
    std::filesystem::path final_directory;
    // 后构造先析构：执行对象/历史场先于源场，源场先于 Mesh 销毁。
    RunTime run_time;
    int last_written_step = -1;
    bool started = false;
    bool finished = false;
};

Case::Case(const std::filesystem::path& directory, const std::string& run_name)
    : m_implementation(std::make_unique<Implementation>(directory, run_name)) {}
Case::~Case() = default;

const std::string& Case::solver() const { return m_implementation->definition.solver; }
const Mesh& Case::mesh() const { return m_implementation->mesh; }
const Parameters& Case::physics() const { return m_implementation->physics; }
const Parameters& Case::solution() const { return m_implementation->solution; }
double Case::time() const { return m_implementation->run_time.time(); }
int Case::step() const { return m_implementation->run_time.step(); }

ScalarField& Case::scalarField(const std::string& name) {
    return m_implementation->field(m_implementation->scalars, name, FieldLocation::Cell);
}
VectorField& Case::vectorField(const std::string& name) {
    return m_implementation->field(m_implementation->vectors, name, FieldLocation::Cell);
}
ScalarField& Case::scalarField(const std::string& name, double initial) {
    return m_implementation->field(m_implementation->scalars, name, FieldLocation::Cell, false, initial);
}
VectorField& Case::vectorField(const std::string& name, Vec3 initial) {
    return m_implementation->field(m_implementation->vectors, name, FieldLocation::Cell, false, initial);
}
TensorField& Case::tensorField(const std::string& name) {
    return m_implementation->field(m_implementation->tensors, name, FieldLocation::Cell);
}
TensorField& Case::tensorField(const std::string& name, Tensor3 initial) {
    return m_implementation->field(m_implementation->tensors, name, FieldLocation::Cell, false, initial);
}
ScalarField& Case::faceField(const std::string& name) {
    return m_implementation->field(m_implementation->scalars, name, FieldLocation::Face);
}
VectorField& Case::faceVectorField(const std::string& name) {
    return m_implementation->field(m_implementation->vectors, name, FieldLocation::Face);
}
TensorField& Case::faceTensorField(const std::string& name) {
    return m_implementation->field(m_implementation->tensors, name, FieldLocation::Face);
}
ScalarField& Case::faceFlux(const std::string& name, const VectorField& velocity) {
    ScalarField& result = faceField(name);
    math::evaluate(math::flux(velocity), result);
    return result;
}

void Case::selectOutput(const std::string& name, const void* field, bool enabled) {
    Implementation& state = *m_implementation;
    bool owned_cell_field = false;
    const auto check = [&](const auto& fields) {
        for (const auto& value : fields)
            if (value.get() == field && value->location() == FieldLocation::Cell)
                owned_cell_field = true;
    };
    check(state.scalars);
    check(state.vectors);
    check(state.tensors);
    if (!owned_cell_field) throw std::invalid_argument("output selection requires a Case-owned cell field");
    auto found = std::find(state.output_names.begin(), state.output_names.end(), name);
    if (enabled && found == state.output_names.end()) state.output_names.push_back(name);
    if (!enabled && found != state.output_names.end()) state.output_names.erase(found);
}

void Case::validate() const {
    physics().requireAllUsed();
    solution().requireAllUsed();
}

void Case::start() {
    if (m_implementation->started) return;
    validate();
    m_implementation->started = true;
}

bool Case::loop() {
    Implementation& state = *m_implementation;
    if (state.finished) return false;
    if (state.run_time.methods().time == TimeMethod::Steady)
        throw std::logic_error("steady case requires an algorithm iteration loop");
    start();
    if (step() > 0) state.writeStep(false);
    if (state.run_time.loop()) return true;
    finish();
    return false;
}

void Case::finish() {
    Implementation& state = *m_implementation;
    if (state.finished) return;
    start();
    state.writeStep(true);
    if (state.final_directory != state.series_directory / timeName(time()))
        state.write(state.final_directory);
    state.finished = true;
    if (state.parallel.rank == 0)
        std::cout << "BabelSim result time=" << time()
                  << " steps=" << step() << " saved to " << state.series_directory << '\n';
}

}  // babelsim 命名空间
