#include "babelsim/case.h"

#include "babelsim/field_io.h"
#include "babelsim/mpi_support.h"
#include "babelsim/numerics_io.h"
#include "babelsim/parallel_writer.h"
#include "babelsim/runtime.h"

#include <array>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

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

bool listed(const std::vector<std::string>& names, const std::string& value) {
    return std::find(names.begin(), names.end(), value) != names.end();
}

RuntimeControl runtimeControl(const CaseDefinition& definition, const Methods& methods) {
    RuntimeControl result;
    result.methods = methods;
    result.time = readTimeControlFile(definition.control_file);
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

PerformanceCounters maximumPerformance(
    const PerformanceCounters& local,
    const ParallelContext& parallel)
{
    const std::array<double, 19> local_values{
        static_cast<double>(local.linear_solves),
        static_cast<double>(local.krylov_iterations),
        static_cast<double>(local.sparse_matvecs),
        static_cast<double>(local.halo_exchanges),
        static_cast<double>(local.halo_bytes),
        static_cast<double>(local.global_reductions),
        static_cast<double>(local.equation_assemblies),
        static_cast<double>(local.preconditioner_setups),
        static_cast<double>(local.preconditioner_applications),
        static_cast<double>(local.output_writes),
        local.elapsed_seconds,
        local.assembly_seconds,
        local.preconditioner_seconds,
        local.preconditioner_apply_seconds,
        local.linear_solve_seconds,
        local.sparse_matvec_seconds,
        local.halo_seconds,
        local.global_reduction_seconds,
        local.output_seconds,
    };
    std::array<double, local_values.size()> global_values{};
    parallel.maximum(
        local_values.data(), global_values.data(),
        detail::mpiCount(global_values.size(), "performance counters"));
    PerformanceCounters global;
    global.linear_solves = static_cast<std::uint64_t>(global_values[0]);
    global.krylov_iterations = static_cast<std::uint64_t>(global_values[1]);
    global.sparse_matvecs = static_cast<std::uint64_t>(global_values[2]);
    global.halo_exchanges = static_cast<std::uint64_t>(global_values[3]);
    global.halo_bytes = static_cast<std::uint64_t>(global_values[4]);
    global.global_reductions = static_cast<std::uint64_t>(global_values[5]);
    global.equation_assemblies = static_cast<std::uint64_t>(global_values[6]);
    global.preconditioner_setups = static_cast<std::uint64_t>(global_values[7]);
    global.preconditioner_applications = static_cast<std::uint64_t>(global_values[8]);
    global.output_writes = static_cast<std::uint64_t>(global_values[9]);
    global.elapsed_seconds = global_values[10];
    global.assembly_seconds = global_values[11];
    global.preconditioner_seconds = global_values[12];
    global.preconditioner_apply_seconds = global_values[13];
    global.linear_solve_seconds = global_values[14];
    global.sparse_matvec_seconds = global_values[15];
    global.halo_seconds = global_values[16];
    global.global_reduction_seconds = global_values[17];
    global.output_seconds = global_values[18];
    return global;
}

}  // 匿名命名空间

struct Case::Implementation {
    mutable std::map<std::string, std::pair<const void*, bool>> equation_bindings;
    Implementation(const std::filesystem::path& directory, const std::string& run_name)
        : definition(readCase(directory)), output(readOutputControl(definition)),
          physics(definition.physics_file), solution(definition.solution_file),
          parallel(activeParallel()), mesh(readDistributedMesh(
              definition.mesh_file, parallel, definition.ghost_layers)),
          run_time(RunTime::forMesh(mesh, runtimeControl(
              definition, readMethodsFile(definition.methods_file))))
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
            if (value->name() == name && value->location() == location) {
                if (!read_file) {
                    throw std::logic_error("field " + name +
                        " was already created; initialization is only accepted on the first create*Field call");
                }
                for (const auto& source : field_sources) {
                    if (source.first != name || source.second == read_file) continue;
                    throw std::logic_error("field " + name + " was already " +
                        (source.second ? "loaded; use scalarField/vectorField/tensorField" :
                                         "created; use the create*Field API"));
                }
                return *value;
            }
        }
        for (const auto& info : field_names) {
            if (info == name) throw std::invalid_argument("field name reused with another type/location: " + name);
        }
        if (started) throw std::logic_error("create fields before starting the time/algorithm loop");
        auto value = std::make_unique<Field<T>>(mesh, location, name, initial);
        if (location == FieldLocation::Cell) {
            if (read_file)
                readFieldFile(definition.fields_directory / (name + ".field"), *value);
            // File fields follow the default all-fields policy. A program field
            // is opt-in, unless output.bs explicitly lists it in writeFields.
            if (((read_file && output.write_fields.empty()) ||
                 listed(output.write_fields, name)) &&
                !listed(output.exclude_fields, name)) output_names.push_back(name);
        }
        fields.push_back(std::move(value));
        field_names.push_back(name);
        field_sources.push_back({name, read_file});
        return *fields.back();
    }

    template <typename T>
    Field<T>& existing(std::vector<std::unique_ptr<Field<T>>>& fields,
                       const std::string& name, FieldLocation location) {
        requireName(name);
        for (const auto& value : fields)
            if (value->name() == name && value->location() == location) return *value;
        throw std::logic_error("field " + name +
            " has not been declared; use the load or create API first");
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
        const auto output_start = std::chrono::steady_clock::now();
        const auto cellField = [&](const std::string& name) {
            const auto hasCell = [&](const auto& fields) {
                for (const auto& field : fields)
                    if (field->name() == name && field->location() == FieldLocation::Cell) return true;
                return false;
            };
            return hasCell(scalars) || hasCell(vectors) || hasCell(tensors);
        };
        for (const auto& name : output.write_fields) {
            if (!cellField(name))
                throw std::invalid_argument("output writeFields requires an available cell field: " + name);
        }
        for (const auto& name : output.exclude_fields) {
            if (!cellField(name))
                throw std::invalid_argument("output excludeFields requires an available cell field: " + name);
        }
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
        run_time.recordOutput(std::chrono::duration<double>(
            std::chrono::steady_clock::now() - output_start).count());
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
    std::vector<std::pair<std::string, bool>> field_sources;
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
const Methods& Case::methods() const { return m_implementation->run_time.methods(); }
double Case::time() const { return m_implementation->run_time.time(); }
int Case::step() const { return m_implementation->run_time.step(); }

ScalarField& Case::scalarField(const std::string& name) {
    return m_implementation->field(m_implementation->scalars, name, FieldLocation::Cell);
}
VectorField& Case::vectorField(const std::string& name) {
    return m_implementation->field(m_implementation->vectors, name, FieldLocation::Cell);
}
ScalarField& Case::existingScalarField(const std::string& name) {
    return m_implementation->existing(m_implementation->scalars, name, FieldLocation::Cell);
}
VectorField& Case::existingVectorField(const std::string& name) {
    return m_implementation->existing(m_implementation->vectors, name, FieldLocation::Cell);
}
ScalarField& Case::createScalarField(const std::string& name, double initial) {
    return m_implementation->field(m_implementation->scalars, name, FieldLocation::Cell, false, initial);
}
VectorField& Case::createVectorField(const std::string& name, Vec3 initial) {
    return m_implementation->field(m_implementation->vectors, name, FieldLocation::Cell, false, initial);
}
TensorField& Case::tensorField(const std::string& name) {
    return m_implementation->field(m_implementation->tensors, name, FieldLocation::Cell);
}
TensorField& Case::existingTensorField(const std::string& name) {
    return m_implementation->existing(m_implementation->tensors, name, FieldLocation::Cell);
}
TensorField& Case::createTensorField(const std::string& name, Tensor3 initial) {
    return m_implementation->field(m_implementation->tensors, name, FieldLocation::Cell, false, initial);
}
ScalarField& Case::createFaceScalarField(const std::string& name) {
    return m_implementation->field(m_implementation->scalars, name, FieldLocation::Face, false);
}
VectorField& Case::createFaceVectorField(const std::string& name) {
    return m_implementation->field(m_implementation->vectors, name, FieldLocation::Face, false);
}
TensorField& Case::createFaceTensorField(const std::string& name) {
    return m_implementation->field(m_implementation->tensors, name, FieldLocation::Face, false);
}
ScalarField& Case::existingFaceField(const std::string& name) {
    return m_implementation->existing(m_implementation->scalars, name, FieldLocation::Face);
}
VectorField& Case::existingFaceVectorField(const std::string& name) {
    return m_implementation->existing(m_implementation->vectors, name, FieldLocation::Face);
}
TensorField& Case::existingFaceTensorField(const std::string& name) {
    return m_implementation->existing(m_implementation->tensors, name, FieldLocation::Face);
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
    methods().requireAllUsed();
}

void Case::start() {
    if (m_implementation->started) return;
    validate();
    m_implementation->started = true;
}

TimeStepper::TimeStepper(Case& problem):case_(&problem),options_(problem.timeControl()),
    value_(options_.start_time),dt_(options_.delta_t) { options_.validate(); }
TimeStepper time::start(Case& problem) { return TimeStepper(problem); }
void TimeStepper::advance() {
    auto& time = *this;
    if(time.finished()) throw std::logic_error("cannot advance beyond endTime");
    const double tolerance=32*std::numeric_limits<double>::epsilon()*
        std::max({std::abs(time.options_.start_time),std::abs(time.options_.end_time),time.options_.delta_t});
    double next=static_cast<double>(static_cast<long double>(time.options_.start_time)+
        static_cast<long double>(time.step_+1)*time.options_.delta_t);
    if(next>=time.options_.end_time-tolerance) next=time.options_.end_time;
    const double dt=next-time.value_;
    if(!(dt>0)) throw std::runtime_error("time step is below representable time precision");
    time.case_->setTime(next,time.step_+1,dt);
    time.value_=next; time.dt_=dt; ++time.step_;
}
void write(Case& problem,const TimeStepper& time) {
    if(&problem!=&time.owner()) throw std::invalid_argument("time service belongs to another case");
    problem.setTime(time.value(),time.step(),time.dt());
    problem.write();
}

TimeOptions readTimeControl(const Case& problem) {
    const auto& c=problem.timeControl(); return {c.start_time,c.end_time,c.delta_t};
}
namespace {
void numericalName(const std::string& name) {
    if (name.empty() || name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos)
        throw std::invalid_argument("invalid numerical configuration name: " + name);
}
EquationControl equationControl(const Case& problem, const std::string& name,
    std::initializer_list<std::string> terms) {
    numericalName(name);
    EquationControl result;
    result.name = name;
    result.spatial = problem.methods().equationOptions(name);
    result.spatial.requireComplete(name);
    for (const auto& term : terms) {
        numericalName(term);
        auto options = result.spatial;
        options.overlay(problem.methods().equationTermOptions(name, term));
        options.requireComplete(name + ".term." + term);
        if (!result.terms.emplace(term, options).second)
            throw std::invalid_argument("duplicate equation term " + name + "." + term);
    }
    const auto& p = problem.solution();
    auto& c = result.linear;
    const auto prefix = "equation." + name + ".";
    const auto solver = p.word(prefix + "solver");
    if (solver == "cg") c.solver = LinearSolverType::ConjugateGradient;
    else if (solver == "bicgstab") c.solver = LinearSolverType::BiCGSTAB;
    else throw std::invalid_argument(p.sourcePath().string() + ": unknown solver for " + name + ": " + solver);
    const auto preconditioner = p.word(prefix + "preconditioner");
    if (preconditioner == "none") c.preconditioner = PreconditionerType::None;
    else if (preconditioner == "incompleteCholesky") c.preconditioner = PreconditionerType::IncompleteCholesky;
    else if (preconditioner == "ilut") c.preconditioner = PreconditionerType::ILUT;
    else if (preconditioner == "amg") c.preconditioner = PreconditionerType::AlgebraicMultigrid;
    else throw std::invalid_argument(p.sourcePath().string() + ": unknown preconditioner for " + name + ": " + preconditioner);
    c.absolute_tolerance = p.number(prefix + "absoluteTolerance");
    c.relative_tolerance = p.number(prefix + "relativeTolerance");
    c.max_iterations = p.integer(prefix + "maxIterations");
    c.warm_start = p.boolean(prefix + "warmStart", c.warm_start);
    c.ilut_drop_tolerance = p.number(prefix + "ilutDropTolerance", c.ilut_drop_tolerance);
    c.ilut_fill_factor = p.integer(prefix + "ilutFillFactor", c.ilut_fill_factor);
    c.amg_max_levels = p.integer(prefix + "amgMaxLevels", c.amg_max_levels);
    c.amg_coarse_size = p.integer(prefix + "amgCoarseSize", c.amg_coarse_size);
    c.amg_smoothing_steps = p.integer(prefix + "amgSmoothingSteps", c.amg_smoothing_steps);
    c.amg_refresh_interval = p.integer(prefix + "amgRefreshInterval", c.amg_refresh_interval);
    c.validate();
    return result;
}
}
void Case::bindEquationIdentity(const std::string& name, const void* field, bool vector) const {
    auto& bindings = m_implementation->equation_bindings;
    const auto found = bindings.find(name);
    if (found != bindings.end() && found->second != std::make_pair(field, vector))
        throw std::invalid_argument("equation name bound to a different unknown: " + name);
    bindings.emplace(name, std::make_pair(field, vector));
}
EquationControl readEquationControl(const Case& problem, const std::string& name,
    const ScalarField& field, std::initializer_list<std::string> terms) {
    if (&field.mesh() != &problem.mesh()) throw std::invalid_argument("equation field belongs to a different case");
    problem.bindEquationIdentity(name, &field, false);
    return equationControl(problem,name,terms);
}
EquationControl readEquationControl(const Case& problem, const std::string& name,
    const VectorField& field, std::initializer_list<std::string> terms) {
    if (&field.mesh() != &problem.mesh()) throw std::invalid_argument("equation field belongs to a different case");
    problem.bindEquationIdentity(name, &field, true);
    return equationControl(problem,name,terms);
}
int readWriteInterval(const Case& problem) { return problem.outputControl().write_interval; }
void write(Case& problem,double value,int step) {
    problem.setTime(value,step,problem.timeControl().delta_t);
    problem.write();
}

const TimeControl& Case::timeControl() const { return m_implementation->run_time.timeControl(); }
const OutputControl& Case::outputControl() const { return m_implementation->output; }
void Case::setTime(double value, int step_value, double dt) {
    if (m_implementation->finished) throw std::logic_error("cannot change a finished case");
    start();
    m_implementation->run_time.setTime(value, step_value, dt);
}
void Case::write() {
    start();
    m_implementation->writeStep(true);
    // Explicit writes update both the time series and the latest written snapshot.
    // This does not mark the physical calculation converged or complete.
    m_implementation->write(m_implementation->final_directory);
}

void Case::finish() {
    Implementation& state = *m_implementation;
    if (state.finished) return;
    start();
    state.writeStep(true);
    if (state.final_directory != state.series_directory / timeName(time()))
        state.write(state.final_directory);
    state.finished = true;

}

PerformanceCounters Case::performance() const {
    const Implementation& state = *m_implementation;
    return maximumPerformance(state.run_time.performance(), state.parallel);
}

PerformanceCounters Case::localPerformance() const {
    return m_implementation->run_time.performance();
}

} // namespace babelsim
