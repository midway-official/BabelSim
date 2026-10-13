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
    if (!initialized)
        return {};
    int finalized = 0;
    detail::checkMpi(MPI_Finalized(&finalized), "MPI_Finalized");
    if (finalized)
        throw std::logic_error("cannot create Case after MPI_Finalize");
    return ParallelContext::world();
}

const char* boundaryTypeName(BoundaryType type) {
    switch (type) {
    case BoundaryType::FixedValue:
        return "fixedValue";
    case BoundaryType::FixedGradient:
        return "fixedGradient";
    case BoundaryType::ZeroGradient:
        return "zeroGradient";
    case BoundaryType::InletOutlet:
        return "inletOutlet";
    case BoundaryType::Symmetry:
        return "symmetry";
    }
    throw std::invalid_argument("invalid field boundary type");
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
    if (path.empty() || path.is_absolute())
        throw std::invalid_argument("output path must be relative");
    for (const auto& part : path) {
        if (part == ".." || part == ".")
            throw std::invalid_argument("output path cannot contain . or ..");
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
    try {
        std::stod(name, &consumed);
    } catch (const std::exception&) {
    }
    if (consumed == name.size() || name == "all" || name == "latest")
        throw std::invalid_argument("run/final label must not be a time or reserved selection: " +
                                    name);
}

std::string timeName(double time) {
    std::ostringstream text;
    text << std::setprecision(15) << time;
    return text.str();
}

PerformanceCounters maximumPerformance(const PerformanceCounters& local,
                                       const ParallelContext& parallel) {
    // Member pointers keep packing and unpacking tied to the same field list.
    // Counters retain the existing double reduction representation.
    static constexpr auto counters = std::array{&PerformanceCounters::linear_solves,
                                                &PerformanceCounters::krylov_iterations,
                                                &PerformanceCounters::sparse_matvecs,
                                                &PerformanceCounters::halo_exchanges,
                                                &PerformanceCounters::halo_bytes,
                                                &PerformanceCounters::global_reductions,
                                                &PerformanceCounters::equation_assemblies,
                                                &PerformanceCounters::preconditioner_setups,
                                                &PerformanceCounters::preconditioner_applications,
                                                &PerformanceCounters::output_writes,
                                                &PerformanceCounters::matrix_pattern_builds,
                                                &PerformanceCounters::matrix_value_updates,
                                                &PerformanceCounters::rhs_only_solves,
                                                &PerformanceCounters::true_residual_checks,
                                                &PerformanceCounters::preconditioner_reuses,
                                                &PerformanceCounters::preconditioner_retries,
                                                &PerformanceCounters::reported_iteration_intervals};
    static constexpr auto timers = std::array{&PerformanceCounters::elapsed_seconds,
                                              &PerformanceCounters::assembly_seconds,
                                              &PerformanceCounters::preconditioner_seconds,
                                              &PerformanceCounters::preconditioner_apply_seconds,
                                              &PerformanceCounters::linear_solve_seconds,
                                              &PerformanceCounters::sparse_matvec_seconds,
                                              &PerformanceCounters::halo_seconds,
                                              &PerformanceCounters::global_reduction_seconds,
                                              &PerformanceCounters::output_seconds,
                                              &PerformanceCounters::matrix_update_seconds,
                                              &PerformanceCounters::rhs_seconds,
                                              &PerformanceCounters::residual_check_seconds,
                                              &PerformanceCounters::pattern_build_seconds,
                                              &PerformanceCounters::iteration_compute_seconds};
    std::array<double, counters.size() + timers.size()> local_values{};
    std::size_t index = 0;
    for (auto member : counters)
        local_values[index++] = static_cast<double>(local.*member);
    for (auto member : timers)
        local_values[index++] = local.*member;
    std::array<double, local_values.size()> global_values{};
    parallel.maximum(local_values.data(), global_values.data(),
                     detail::mpiCount(global_values.size(), "performance counters"));
    PerformanceCounters global;
    index = 0;
    for (auto member : counters)
        global.*member = static_cast<std::uint64_t>(global_values[index++]);
    for (auto member : timers)
        global.*member = global_values[index++];
    global.equation_systems = local.equation_systems;
    return global;
}

} // 匿名命名空间

struct Case::Implementation {
    struct EquationBinding {
        const void* field;
        bool vector;
        std::vector<std::string> terms;
    };
    mutable std::map<std::string, EquationBinding> equation_bindings;
    mutable std::map<std::string, ModuleDescription> modules;
    std::map<std::string, std::string> correction_sources;
    void requireInitialization() const {
        if (started)
            throw std::logic_error("declarations are closed by finalizeInitialization()");
    }
    void requireReady() const {
        if (!started)
            throw std::logic_error(
                "call finalizeInitialization() before advancing time or writing results");
    }
    Implementation(const std::filesystem::path& directory, const std::string& run_name)
        : definition(readCase(directory)), output(readOutputControl(definition)),
          physics(definition.physics_file), solution(definition.solution_file),
          parallel(activeParallel()),
          mesh(readDistributedMesh(definition.mesh_file, parallel, definition.ghost_layers)),
          run_time(RunTime::forMesh(
              mesh, runtimeControl(definition, readMethodsFile(definition.methods_file)))) {
        requireRelativePath(output.directory);
        requireRunName(output.time_name);
        if (!run_name.empty())
            requireRunName(run_name);
        const auto base = definition.root / output.directory;
        // 命名运行隔离时间序列，避免不同进程数/参数的实验互相覆盖。
        series_directory = run_name.empty() ? base : base / run_name;
        final_directory = base / (run_name.empty() ? output.time_name : run_name);
    }

    template <typename T>
    Field<T>& field(std::vector<std::unique_ptr<Field<T>>>& fields, const std::string& name,
                    FieldLocation location, bool read_file = true, T initial = T{}) {
        requireName(name);
        for (const auto& value : fields) {
            if (value->name() == name && value->location() == location) {
                if (!read_file) {
                    throw std::logic_error(
                        "field " + name +
                        " already exists; use problem.fields to access it without reinitializing");
                }
                for (const auto& source : field_sources) {
                    if (source.first == name && !source.second)
                        throw std::logic_error("field " + name +
                                               " was created in code; use problem.fields to access "
                                               "it without loading a file");
                }
                return *value;
            }
        }
        for (const auto& info : field_names) {
            if (info == name)
                throw std::invalid_argument("field name reused with another type/location: " +
                                            name);
        }
        requireInitialization();
        auto value = std::make_unique<Field<T>>(mesh, location, name, initial);
        if (location == FieldLocation::Cell) {
            if (read_file)
                readFieldFile(definition.fields_directory / (name + ".field"), *value);
            // File fields follow the default all-fields policy. A program field
            // is opt-in, unless output.bs explicitly lists it in writeFields.
            if (((read_file && output.write_fields.empty()) || listed(output.write_fields, name)) &&
                !listed(output.exclude_fields, name))
                output_names.push_back(name);
        }
        fields.push_back(std::move(value));
        field_names.push_back(name);
        field_sources.push_back({name, read_file});
        return *fields.back();
    }

    template <typename T>
    Field<T>& existing(std::vector<std::unique_ptr<Field<T>>>& fields, const std::string& name,
                       FieldLocation location) {
        requireName(name);
        for (const auto& value : fields)
            if (value->name() == name && value->location() == location)
                return *value;
        throw std::logic_error("field " + name +
                               " has not been declared; use problem.load or problem.create first");
    }

    template <typename T>
    void writeFields(const std::filesystem::path& directory,
                     const std::vector<std::unique_ptr<Field<T>>>& fields, const char* type,
                     std::vector<FieldOutputInfo>& info) {
        for (const auto& field : fields) {
            if (std::find(output_names.begin(), output_names.end(), field->name()) ==
                output_names.end())
                continue;
            writeOwnedFieldCsv(directory, *field, parallel);
            info.push_back({field->name(), type, field->location()});
        }
    }

    void write(const std::filesystem::path& directory) {
        const auto output_start = std::chrono::steady_clock::now();
        const auto cellField = [&](const std::string& name) {
            const auto hasCell = [&](const auto& fields) {
                for (const auto& field : fields)
                    if (field->name() == name && field->location() == FieldLocation::Cell)
                        return true;
                return false;
            };
            return hasCell(scalars) || hasCell(vectors) || hasCell(tensors);
        };
        for (const auto& name : output.write_fields) {
            if (!cellField(name))
                throw std::invalid_argument(
                    "output writeFields requires an available cell field: " + name);
        }
        for (const auto& name : output.exclude_fields) {
            if (!cellField(name))
                throw std::invalid_argument(
                    "output excludeFields requires an available cell field: " + name);
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
                            entry.path().filename().string().rfind("rank-", 0) == 0)
                            ++ranks;
                    }
                }
                if (ranks != 0 && ranks != parallel.size)
                    incompatible = 1;
            } catch (const std::exception&) {
                incompatible = 1;
            }
        }
        if (parallel.maximum(incompatible))
            throw std::runtime_error(
                "existing result has incompatible partitions; use a new run label: " +
                directory.string());
        std::vector<FieldOutputInfo> info;
        writeFields(directory, scalars, "scalar", info);
        writeFields(directory, vectors, "vector", info);
        writeFields(directory, tensors, "tensor", info);
        writeOwnedResultMetadata(directory, mesh, parallel, timeName(run_time.time()), info);
        run_time.recordOutput(
            std::chrono::duration<double>(std::chrono::steady_clock::now() - output_start).count());
    }

    void writeStep(bool force) {
        if (run_time.step() == last_written_step)
            return;
        if (!force && run_time.step() % output.write_interval != 0)
            return;
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
    : load(*this), create(*this), fields(*this), read(*this),
      m_implementation(std::make_unique<Implementation>(directory, run_name)) {}
Case::~Case() = default;

const std::string& Case::solver() const {
    return m_implementation->definition.solver;
}
const CaseDefinition& Case::definition() const {
    return m_implementation->definition;
}
ModuleParameters Case::Read::parameters(const ModuleDescription& description) const {
    auto& state = *owner_.m_implementation;
    state.requireInitialization();
    auto result = readModuleParameters(owner_, description);
    if (!description.time_methods.empty() &&
        std::find(description.time_methods.begin(), description.time_methods.end(),
                  methods().time) == description.time_methods.end())
        throw std::invalid_argument("module " + description.name +
                                    " does not support the configured time scheme");
    if (!description.compatible_solvers.empty() &&
        !listed(description.compatible_solvers, owner_.solver()))
        throw std::invalid_argument("module " + description.name + " is incompatible with solver " +
                                    owner_.solver());
    state.modules.emplace(description.name, description);
    return result;
}
const Mesh& Case::mesh() const {
    return m_implementation->mesh;
}
const Parameters& Case::Read::physics() const {
    return owner_.m_implementation->physics;
}
const Parameters& Case::Read::solution() const {
    return owner_.m_implementation->solution;
}
const Methods& Case::Read::methods() const {
    return owner_.m_implementation->run_time.methods();
}
double Case::time() const {
    return m_implementation->run_time.time();
}
int Case::step() const {
    return m_implementation->run_time.step();
}

ScalarField& Case::Load::scalarField(const std::string& name) {
    return owner_.m_implementation->field(owner_.m_implementation->scalars, name,
                                          FieldLocation::Cell);
}
VectorField& Case::Load::vectorField(const std::string& name) {
    return owner_.m_implementation->field(owner_.m_implementation->vectors, name,
                                          FieldLocation::Cell);
}
ScalarField& Case::Fields::scalarField(const std::string& name) {
    return owner_.m_implementation->existing(owner_.m_implementation->scalars, name,
                                             FieldLocation::Cell);
}
VectorField& Case::Fields::vectorField(const std::string& name) {
    return owner_.m_implementation->existing(owner_.m_implementation->vectors, name,
                                             FieldLocation::Cell);
}
ScalarField& Case::Create::scalarField(const std::string& name, double initial) {
    return owner_.m_implementation->field(owner_.m_implementation->scalars, name,
                                          FieldLocation::Cell, false, initial);
}
VectorField& Case::Create::vectorField(const std::string& name, Vec3 initial) {
    return owner_.m_implementation->field(owner_.m_implementation->vectors, name,
                                          FieldLocation::Cell, false, initial);
}
TensorField& Case::Load::tensorField(const std::string& name) {
    return owner_.m_implementation->field(owner_.m_implementation->tensors, name,
                                          FieldLocation::Cell);
}
TensorField& Case::Fields::tensorField(const std::string& name) {
    return owner_.m_implementation->existing(owner_.m_implementation->tensors, name,
                                             FieldLocation::Cell);
}
TensorField& Case::Create::tensorField(const std::string& name, Tensor3 initial) {
    return owner_.m_implementation->field(owner_.m_implementation->tensors, name,
                                          FieldLocation::Cell, false, initial);
}
ScalarField& Case::Create::faceScalarField(const std::string& name, double initial) {
    return owner_.m_implementation->field(owner_.m_implementation->scalars, name,
                                          FieldLocation::Face, false, initial);
}
VectorField& Case::Create::faceVectorField(const std::string& name, Vec3 initial) {
    return owner_.m_implementation->field(owner_.m_implementation->vectors, name,
                                          FieldLocation::Face, false, initial);
}
TensorField& Case::Create::faceTensorField(const std::string& name, Tensor3 initial) {
    return owner_.m_implementation->field(owner_.m_implementation->tensors, name,
                                          FieldLocation::Face, false, initial);
}
ScalarField& Case::Fields::faceScalarField(const std::string& name) {
    return owner_.m_implementation->existing(owner_.m_implementation->scalars, name,
                                             FieldLocation::Face);
}
VectorField& Case::Fields::faceVectorField(const std::string& name) {
    return owner_.m_implementation->existing(owner_.m_implementation->vectors, name,
                                             FieldLocation::Face);
}
TensorField& Case::Fields::faceTensorField(const std::string& name) {
    return owner_.m_implementation->existing(owner_.m_implementation->tensors, name,
                                             FieldLocation::Face);
}

ScalarField& Case::Create::homogeneousCorrection(const ScalarField& original,
                                                 const std::string& name) {
    if (&original.mesh() != &owner_.mesh() || original.location() != FieldLocation::Cell)
        throw std::invalid_argument(
            "homogeneous correction requires a cell field on this Case mesh");
    original.validateStorage();
    auto& correction = scalarField(name);
    setHomogeneousCorrectionBoundaries(correction, original);
    owner_.m_implementation->correction_sources.emplace(name, original.name());
    return correction;
}

void Case::selectOutput(const std::string& name, const void* field, bool enabled) {
    Implementation& state = *m_implementation;
    bool owned_cell_field = false;
    const auto check = [&](const auto& registered) {
        for (const auto& value : registered)
            if (value.get() == field && value->location() == FieldLocation::Cell)
                owned_cell_field = true;
    };
    check(state.scalars);
    check(state.vectors);
    check(state.tensors);
    if (!owned_cell_field)
        throw std::invalid_argument("output selection requires a Case-owned cell field");
    auto found = std::find(state.output_names.begin(), state.output_names.end(), name);
    if (enabled && found == state.output_names.end())
        state.output_names.push_back(name);
    if (!enabled && found != state.output_names.end())
        state.output_names.erase(found);
}

void Case::finalizeInitialization() {
    auto& state = *m_implementation;
    if (state.started)
        return;
    for (const auto& item : state.modules) {
        const auto& module = item.second;
        for (const auto& described : module.fields) {
            const void* identity = nullptr;
            const auto check = [&](const auto& registered, DescribedFieldType type) {
                for (const auto& field : registered) {
                    if (field->name() != described.name)
                        continue;
                    const bool cell = field->location() == FieldLocation::Cell;
                    if (type != described.type ||
                        cell != (described.location == DescribedFieldLocation::Cell))
                        throw std::invalid_argument(
                            "module " + module.name +
                            ": field type/location differs: " + described.name);
                    field->validateStorage();
                    if (cell && !field->calculatedBoundary() &&
                        described.initialization != FieldInitialization::HomogeneousCorrection)
                        for (Index patch = 0; patch < state.mesh.patchCount(); ++patch)
                            if (!module.boundary_types.empty() &&
                                !listed(module.boundary_types,
                                        boundaryTypeName(field->boundary(patch).type)))
                                throw std::invalid_argument("module " + module.name +
                                                            ": unsupported boundary for " +
                                                            described.name);
                    identity = field.get();
                }
            };
            check(state.scalars, DescribedFieldType::Scalar);
            check(state.vectors, DescribedFieldType::Vector);
            check(state.tensors, DescribedFieldType::Tensor);
            if (!identity)
                throw std::logic_error("module " + module.name +
                                       ": required field is not declared: " + described.name);
            const auto source =
                std::find_if(state.field_sources.begin(), state.field_sources.end(),
                             [&](const auto& value) { return value.first == described.name; });
            if ((described.initialization == FieldInitialization::File && !source->second) ||
                (described.initialization != FieldInitialization::File &&
                 described.initialization != FieldInitialization::Provided && source->second))
                throw std::invalid_argument(
                    "module " + module.name +
                    ": field initialization source differs: " + described.name);
            if (described.initialization == FieldInitialization::HomogeneousCorrection) {
                const auto found = state.correction_sources.find(described.name);
                if (found == state.correction_sources.end() ||
                    found->second != described.source_field)
                    throw std::invalid_argument("module " + module.name +
                                                ": use create.homogeneousCorrection for " +
                                                described.name);
            }
            for (const auto& equation : module.equations) {
                if (equation.unknown != described.name)
                    continue;
                const auto found = state.equation_bindings.find(equation.name);
                if (found == state.equation_bindings.end())
                    throw std::logic_error("module " + module.name +
                                           ": equation is not bound: " + equation.name);
                if (found->second.field != identity)
                    throw std::invalid_argument("module " + module.name +
                                                ": equation unknown differs: " + equation.name);
                for (const auto& term : equation.terms)
                    if (!listed(found->second.terms, term))
                        throw std::invalid_argument(
                            "module " + module.name +
                            ": equation term is not bound: " + equation.name + "." + term);
            }
        }
    }
    read.physics().requireAllUsed();
    read.solution().requireAllUsed();
    read.methods().requireAllUsed();
    state.started = true;
}

TimeStepper::TimeStepper(Case& problem)
    : case_(&problem), options_(problem.read.timeControl()), value_(options_.start_time),
      dt_(options_.delta_t) {
    options_.validate();
}
TimeStepper time::makeStepper(Case& problem) {
    return TimeStepper(problem);
}
const Mesh& TimeStepper::mesh() const {
    return case_->mesh();
}
void TimeStepper::advance() {
    advanceImpl(options_.delta_t, false);
}
void TimeStepper::advance(double requestedDt) {
    advanceImpl(requestedDt, true);
}
void TimeStepper::advanceImpl(double requestedDt, bool explicitStep) {
    auto& time = *this;
    if (time.attempt_active_)
        throw std::logic_error("resolve the pending step before advancing the clock");
    if (!(requestedDt > 0) || !std::isfinite(requestedDt))
        throw std::invalid_argument("requested step requires positive finite dt");
    if (time.finished())
        throw std::logic_error("cannot advance beyond endTime");
    const double tolerance = 32 * std::numeric_limits<double>::epsilon() *
                             std::max({std::abs(time.options_.start_time),
                                       std::abs(time.options_.end_time), time.options_.delta_t});
    const bool variable = time.variable_steps_ || explicitStep;
    double next =
        variable
            ? static_cast<double>(static_cast<long double>(time.value_) + requestedDt)
            : static_cast<double>(static_cast<long double>(time.options_.start_time) +
                                  static_cast<long double>(time.step_ + 1) * time.options_.delta_t);
    if (next >= time.options_.end_time - tolerance)
        next = time.options_.end_time;
    const double dt = next - time.value_;
    if (!(dt > 0))
        throw std::runtime_error("time step is below representable time precision");
    time.case_->setTime(next, time.step_ + 1, dt);
    time.value_ = next;
    time.dt_ = dt;
    time.variable_steps_ = variable;
    ++time.step_;
}
void write(Case& problem, const TimeStepper& time) {
    if (&problem != &time.owner())
        throw std::invalid_argument("time service belongs to another case");
    problem.setTime(time.value(), time.step(), time.dt());
    problem.write();
}

TimeOptions readTimeControl(const Case& problem) {
    const auto& c = problem.read.timeControl();
    return {c.start_time, c.end_time, c.delta_t};
}
namespace {
void numericalName(const std::string& name) {
    if (name.empty() || name.find_first_not_of(
                            "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") !=
                            std::string::npos)
        throw std::invalid_argument("invalid numerical configuration name: " + name);
}
template <class Terms>
EquationControl equationControl(const Case& problem, const std::string& name, const Terms& terms) {
    numericalName(name);
    EquationControl result;
    result.name = name;
    result.spatial = problem.read.methods().equationOptions(name);
    result.spatial.requireComplete(name);
    for (const auto& term : terms) {
        numericalName(term);
        auto options = result.spatial;
        options.overlay(problem.read.methods().equationTermOptions(name, term));
        options.requireComplete(name + ".term." + term);
        if (!result.terms.emplace(term, options).second)
            throw std::invalid_argument("duplicate equation term " + name + "." + term);
    }
    const auto& p = problem.read.solution();
    auto& c = result.linear;
    const auto prefix = "equation." + name + ".";
    const auto solver = p.word(prefix + "kspType");
    if (solver == "cg")
        c.solver = LinearSolverType::ConjugateGradient;
    else if (solver == "bcgs")
        c.solver = LinearSolverType::BiCGSTAB;
    else if (solver == "gmres")
        c.solver = LinearSolverType::GMRES;
    else if (solver == "fgmres")
        c.solver = LinearSolverType::FGMRES;
    else
        throw std::invalid_argument(p.sourcePath().string() + ": unknown PETSc KSP type for " +
                                    name + ": " + solver);
    const auto preconditioner = p.word(prefix + "pcType");
    if (preconditioner == "none")
        c.preconditioner = PreconditionerType::None;
    else if (preconditioner == "icc")
        c.preconditioner = PreconditionerType::IncompleteCholesky;
    else if (preconditioner == "hypre")
        c.preconditioner = PreconditionerType::Hypre;
    else if (preconditioner == "gamg")
        c.preconditioner = PreconditionerType::GAMG;
    else if (preconditioner == "bjacobi")
        c.preconditioner = PreconditionerType::BlockJacobi;
    else if (preconditioner == "asm")
        c.preconditioner = PreconditionerType::ASM;
    else if (preconditioner == "jacobi")
        c.preconditioner = PreconditionerType::Jacobi;
    else
        throw std::invalid_argument(p.sourcePath().string() + ": unknown PETSc PC type for " +
                                    name + ": " + preconditioner);
    c.absolute_tolerance = p.number(prefix + "absoluteTolerance");
    c.relative_tolerance = p.number(prefix + "relativeTolerance");
    c.max_iterations = p.integer(prefix + "maxIterations");
    c.warm_start = p.boolean(prefix + "warmStart", c.warm_start);
    c.preconditioner_rebuild_interval =
        p.integer(prefix + "preconditionerRebuildInterval", c.preconditioner_rebuild_interval);
    c.preconditioner_reuse_iteration_growth = p.number(
        prefix + "preconditionerReuseIterationGrowth", c.preconditioner_reuse_iteration_growth);
    c.amg_max_levels = p.integer(prefix + "amgMaxLevels", c.amg_max_levels);
    c.amg_coarse_size = p.integer(prefix + "amgCoarseSize", c.amg_coarse_size);
    c.amg_smoothing_steps = p.integer(prefix + "amgSmoothingSteps", c.amg_smoothing_steps);
    c.validate();
    return result;
}
}
void Case::bindEquationIdentity(const std::string& name, const void* field, bool vector,
                                const std::vector<std::string>& terms) const {
    auto& state = *m_implementation;
    state.requireInitialization();
    auto& bindings = state.equation_bindings;
    const auto found = bindings.find(name);
    if (found != bindings.end()) {
        if (found->second.field != field || found->second.vector != vector)
            throw std::invalid_argument("equation name bound to a different unknown: " + name);
        for (const auto& term : terms)
            if (!listed(found->second.terms, term))
                found->second.terms.push_back(term);
    } else
        bindings.emplace(name, Implementation::EquationBinding{field, vector, terms});
}
template <class T, class Terms>
EquationControl Case::bindEquationControl(const std::string& name, const Field<T>& field,
                                          const Terms& terms,
                                          const std::string& expectedUnknown) const {
    m_implementation->requireInitialization();
    if (&field.mesh() != &mesh() || field.location() != FieldLocation::Cell)
        throw std::invalid_argument("equation " + name +
                                    " requires a cell field on this Case mesh");
    if (!expectedUnknown.empty() && field.name() != expectedUnknown)
        throw std::invalid_argument("equation " + name + " requires Case cell unknown " +
                                    expectedUnknown);
    auto result = equationControl(*this, name, terms);
    bindEquationIdentity(name, &field, std::is_same_v<T, Vec3>, terms);
    return result;
}
EquationControl readEquationControl(const Case& problem, const std::string& name,
                                    const ScalarField& field,
                                    std::initializer_list<std::string> terms) {
    return problem.bindEquationControl(name, field, terms);
}
EquationControl readEquationControl(const Case& problem, const std::string& name,
                                    const VectorField& field,
                                    std::initializer_list<std::string> terms) {
    return problem.bindEquationControl(name, field, terms);
}
EquationControl readEquationControl(const Case& problem, const EquationDescription& description,
                                    const ScalarField& field) {
    return problem.bindEquationControl(description.name, field, description.terms,
                                       description.unknown);
}
EquationControl readEquationControl(const Case& problem, const EquationDescription& description,
                                    const VectorField& field) {
    return problem.bindEquationControl(description.name, field, description.terms,
                                       description.unknown);
}
int readWriteInterval(const Case& problem) {
    return problem.read.outputControl().write_interval;
}
void write(Case& problem, double value, int step) {
    problem.setTime(value, step, problem.read.timeControl().delta_t);
    problem.write();
}

const TimeControl& Case::Read::timeControl() const {
    return owner_.m_implementation->run_time.timeControl();
}
const OutputControl& Case::Read::outputControl() const {
    return owner_.m_implementation->output;
}
void Case::setTime(double value, int step_value, double dt) {
    if (m_implementation->finished)
        throw std::logic_error("cannot change a finished case");
    m_implementation->requireReady();
    m_implementation->run_time.setTime(value, step_value, dt);
}
void Case::write() {
    m_implementation->requireReady();
    m_implementation->writeStep(true);
    // Explicit writes update both the time series and the latest written snapshot.
    // This does not mark the physical calculation converged or complete.
    m_implementation->write(m_implementation->final_directory);
}

void Case::finish() {
    Implementation& state = *m_implementation;
    if (state.finished)
        return;
    m_implementation->requireReady();
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
