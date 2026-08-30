#include "babelsim/runtime.h"

#include "babelsim/assembly.h"
#include "babelsim/distributed_solver.h"
#include "babelsim/equation.h"
#include "babelsim/linear_solver.h"
#include "babelsim/mpi_support.h"
#include "babelsim/operators.h"
#include "babelsim/parallel.h"
#include "internal/scalar_equation_control.h"
#include "internal/vector_equation_control.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace babelsim {
namespace {

template <typename EquationType>
void prepare(
    const Eigen::SparseMatrix<double>& matrix,
    const EquationType& equation,
    bool& pattern_ready,
    PreparedLinearSolver& serial_solver,
    DistributedLinearSolver* distributed_solver)
{
    if (distributed_solver != nullptr) {
        if (pattern_ready) distributed_solver->factorize(matrix, equation);
        else distributed_solver->compute(matrix, equation);
    } else if (pattern_ready) {
        serial_solver.factorize(matrix);
    } else {
        serial_solver.compute(matrix);
    }
    pattern_ready = true;
}

template <typename T>
void requireCellField(const Field<T>& field, const Mesh& mesh, const char* name) {
    field.validateStorage();
    if (&field.mesh() != &mesh || field.location() != FieldLocation::Cell) {
        throw std::invalid_argument(std::string(name) + " must be a cell field on the run mesh");
    }
}

void requireFaceField(const ScalarField& field, const Mesh& mesh, const char* name) {
    field.validateStorage();
    if (&field.mesh() != &mesh || field.location() != FieldLocation::Face) {
        throw std::invalid_argument(std::string(name) + " must be a face field on the run mesh");
    }
}

int residualSign(bool lhs, int term_sign) {
    return (lhs ? 1 : -1) * term_sign;
}

thread_local RunTime* active_run_time = nullptr;

}  // 匿名命名空间

struct RunTime::Implementation {
    struct ScalarHistory {
        explicit ScalarHistory(const ScalarField& field)
            : source(&field), previous(field), older(field)
        {}

        const ScalarField* source;
        ScalarField previous;
        ScalarField older;
        bool has_older = false;

        void advance(const ScalarField& current) {
            for (std::size_t index = 0; index < previous.size(); ++index) {
                older[static_cast<Index>(index)] = previous[static_cast<Index>(index)];
                previous[static_cast<Index>(index)] = current[static_cast<Index>(index)];
            }
            has_older = true;
        }
    };

    struct VectorHistory {
        explicit VectorHistory(const VectorField& field)
            : source(&field), previous(field), older(field)
        {}

        const VectorField* source;
        VectorField previous;
        VectorField older;
        bool has_older = false;

        void advance(const VectorField& current) {
            for (std::size_t index = 0; index < previous.size(); ++index) {
                older[static_cast<Index>(index)] = previous[static_cast<Index>(index)];
                previous[static_cast<Index>(index)] = current[static_cast<Index>(index)];
            }
            has_older = true;
        }
    };

    Implementation(const Mesh& mesh_value, RuntimeControl settings, ParallelContext parallel_value)
        : mesh(&mesh_value),
          control(std::move(settings)),
          parallel(std::move(parallel_value)),
          scalar_assembly(mesh_value),
          vector_assembly(mesh_value),
          scalar_linear_solver(control.scalar_solver),
          vector_linear_solver(control.vector_solver),
          scalar_source(Eigen::VectorXd::Zero(mesh_value.ownedCellCount())),
          scalar_solution(Eigen::VectorXd::Zero(mesh_value.ownedCellCount())),
          gradient_workspace(mesh_value, FieldLocation::Cell, "grad"),
          face_coefficient_workspace(mesh_value, FieldLocation::Face, "faceCoefficient")
    {
        mesh->validate();
        control.validate();
        parallel.validate();
        if (parallel.distributed() != (mesh->ownedCellCount() < mesh->cellCount())) {
            throw std::invalid_argument("run time and mesh ownership are inconsistent");
        }
        if (parallel.distributed()) {
            halo = std::make_unique<HaloExchange>(*mesh, parallel);
            scalar_distributed_solver = std::make_unique<DistributedLinearSolver>(
                *mesh, parallel, control.scalar_solver);
            vector_distributed_solver = std::make_unique<DistributedLinearSolver>(
                *mesh, parallel, control.vector_solver);
        }
        for (std::size_t component = 0; component < 3; ++component) {
            vector_source[component].resize(mesh->ownedCellCount());
            vector_solution[component].resize(mesh->ownedCellCount());
        }
    }

    ScalarHistory& history(const ScalarField& field) {
        for (ScalarHistory& value : scalar_histories) {
            if (value.source == &field) return value;
        }
        scalar_histories.emplace_back(field);
        return scalar_histories.back();
    }

    VectorHistory& history(const VectorField& field) {
        for (VectorHistory& value : vector_histories) {
            if (value.source == &field) return value;
        }
        vector_histories.emplace_back(field);
        return vector_histories.back();
    }

    void synchronize(ScalarField& field) {
        requireCellOrFace(field);
        if (halo) halo->exchange(field);
    }

    void synchronize(VectorField& field) {
        requireCellOrFace(field);
        if (halo) halo->exchange(field);
    }

    void synchronize(TensorField& field) {
        requireCellOrFace(field);
        if (halo) halo->exchange(field);
    }

    template <typename T>
    void requireCellOrFace(const Field<T>& field) const {
        field.validateStorage();
        if (&field.mesh() != mesh ||
            (field.location() != FieldLocation::Cell && field.location() != FieldLocation::Face)) {
            throw std::invalid_argument("field does not belong to the run mesh");
        }
    }

    const Mesh* mesh;
    RuntimeControl control;
    ParallelContext parallel;
    std::unique_ptr<HaloExchange> halo;
    SparseAssembly scalar_assembly;
    SparseAssembly vector_assembly;
    PreparedLinearSolver scalar_linear_solver;
    PreparedLinearSolver vector_linear_solver;
    std::unique_ptr<DistributedLinearSolver> scalar_distributed_solver;
    std::unique_ptr<DistributedLinearSolver> vector_distributed_solver;
    Eigen::VectorXd scalar_source;
    Eigen::VectorXd scalar_solution;
    std::array<Eigen::VectorXd, 3> vector_source;
    std::array<Eigen::VectorXd, 3> vector_solution;
    bool scalar_pattern_ready = false;
    bool vector_pattern_ready = false;
    std::vector<ScalarHistory> scalar_histories;
    std::vector<VectorHistory> vector_histories;
    VectorField gradient_workspace;
    ScalarField face_coefficient_workspace;
    double current_time = 0.0;
    int current_step = 0;
};

void RuntimeControl::validate() const {
    time.validate();
    scalar_solver.validate();
    vector_solver.validate();
}

RunTime::RunTime(const Mesh& mesh, RuntimeControl control)
    : m_implementation(nullptr)
{
    if (active_run_time != nullptr) {
        throw std::logic_error(
            "only one RunTime may be active in a solver thread; destroy the previous run first");
    }
    int initialized = 0;
    detail::checkMpi(MPI_Initialized(&initialized), "MPI_Initialized");
    ParallelContext parallel;
    if (initialized != 0) {
        int finalized = 0;
        detail::checkMpi(MPI_Finalized(&finalized), "MPI_Finalized");
        if (finalized != 0) throw std::logic_error("cannot create RunTime after MPI_Finalize");
        parallel = ParallelContext::world();
    }
    m_implementation = std::make_unique<Implementation>(mesh, std::move(control), parallel);
    m_implementation->current_time = m_implementation->control.time.start_time;
    active_run_time = this;
}

RunTime RunTime::forMesh(const Mesh& mesh, RuntimeControl control) {
    return RunTime(mesh, std::move(control));
}

RunTime::~RunTime() {
    if (active_run_time == this) active_run_time = nullptr;
}

RunTime& RunTime::current() {
    if (active_run_time == nullptr) {
        throw std::logic_error("solve/fvc/diagnostics require an active RunTime");
    }
    return *active_run_time;
}

const Mesh& RunTime::mesh() const { return *m_implementation->mesh; }
const Methods& RunTime::methods() const { return m_implementation->control.methods; }
double RunTime::time() const { return m_implementation->current_time; }
double RunTime::deltaT() const { return m_implementation->control.time.delta_t; }
int RunTime::step() const { return m_implementation->current_step; }
bool RunTime::primary() const { return m_implementation->parallel.rank == 0; }

void RunTime::synchronize(ScalarField& field) { m_implementation->synchronize(field); }
void RunTime::synchronize(VectorField& field) { m_implementation->synchronize(field); }
void RunTime::synchronize(TensorField& field) { m_implementation->synchronize(field); }

double RunTime::relativeChange(
    const VectorField& current,
    const VectorField& previous) const
{
    const Implementation& state = *m_implementation;
    requireCellField(current, *state.mesh, "relative-change current");
    requireCellField(previous, *state.mesh, "relative-change previous");
    double difference_squared = 0.0;
    double current_squared = 0.0;
    for (Index cell : state.mesh->owned_cells) {
        difference_squared += squaredNorm(current[cell] - previous[cell]);
        current_squared += squaredNorm(current[cell]);
    }
    const double local[2] = {difference_squared, current_squared};
    double global[2]{};
    state.parallel.sum(local, global, 2);
    return std::sqrt(global[0]) / std::max(std::sqrt(global[1]), 1e-30);
}

double RunTime::relativeChange(
    const ScalarField& current,
    const ScalarField& previous) const
{
    const Implementation& state = *m_implementation;
    requireCellField(current, *state.mesh, "relative-change current");
    requireCellField(previous, *state.mesh, "relative-change previous");
    double difference_squared = 0.0;
    double current_squared = 0.0;
    for (Index cell : state.mesh->owned_cells) {
        const double difference = current[cell] - previous[cell];
        difference_squared += difference * difference;
        current_squared += current[cell] * current[cell];
    }
    const double local[2] = {difference_squared, current_squared};
    double global[2]{};
    state.parallel.sum(local, global, 2);
    return std::sqrt(global[0]) / std::max(std::sqrt(global[1]), 1e-30);
}

double RunTime::relativeMagnitude(
    const ScalarField& value,
    const ScalarField& reference) const
{
    const Implementation& state = *m_implementation;
    requireCellField(value, *state.mesh, "relative-magnitude value");
    requireCellField(reference, *state.mesh, "relative-magnitude reference");
    double value_squared = 0.0;
    double reference_squared = 0.0;
    for (Index cell : state.mesh->owned_cells) {
        value_squared += value[cell] * value[cell];
        reference_squared += reference[cell] * reference[cell];
    }
    const double local[2] = {value_squared, reference_squared};
    double global[2]{};
    state.parallel.sum(local, global, 2);
    return std::sqrt(global[0]) / std::max(std::sqrt(global[1]), 1e-30);
}

FluxBalance RunTime::fluxBalance(const ScalarField& face_flux) const {
    const Implementation& state = *m_implementation;
    requireFaceField(face_flux, *state.mesh, "face flux");
    FluxBalance result;
    double squared = 0.0;
    double scale = 0.0;
    for (Index cell : state.mesh->owned_cells) {
        double imbalance = 0.0;
        for (Index face : state.mesh->cell_faces[static_cast<std::size_t>(cell)]) {
            const std::size_t index = static_cast<std::size_t>(face);
            const double outward = state.mesh->face_owner[index] == cell
                ? face_flux[face] : -face_flux[face];
            imbalance += outward;
            scale += std::abs(outward);
        }
        result.l1 += std::abs(imbalance);
        squared += imbalance * imbalance;
        result.maximum = std::max(result.maximum, std::abs(imbalance));
    }
    const double local[3] = {result.l1, squared, scale};
    double global[3]{};
    state.parallel.sum(local, global, 3);
    const double local_maximum = result.maximum;
    state.parallel.maximum(&local_maximum, &result.maximum, 1);
    result.l1 = global[0];
    result.l2 = std::sqrt(global[1]);
    result.relative = global[0] / std::max(global[2], 1e-30);
    return result;
}

bool RunTime::all(bool local_condition) const {
    const Implementation& state = *m_implementation;
    return state.parallel.sum(local_condition ? 1 : 0) == state.parallel.size;
}

bool RunTime::loop() {
    Implementation& state = *m_implementation;
    const double next_time = state.current_time + state.control.time.delta_t;
    if (next_time > state.control.time.end_time +
                        0.5 * state.control.time.delta_t) {
        return false;
    }
    state.current_time = next_time;
    ++state.current_step;
    return true;
}

namespace {

void addScalarSource(
    ScalarEquation& equation, const Mesh& mesh, const ScalarFvmTerm& term, int canonical)
{
    const double scale = -static_cast<double>(canonical) * term.coefficient;
    if (term.field == nullptr) {
        for (Index cell : mesh.owned_cells) {
            equation.source[static_cast<std::size_t>(cell)] +=
                scale * mesh.cell_volumes[static_cast<std::size_t>(cell)];
        }
        return;
    }
    requireCellField(*term.field, mesh, "scalar source");
    for (Index cell : mesh.owned_cells) {
        equation.source[static_cast<std::size_t>(cell)] +=
            scale * (*term.field)[cell] * mesh.cell_volumes[static_cast<std::size_t>(cell)];
    }
}

void addVectorSource(
    VectorEquation& equation, const Mesh& mesh, const VectorFvmTerm& term, int canonical)
{
    const double scale = -static_cast<double>(canonical);
    for (Index cell : mesh.owned_cells) {
        equation.source[static_cast<std::size_t>(cell)] +=
            scale * term.vector_source * mesh.cell_volumes[static_cast<std::size_t>(cell)];
    }
}

}  // 匿名命名空间

SolveResult RunTime::solve(const ScalarEquationDefinition& expression) {
    return solve(expression, {});
}

SolveResult RunTime::solve(
    const ScalarEquationDefinition& expression,
    ScalarEquationControl equation_control)
{
    Implementation& state = *m_implementation;
    const ScalarField* unknown_pointer = nullptr;
    const auto inspect_unknown = [&unknown_pointer](const std::vector<ScalarFvmTerm>& terms) {
        for (const ScalarFvmTerm& term : terms) {
            if (term.kind == FvmTermKind::Source || term.field == nullptr) continue;
            if (unknown_pointer != nullptr && unknown_pointer != term.field) {
                throw std::invalid_argument("a scalar equation must have one transported field");
            }
            unknown_pointer = term.field;
        }
    };
    inspect_unknown(expression.lhs.m_terms);
    inspect_unknown(expression.rhs.m_terms);
    if (unknown_pointer == nullptr) {
        throw std::invalid_argument("scalar equation has no unknown field");
    }
    ScalarField& unknown = const_cast<ScalarField&>(*unknown_pointer);
    requireCellField(unknown, *state.mesh, "scalar unknown");
    ScalarEquation equation(*state.mesh);
    bool has_time_derivative = false;

    const auto add = [&](const std::vector<ScalarFvmTerm>& terms, bool lhs) {
        for (const ScalarFvmTerm& term : terms) {
            const int canonical = residualSign(lhs, term.sign);
            switch (term.kind) {
                case FvmTermKind::TimeDerivative: {
                    if (canonical != 1 || term.field != &unknown) {
                        throw std::invalid_argument("ddt must be a positive left-hand-side term");
                    }
                    Implementation::ScalarHistory& history = state.history(unknown);
                    has_time_derivative = true;
                    if (term.coefficient_field != nullptr) {
                        requireCellField(
                            *term.coefficient_field, *state.mesh, "time coefficient");
                        state.synchronize(const_cast<ScalarField&>(*term.coefficient_field));
                        addTimeDerivative(
                            equation, history.previous, deltaT(), *term.coefficient_field,
                            state.control.methods.time,
                            history.has_older ? &history.older : nullptr);
                    } else {
                        addTimeDerivative(
                            equation, history.previous, deltaT(), term.coefficient,
                            state.control.methods.time,
                            history.has_older ? &history.older : nullptr);
                    }
                    break;
                }
                case FvmTermKind::Convection:
                    if (canonical != 1 || term.field != &unknown || term.flux == nullptr) {
                        throw std::invalid_argument("implicit convection must be a positive left-hand-side term");
                    }
                    requireFaceField(*term.flux, *state.mesh, "convection flux");
                    state.synchronize(const_cast<ScalarField&>(*term.flux));
                    state.synchronize(unknown);
                    addConvection(
                        equation, *term.flux, unknown, state.control.methods.convection,
                        state.control.methods.interpolation, state.control.methods.gradient);
                    break;
                case FvmTermKind::Laplacian:
                    if (canonical != -1 || term.field != &unknown) {
                        throw std::invalid_argument("laplacian must appear on the right-hand side or negated on the left");
                    }
                    state.synchronize(unknown);
                    if (term.coefficient_field != nullptr) {
                        const ScalarField& coefficient = *term.coefficient_field;
                        if (&coefficient.mesh() != state.mesh ||
                            (coefficient.location() != FieldLocation::Cell &&
                             coefficient.location() != FieldLocation::Face)) {
                            throw std::invalid_argument("diffusivity must be a cell or face field");
                        }
                        state.synchronize(const_cast<ScalarField&>(coefficient));
                        const ScalarField* face_coefficient = &coefficient;
                        if (coefficient.location() == FieldLocation::Cell) {
                            interpolate(
                                coefficient, state.face_coefficient_workspace,
                                state.control.methods.interpolation,
                                state.control.methods.gradient);
                            state.synchronize(state.face_coefficient_workspace);
                            face_coefficient = &state.face_coefficient_workspace;
                        }
                        addDiffusion(
                            equation, *face_coefficient, unknown,
                            state.control.methods.gradient, state.control.methods.diffusion);
                    } else {
                        addDiffusion(
                            equation, term.coefficient, unknown,
                            state.control.methods.gradient, state.control.methods.diffusion);
                    }
                    break;
                case FvmTermKind::Source:
                    addScalarSource(equation, *state.mesh, term, canonical);
                    break;
                case FvmTermKind::Gradient:
                    throw std::invalid_argument("gradient is not a scalar implicit term");
            }
        }
    };
    add(expression.lhs.m_terms, true);
    add(expression.rhs.m_terms, false);

    if (equation_control.fix_reference) {
        for (Index cell : state.mesh->owned_cells) {
            if (state.mesh->globalCellId(cell) != 0) continue;
            const std::size_t index = static_cast<std::size_t>(cell);
            if (!(equation.diagonal[index] > 0.0)) {
                throw std::runtime_error("reference equation diagonal is invalid");
            }
            equation.diagonal[index] += equation.diagonal[index];
            break;
        }
    }

    state.scalar_assembly.update(equation);
    assembleSource(equation, state.scalar_source);
    prepare(
        state.scalar_assembly.matrix(), equation, state.scalar_pattern_ready,
        state.scalar_linear_solver, state.scalar_distributed_solver.get());
    for (Index cell : state.mesh->owned_cells) {
        state.scalar_solution[state.mesh->ownedIndex(cell)] = unknown[cell];
    }
    const SolveResult result = state.scalar_distributed_solver
        ? state.scalar_distributed_solver->solve(state.scalar_source, state.scalar_solution)
        : state.scalar_linear_solver.solve(state.scalar_source, state.scalar_solution);
    for (Index cell : state.mesh->owned_cells) {
        unknown[cell] = state.scalar_solution[state.mesh->ownedIndex(cell)];
    }
    state.synchronize(unknown);
    if (has_time_derivative) state.history(unknown).advance(unknown);
    return result;
}

std::array<SolveResult, 3> RunTime::solve(const VectorEquationDefinition& expression) {
    return solve(expression, {});
}

std::array<SolveResult, 3> RunTime::solve(
    const VectorEquationDefinition& expression,
    VectorEquationControl equation_control)
{
    Implementation& state = *m_implementation;
    const VectorField* unknown_pointer = nullptr;
    const auto inspect_unknown = [&unknown_pointer](const std::vector<VectorFvmTerm>& terms) {
        for (const VectorFvmTerm& term : terms) {
            if (term.vector_field == nullptr) continue;
            if (unknown_pointer != nullptr && unknown_pointer != term.vector_field) {
                throw std::invalid_argument("a vector equation must have one transported field");
            }
            unknown_pointer = term.vector_field;
        }
    };
    inspect_unknown(expression.lhs.m_terms);
    inspect_unknown(expression.rhs.m_terms);
    if (unknown_pointer == nullptr) {
        throw std::invalid_argument("vector equation has no unknown field");
    }
    VectorField& unknown = const_cast<VectorField&>(*unknown_pointer);
    requireCellField(unknown, *state.mesh, "vector unknown");
    VectorEquation equation(*state.mesh);
    bool has_time_derivative = false;

    const auto add = [&](const std::vector<VectorFvmTerm>& terms, bool lhs) {
        for (const VectorFvmTerm& term : terms) {
            const int canonical = residualSign(lhs, term.sign);
            switch (term.kind) {
                case FvmTermKind::TimeDerivative: {
                    if (canonical != 1 || term.vector_field != &unknown) {
                        throw std::invalid_argument("ddt must be a positive left-hand-side term");
                    }
                    Implementation::VectorHistory& history = state.history(unknown);
                    has_time_derivative = true;
                    if (term.coefficient_field != nullptr) {
                        requireCellField(
                            *term.coefficient_field, *state.mesh, "time coefficient");
                        state.synchronize(const_cast<ScalarField&>(*term.coefficient_field));
                        addTimeDerivative(
                            equation, history.previous, deltaT(), *term.coefficient_field,
                            state.control.methods.time,
                            history.has_older ? &history.older : nullptr);
                    } else {
                        addTimeDerivative(
                            equation, history.previous, deltaT(), term.coefficient,
                            state.control.methods.time,
                            history.has_older ? &history.older : nullptr);
                    }
                    break;
                }
                case FvmTermKind::Convection:
                    if (canonical != 1 || term.vector_field != &unknown || term.flux == nullptr) {
                        throw std::invalid_argument("implicit convection must be a positive left-hand-side term");
                    }
                    requireFaceField(*term.flux, *state.mesh, "convection flux");
                    state.synchronize(const_cast<ScalarField&>(*term.flux));
                    state.synchronize(unknown);
                    addConvection(
                        equation, *term.flux, unknown, state.control.methods.convection,
                        state.control.methods.interpolation, state.control.methods.gradient);
                    break;
                case FvmTermKind::Laplacian:
                    if (canonical != -1 || term.vector_field != &unknown) {
                        throw std::invalid_argument("laplacian must appear on the right-hand side or negated on the left");
                    }
                    state.synchronize(unknown);
                    if (term.coefficient_field != nullptr) {
                        const ScalarField& coefficient = *term.coefficient_field;
                        if (&coefficient.mesh() != state.mesh ||
                            (coefficient.location() != FieldLocation::Cell &&
                             coefficient.location() != FieldLocation::Face)) {
                            throw std::invalid_argument("diffusivity must be a cell or face field");
                        }
                        state.synchronize(const_cast<ScalarField&>(coefficient));
                        const ScalarField* face_coefficient = &coefficient;
                        if (coefficient.location() == FieldLocation::Cell) {
                            interpolate(
                                coefficient, state.face_coefficient_workspace,
                                state.control.methods.interpolation,
                                state.control.methods.gradient);
                            state.synchronize(state.face_coefficient_workspace);
                            face_coefficient = &state.face_coefficient_workspace;
                        }
                        addDiffusion(
                            equation, *face_coefficient, unknown,
                            state.control.methods.gradient, state.control.methods.diffusion);
                    } else {
                        addDiffusion(
                            equation, term.coefficient, unknown,
                            state.control.methods.gradient, state.control.methods.diffusion);
                    }
                    break;
                case FvmTermKind::Gradient:
                    if (term.scalar_field == nullptr) {
                        throw std::invalid_argument("gradient term has no scalar field");
                    }
                    requireCellField(*term.scalar_field, *state.mesh, "gradient field");
                    state.synchronize(const_cast<ScalarField&>(*term.scalar_field));
                    gradient(
                        *term.scalar_field, state.gradient_workspace,
                        state.control.methods.gradient);
                    state.synchronize(state.gradient_workspace);
                    for (Index cell : state.mesh->owned_cells) {
                        equation.source[static_cast<std::size_t>(cell)] -=
                            static_cast<double>(canonical) *
                            state.mesh->cell_volumes[static_cast<std::size_t>(cell)] *
                            state.gradient_workspace[cell];
                    }
                    break;
                case FvmTermKind::Source:
                    addVectorSource(equation, *state.mesh, term, canonical);
                    break;
            }
        }
    };
    add(expression.lhs.m_terms, true);
    add(expression.rhs.m_terms, false);

    if (!(equation_control.relaxation > 0.0 &&
          equation_control.relaxation <= 1.0)) {
        throw std::invalid_argument("vector equation relaxation must be in (0, 1]");
    }
    if (equation_control.mobility != nullptr) {
        requireCellField(*equation_control.mobility, *state.mesh, "equation mobility");
    }
    if (equation_control.relaxation != 1.0) {
        for (Index face : state.mesh->owned_faces) {
            equation.upper[static_cast<std::size_t>(face)] *= equation_control.relaxation;
            equation.lower[static_cast<std::size_t>(face)] *= equation_control.relaxation;
        }
    }
    for (Index cell : state.mesh->owned_cells) {
        const std::size_t index = static_cast<std::size_t>(cell);
        if (!(equation.diagonal[index] > 0.0) || !std::isfinite(equation.diagonal[index])) {
            throw std::runtime_error("vector equation diagonal is not positive and finite");
        }
        if (equation_control.relaxation != 1.0) {
            equation.source[index] = equation_control.relaxation * equation.source[index] +
                (1.0 - equation_control.relaxation) *
                equation.diagonal[index] * unknown[cell];
        }
        if (equation_control.mobility != nullptr) {
            (*equation_control.mobility)[cell] =
                state.mesh->cell_volumes[index] / equation.diagonal[index];
        }
    }
    if (equation_control.mobility != nullptr) {
        state.synchronize(*equation_control.mobility);
    }

    state.vector_assembly.update(equation);
    assembleSource(equation, state.vector_source);
    prepare(
        state.vector_assembly.matrix(), equation, state.vector_pattern_ready,
        state.vector_linear_solver, state.vector_distributed_solver.get());
    for (std::size_t component = 0; component < 3; ++component) {
        for (Index cell : state.mesh->owned_cells) {
            state.vector_solution[component][state.mesh->ownedIndex(cell)] =
                unknown[cell][component];
        }
    }
    std::array<SolveResult, 3> results;
    for (std::size_t component = 0; component < 3; ++component) {
        results[component] = state.vector_distributed_solver
            ? state.vector_distributed_solver->solve(
                  state.vector_source[component], state.vector_solution[component])
            : state.vector_linear_solver.solve(
                  state.vector_source[component], state.vector_solution[component]);
    }
    for (Index cell : state.mesh->owned_cells) {
        for (std::size_t component = 0; component < 3; ++component) {
            unknown[cell][component] =
                state.vector_solution[component][state.mesh->ownedIndex(cell)];
        }
    }
    state.synchronize(unknown);
    if (has_time_derivative) state.history(unknown).advance(unknown);
    return results;
}

void RunTime::evaluate(fvc::ScalarGradient operation, VectorField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "gradient input");
    requireCellField(result, *state.mesh, "gradient result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    gradient(operation.field, result, state.control.methods.gradient);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::VectorGradient operation, TensorField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "gradient input");
    requireCellField(result, *state.mesh, "gradient result");
    state.synchronize(const_cast<VectorField&>(operation.field));
    gradient(operation.field, result, state.control.methods.gradient);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::FaceFlux operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.velocity, *state.mesh, "flux input");
    requireFaceField(result, *state.mesh, "flux result");
    state.synchronize(const_cast<VectorField&>(operation.velocity));
    flux(operation.velocity, result, state.control.methods.interpolation, state.control.methods.gradient);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::FaceDivergence operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireFaceField(operation.flux, *state.mesh, "divergence input");
    requireCellField(result, *state.mesh, "divergence result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    divergence(operation.flux, result);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::VectorDivergence operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "divergence input");
    requireCellField(result, *state.mesh, "divergence result");
    state.synchronize(const_cast<VectorField&>(operation.field));
    divergence(
        operation.field, result, state.control.methods.interpolation,
        state.control.methods.gradient);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::ScalarConvection operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireFaceField(operation.flux, *state.mesh, "convection flux");
    requireCellField(operation.field, *state.mesh, "convection field");
    requireCellField(result, *state.mesh, "convection result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    state.synchronize(const_cast<ScalarField&>(operation.field));
    convection(
        operation.flux, operation.field, result, state.control.methods.convection,
        state.control.methods.interpolation, state.control.methods.gradient);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::VectorConvection operation, VectorField& result) {
    Implementation& state = *m_implementation;
    requireFaceField(operation.flux, *state.mesh, "convection flux");
    requireCellField(operation.field, *state.mesh, "convection field");
    requireCellField(result, *state.mesh, "convection result");
    state.synchronize(const_cast<ScalarField&>(operation.flux));
    state.synchronize(const_cast<VectorField&>(operation.field));
    convection(
        operation.flux, operation.field, result, state.control.methods.convection,
        state.control.methods.interpolation, state.control.methods.gradient);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::ScalarInterpolation operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "interpolation input");
    requireFaceField(result, *state.mesh, "interpolation result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    interpolate(
        operation.field, result, state.control.methods.interpolation,
        state.control.methods.gradient);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::VectorInterpolation operation, VectorField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "interpolation input");
    if (&result.mesh() != state.mesh || result.location() != FieldLocation::Face) {
        throw std::invalid_argument("vector interpolation result must be a face field");
    }
    state.synchronize(const_cast<VectorField&>(operation.field));
    interpolate(
        operation.field, result, state.control.methods.interpolation,
        state.control.methods.gradient);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::ScalarReconstruction operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "reconstruction input");
    requireCellField(operation.gradient, *state.mesh, "reconstruction gradient");
    requireFaceField(result, *state.mesh, "reconstruction result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    state.synchronize(const_cast<VectorField&>(operation.gradient));
    reconstruct(operation.field, operation.gradient, result);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::VectorReconstruction operation, VectorField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "reconstruction input");
    requireCellField(operation.gradient, *state.mesh, "reconstruction gradient");
    if (&result.mesh() != state.mesh || result.location() != FieldLocation::Face) {
        throw std::invalid_argument("vector reconstruction result must be a face field");
    }
    state.synchronize(const_cast<VectorField&>(operation.field));
    state.synchronize(const_cast<TensorField&>(operation.gradient));
    reconstruct(operation.field, operation.gradient, result);
    state.synchronize(result);
}

void RunTime::evaluate(fvc::ScalarLaplacian operation, ScalarField& result) {
    Implementation& state = *m_implementation;
    requireCellField(operation.field, *state.mesh, "laplacian input");
    requireCellField(result, *state.mesh, "laplacian result");
    state.synchronize(const_cast<ScalarField&>(operation.field));
    if (operation.coefficient_field == nullptr) {
        laplacian(
            operation.coefficient, operation.field, result,
            state.control.methods.gradient, state.control.methods.diffusion);
    } else {
        const ScalarField& coefficient = *operation.coefficient_field;
        if (&coefficient.mesh() != state.mesh ||
            (coefficient.location() != FieldLocation::Cell &&
             coefficient.location() != FieldLocation::Face)) {
            throw std::invalid_argument("laplacian diffusivity must be a cell or face field");
        }
        state.synchronize(const_cast<ScalarField&>(coefficient));
        if (coefficient.location() == FieldLocation::Cell) {
            interpolate(
                coefficient, state.face_coefficient_workspace,
                state.control.methods.interpolation,
                state.control.methods.gradient);
            state.synchronize(state.face_coefficient_workspace);
            laplacian(
                state.face_coefficient_workspace, operation.field, result,
                state.control.methods.gradient, state.control.methods.diffusion);
        } else {
            laplacian(
                coefficient, operation.field, result, state.control.methods.gradient,
                state.control.methods.diffusion);
        }
    }
    state.synchronize(result);
}

SolveResult solve(const ScalarEquationDefinition& equation) {
    return RunTime::current().solve(equation);
}

std::array<SolveResult, 3> solve(const VectorEquationDefinition& equation) {
    return RunTime::current().solve(equation);
}

namespace fvc {

void evaluate(ScalarGradient operation, VectorField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(VectorGradient operation, TensorField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(FaceFlux operation, ScalarField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(FaceDivergence operation, ScalarField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(VectorDivergence operation, ScalarField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(ScalarConvection operation, ScalarField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(VectorConvection operation, VectorField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(ScalarInterpolation operation, ScalarField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(VectorInterpolation operation, VectorField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(ScalarReconstruction operation, ScalarField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(VectorReconstruction operation, VectorField& result) {
    RunTime::current().evaluate(operation, result);
}
void evaluate(ScalarLaplacian operation, ScalarField& result) {
    RunTime::current().evaluate(operation, result);
}

}  // fvc 命名空间

namespace diagnostics {

double relativeChange(const VectorField& current, const VectorField& previous) {
    return RunTime::current().relativeChange(current, previous);
}

double relativeChange(const ScalarField& current, const ScalarField& previous) {
    return RunTime::current().relativeChange(current, previous);
}

double relativeMagnitude(const ScalarField& value, const ScalarField& reference) {
    return RunTime::current().relativeMagnitude(value, reference);
}

FluxBalance fluxBalance(const ScalarField& face_flux) {
    return RunTime::current().fluxBalance(face_flux);
}

bool all(bool local_condition) {
    return RunTime::current().all(local_condition);
}

}  // diagnostics 命名空间

}  // babelsim 命名空间
