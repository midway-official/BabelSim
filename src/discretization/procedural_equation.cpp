#include "babelsim/equ.h"
#include "babelsim/operators.h"
#include "internal/fvm_execution.h"
#include "internal/compute_backend.h"
#include "internal/field_access.h"
#include "internal/mesh_access.h"
#include "operator_support.h"
#include <cmath>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace babelsim::equ {
namespace impl = babelsim::detail;
template <class T>
struct Equation<T>::Storage {
    explicit Storage(Field<T>& field) : unknown(&field), coefficients(field.mesh()) {}
    Field<T>* unknown;
    DiscreteEquation<T> coefficients;
    EquationControl control;
    // 本方程已绑定的唯一面通量上下文。不同通量会给出不同的边界迹，
    // 一个方程只能用一个；记录首个指针以便不一致时立即报错。
    const ScalarField* boundary_flux = nullptr;
    struct DiffusionSnapshot {
        ScalarField coefficient;
        ScalarField iterate;
        std::optional<VectorField> gradient;
        DiffusionMethod method;
    };
    struct DiffusionContribution {
        std::shared_ptr<const DiffusionSnapshot> snapshot;
        double factor;
    };
    std::vector<DiffusionContribution> diffusion;
};
namespace detail {
struct EquationAccess {
    template <class T>
    static auto& get(Equation<T>& equation) {
        return *equation.storage_;
    }
    template <class T>
    static const auto& get(const Equation<T>& equation) {
        return *equation.storage_;
    }
};
}
using detail::EquationAccess;
namespace {
auto& backend() {
    return impl::execution().backend();
}
template <class T>
void requireCellField(const Field<T>& f, const Mesh& mesh) {
    f.validateStorage();
    if (&f.mesh() != &mesh || f.location() != FieldLocation::Cell)
        throw std::invalid_argument("equ requires compatible cell fields");
}
template <class T>
auto& equationState(const Equation<T>& equation) {
    auto& equation_state = EquationAccess::get(equation);
    if (&equation_state.unknown->mesh() != &impl::execution().mesh())
        throw std::invalid_argument("equ matrix belongs to a different run mesh");
    return equation_state;
}
template <class T>
auto& equationState(Equation<T>& equation) {
    equationState(static_cast<const Equation<T>&>(equation));
    return EquationAccess::get(equation);
}
void requireFinite(double v) {
    if (!std::isfinite(v))
        throw std::invalid_argument("nonfinite equ coefficient");
}
template <class T>
void accumulate(DiscreteEquation<T>& equation, const DiscreteEquation<T>& b, double f) {
    requireFinite(f);
    for (std::size_t i = 0; i < equation.diagonal.size(); ++i) {
        equation.diagonal[i] += f * b.diagonal[i];
        equation.source[i] += f * b.source[i];
    }
    for (std::size_t i = 0; i < equation.upper.size(); ++i) {
        equation.upper[i] += f * b.upper[i];
        equation.lower[i] += f * b.lower[i];
    }
}
template <class T>
void synchronize(Field<T>& f) {
    backend().synchronize(f);
}
template <class T>
Field<T> resultField(const Equation<T>& equation, const char* name) {
    return Field<T>(equationState(equation).unknown->mesh(), FieldLocation::Cell, name);
}
// Freeze the actual explicit flux, including any limiter evaluated at assembly.
// A frozen gradient alone is insufficient when the limiter also depends on x.
ScalarField diffusionCorrection(const ScalarField& diffusivity, const ScalarField& field,
                                const VectorField* reconstructed_gradient,
                                DiffusionMethod diffusionMethod) {
    ScalarField correction(field.mesh(), FieldLocation::Face, "equ.frozenCorrection");
    if (diffusionMethod == DiffusionMethod::Orthogonal)
        return correction;
    if (!reconstructed_gradient)
        throw std::logic_error("corrected diffusion requires a reconstructed gradient");
    diffusionFlux(diffusivity, field, *reconstructed_gradient, correction, diffusionMethod);
    ScalarField orthogonal(field.mesh(), FieldLocation::Face, "equ.orthogonalFlux");
    diffusionFlux(diffusivity, field, *reconstructed_gradient, orthogonal,
                  DiffusionMethod::Orthogonal);
    correction -= orthogonal;
    return correction;
}
}
template <class T>
Equation<T>::Equation(Field<T>& field, const EquationControl& control)
    : storage_(std::make_unique<Storage>(field)) {
    requireCellField(field, impl::execution().mesh());
    if (field.calculatedBoundary())
        throw std::invalid_argument("equ unknown requires physical boundaries");
    control.linear.validate();
    control.requireSpatial(control.name);
    storage_->control = control;
    for (auto& entry : storage_->control.terms) {
        entry.second.requireComplete(control.name + ".term." + entry.first);
    }
    storage_->coefficients.numerical_identity = control.name;
}
template <class T>
OperatorOptions Equation<T>::options(const std::string& term) const {
    return storage_->control.options(term);
}
template <class T>
const std::string& Equation<T>::name() const {
    return storage_->control.name;
}
template <class T>
Equation<T>::~Equation() = default;
template <class T>
Equation<T>::Equation(Equation&&) noexcept = default;
template <class T>
Equation<T>& Equation<T>::operator=(Equation&&) noexcept = default;
template <class T>
void clear(Equation<T>& equation) {
    auto& storage = equationState(equation);
    storage.coefficients.reset();
    storage.diffusion.clear();
    storage.boundary_flux = nullptr;
}
template <class T>
Equation<T> copy(const Equation<T>& equation) {
    const auto& original = equationState(equation);
    Equation<T> result(*original.unknown, original.control);
    auto& copied = equationState(result);
    copied.coefficients = original.coefficients;
    copied.diffusion = original.diffusion;
    copied.boundary_flux = original.boundary_flux;
    return result;
}
// Spatial terms: assemble immediately into the bound equation.
template <class T>
void div(Equation<T>& equation, const ScalarField& phi, double factor, const std::string& term,
         const OperatorOptions& overrides) {
    requireFinite(factor);
    auto& equation_state = equationState(equation);
    auto& field = *equation_state.unknown;
    if (&phi.mesh() != &field.mesh() || phi.location() != FieldLocation::Face)
        throw std::invalid_argument("equ::div requires a face flux");
    if (equation_state.boundary_flux != nullptr && equation_state.boundary_flux != &phi)
        throw std::invalid_argument("one equation requires one boundary flux context");
    equation_state.boundary_flux = &phi;
    backend().synchronize(const_cast<ScalarField&>(phi));
    auto options = equation.options(term);
    options.overlay(overrides);
    // Upwind reads only the face flux and the local owner at physical boundaries.
    if (*options.convection != ConvectionMethod::Upwind)
        synchronize(field);
    addConvection(equation_state.coefficients, phi, field, *options.convection,
                  *options.interpolation, *options.gradient, factor);
}
template <class T>
void laplacian(Equation<T>& equation, double diffusivity, double factor, const std::string& term,
               const OperatorOptions& overrides) {
    requireFinite(diffusivity);
    requireFinite(factor);
    auto& equation_state = equationState(equation);
    auto& field = *equation_state.unknown;
    synchronize(field);
    DiscreteEquation<T> contribution(field.mesh());
    auto options = equation.options(term);
    options.overlay(overrides);
    if constexpr (std::is_same_v<T, double>) {
        std::optional<VectorField> reconstructed_gradient;
        if (*options.diffusion != DiffusionMethod::Orthogonal) {
            reconstructed_gradient.emplace(field.mesh(), FieldLocation::Cell, "equ.frozenGradient");
            impl::operators::gradientForExchange(field, *reconstructed_gradient, *options.gradient);
            backend().synchronizeUpdated(*reconstructed_gradient);
            addDiffusion(contribution, diffusivity, field, *reconstructed_gradient,
                         *options.diffusion);
        } else
            addDiffusion(contribution, diffusivity, field, *options.gradient, *options.diffusion);
        ScalarField face_diffusivity(field.mesh(), FieldLocation::Face, "equ.frozenDiffusivity",
                                     diffusivity);
        using Snapshot =
            typename std::remove_reference_t<decltype(equation_state)>::DiffusionSnapshot;
        equation_state.diffusion.push_back(
            {std::make_shared<Snapshot>(
                 Snapshot{std::move(face_diffusivity), field, std::move(reconstructed_gradient),
                          *options.diffusion}),
             factor});
    } else
        addDiffusion(contribution, diffusivity, field, *options.gradient, *options.diffusion);
    // addDiffusion assembles -div(k grad(field)); equ::laplacian represents
    // +div(k grad(field)), hence the negative multiplier here.
    accumulate(equation_state.coefficients, contribution, -factor);
}
template <class T>
void laplacian(Equation<T>& equation, const ScalarField& diffusivity, double factor,
               const std::string& term, const OperatorOptions& overrides) {
    requireFinite(factor);
    auto& equation_state = equationState(equation);
    auto& field = *equation_state.unknown;
    if (&diffusivity.mesh() != &field.mesh() || (diffusivity.location() != FieldLocation::Cell &&
                                                 diffusivity.location() != FieldLocation::Face))
        throw std::invalid_argument("equ diffusivity must be a compatible cell or face field");
    backend().synchronize(const_cast<ScalarField&>(diffusivity));
    synchronize(field);
    std::optional<ScalarField> face_diffusivity;
    auto options = equation.options(term);
    options.overlay(overrides);
    const ScalarField* face = &diffusivity;
    if (diffusivity.location() == FieldLocation::Cell) {
        face_diffusivity.emplace(field.mesh(), FieldLocation::Face, "equ.diffusivity");
        interpolate(diffusivity, *face_diffusivity,
                    options.coefficientInterpolation.value_or(*options.interpolation),
                    options.coefficientGradient.value_or(*options.gradient));
        backend().synchronizeUpdated(*face_diffusivity);
        face = &*face_diffusivity;
    }
    DiscreteEquation<T> contribution(field.mesh());
    if constexpr (std::is_same_v<T, double>) {
        std::optional<VectorField> reconstructed_gradient;
        if (*options.diffusion != DiffusionMethod::Orthogonal) {
            // Flux recovery can read every local face, including ghost faces.
            // Validate the frozen coefficients here, while assembling the term.
            for (double value : impl::fieldValues(*face))
                if (value < 0.0 || !std::isfinite(value))
                    throw std::invalid_argument("face diffusivity must be non-negative and finite");
            reconstructed_gradient.emplace(field.mesh(), FieldLocation::Cell, "equ.frozenGradient");
            impl::operators::gradientForExchange(field, *reconstructed_gradient, *options.gradient);
            backend().synchronizeUpdated(*reconstructed_gradient);
            addDiffusion(contribution, *face, field, *reconstructed_gradient, *options.diffusion);
        } else
            addDiffusion(contribution, *face, field, *options.gradient, *options.diffusion);
        using Snapshot =
            typename std::remove_reference_t<decltype(equation_state)>::DiffusionSnapshot;
        equation_state.diffusion.push_back(
            {std::make_shared<Snapshot>(
                 Snapshot{*face, field, std::move(reconstructed_gradient), *options.diffusion}),
             factor});
    } else
        addDiffusion(contribution, *face, field, *options.gradient, *options.diffusion);
    // addDiffusion assembles -div(k grad(field)); equ::laplacian represents
    // +div(k grad(field)), hence the negative multiplier here.
    accumulate(equation_state.coefficients, contribution, -factor);
}
template <class T>
void reaction(Equation<T>& equation, double c) {
    requireFinite(c);
    auto& equation_state = equationState(equation);
    const auto& mesh = equation_state.unknown->mesh();
    for (Index i : impl::meshData(mesh).owned_cells)
        equation_state.coefficients.diagonal[i] += c * mesh.cellVolume(i);
}
template <class T>
void reaction(Equation<T>& equation, const ScalarField& c, double factor) {
    requireFinite(factor);
    auto& equation_state = equationState(equation);
    const auto& mesh = equation_state.unknown->mesh();
    requireCellField(c, mesh);
    for (Index i : impl::meshData(mesh).owned_cells) {
        requireFinite(impl::fieldData(c)[i]);
        equation_state.coefficients.diagonal[i] +=
            factor * impl::fieldData(c)[i] * mesh.cellVolume(i);
    }
}
template <class T>
void source(Equation<T>& equation, T value) {
    auto& equation_state = equationState(equation);
    const auto& mesh = equation_state.unknown->mesh();
    for (Index i : impl::meshData(mesh).owned_cells)
        equation_state.coefficients.source[i] += mesh.cellVolume(i) * value;
}
template <class T>
void source(Equation<T>& equation, const Field<T>& value, double factor) {
    requireFinite(factor);
    auto& equation_state = equationState(equation);
    const auto& mesh = equation_state.unknown->mesh();
    requireCellField(value, mesh);
    for (Index i : impl::meshData(mesh).owned_cells)
        equation_state.coefficients.source[i] +=
            factor * mesh.cellVolume(i) * impl::fieldData(value)[i];
}
namespace {
// Time terms: history is read here; advancing history belongs to the caller.
template <class T, class Capacity>
void timeTerm(Equation<T>& equation, const Capacity& capacity, const Field<T>& previous, double dt,
              TimeMethod method, const Field<T>* older, double previous_dt) {
    auto& equation_state = equationState(equation);
    const auto& mesh = equation_state.unknown->mesh();
    requireCellField(previous, mesh);
    if (method == TimeMethod::Steady)
        return;
    if (!(dt > 0) || !std::isfinite(dt))
        throw std::invalid_argument("equ::ddt requires positive finite dt");
    if constexpr (!std::is_same_v<Capacity, double>)
        requireCellField(capacity, mesh);
    // V * capacity * (a0 * current - a1 * previous - a2 * older) = 0.
    // a0 enters the matrix; the two known time levels enter the right-hand side.
    double a0 = 1 / dt, a1 = 1 / dt, a2 = 0;
    if (method == TimeMethod::BDF2) {
        if (!older)
            throw std::invalid_argument("BDF2 requires explicit older field; use Euler to start");
        requireCellField(*older, mesh);
        if (previous_dt == 0)
            previous_dt = dt;
        if (!(previous_dt > 0) || !std::isfinite(previous_dt))
            throw std::invalid_argument("invalid previous dt");
        const double r = dt / previous_dt;
        a0 = (1 + 2 * r) / ((1 + r) * dt);
        a1 = (1 + r) / dt;
        a2 = -r * r / ((1 + r) * dt);
    }
    for (Index i : impl::meshData(mesh).owned_cells) {
        double c;
        if constexpr (std::is_same_v<Capacity, double>)
            c = capacity;
        else
            c = impl::fieldData(capacity)[i];
        if (!(c > 0) || !std::isfinite(c))
            throw std::invalid_argument("time capacity must be positive and finite");
        const double w = c * mesh.cellVolume(i);
        equation_state.coefficients.diagonal[i] += w * a0;
        equation_state.coefficients.source[i] += w * a1 * impl::fieldData(previous)[i];
        if (a2 != 0)
            equation_state.coefficients.source[i] += w * a2 * impl::fieldData(*older)[i];
    }
}
}
template <class T>
void ddt(Equation<T>& equation, double c, const Field<T>& old, double dt, TimeMethod options,
         const Field<T>* older, double previous_dt) {
    timeTerm(equation, c, old, dt, options, older, previous_dt);
}
template <class T>
void ddt(Equation<T>& equation, const ScalarField& c, const Field<T>& old, double dt,
         TimeMethod options, const Field<T>* older, double previous_dt) {
    timeTerm(equation, c, old, dt, options, older, previous_dt);
}
template <class T, class Capacity>
void historyTime(Equation<T>& equation, const Capacity& c, const time::History<T>& h) {
    if (numericalMethods().time == TimeMethod::Steady)
        return;
    if (h.levels() == 0)
        throw std::logic_error("saveOld must precede a transient assembly");
    const auto configured = numericalMethods().time;
    const auto method =
        configured == TimeMethod::BDF2 && h.levels() < 2 ? TimeMethod::Euler : configured;
    ddt(equation, c, h.previous(), h.dt(), method, &h.older(), h.previousDt());
}
template <class T>
void ddt(Equation<T>& equation, double c, const time::History<T>& h) {
    historyTime(equation, c, h);
}
template <class T>
void ddt(Equation<T>& equation, const ScalarField& c, const time::History<T>& h) {
    historyTime(equation, c, h);
}
// Algebra operates on integrated coefficients (no extra cell-volume factor).
template <class T>
void scale(Equation<T>& equation, double f) {
    requireFinite(f);
    auto& coefficients = equationState(equation).coefficients;
    for (auto& v : coefficients.diagonal)
        v *= f;
    for (auto& v : coefficients.upper)
        v *= f;
    for (auto& v : coefficients.lower)
        v *= f;
    for (auto& v : coefficients.source)
        v = f * v;
    for (auto& contribution : equationState(equation).diffusion)
        contribution.factor *= f;
}
template <class T>
void add(Equation<T>& equation, const Equation<T>& b, double f) {
    auto& target_state = equationState(equation);
    const auto& source_state = equationState(b);
    if (target_state.unknown != source_state.unknown)
        throw std::invalid_argument("matrix addition requires same unknown binding");
    if (source_state.boundary_flux != nullptr) {
        if (target_state.boundary_flux != nullptr &&
            target_state.boundary_flux != source_state.boundary_flux)
            throw std::invalid_argument("one equation requires one boundary flux context");
        target_state.boundary_flux = source_state.boundary_flux;
    }
    // Copy first so self-addition does not invalidate iterators.
    auto contributions = source_state.diffusion;
    accumulate(target_state.coefficients, source_state.coefficients, f);
    for (auto& contribution : contributions) {
        contribution.factor *= f;
        target_state.diffusion.push_back(std::move(contribution));
    }
}
template <class T>
void addDiagonal(Equation<T>& equation, const ScalarField& d) {
    auto& equation_state = equationState(equation);
    requireCellField(d, equation_state.unknown->mesh());
    for (Index i : impl::meshData(d.mesh()).owned_cells)
        equation_state.coefficients.diagonal[i] += impl::fieldData(d)[i];
}
template <class T>
void addRhs(Equation<T>& equation, const Field<T>& b) {
    auto& equation_state = equationState(equation);
    requireCellField(b, equation_state.unknown->mesh());
    for (Index i : impl::meshData(b.mesh()).owned_cells)
        equation_state.coefficients.source[i] += impl::fieldData(b)[i];
}
template <class T>
ScalarField diagonal(const Equation<T>& equation) {
    const auto& equation_state = equationState(equation);
    ScalarField d(equation_state.unknown->mesh(), FieldLocation::Cell, "equ.diagonal");
    for (Index i : impl::meshData(d.mesh()).owned_cells)
        impl::fieldData(d)[i] = equation_state.coefficients.diagonal[i];
    backend().synchronizeUpdated(d);
    return d;
}
template <class T>
Field<T> rhs(const Equation<T>& equation) {
    auto b = resultField(equation, "equ.rhs");
    const auto& coefficients = equationState(equation).coefficients;
    for (Index i : impl::meshData(b.mesh()).owned_cells)
        impl::fieldData(b)[i] = coefficients.source[i];
    backend().synchronizeUpdated(b);
    return b;
}
namespace {
template <class T>
Field<T> applyOwned(const Equation<T>& equation, const Field<T>& field) {
    const auto& equation_state = equationState(equation);
    const auto& mesh = equation_state.unknown->mesh();
    requireCellField(field, mesh);
    backend().synchronize(const_cast<Field<T>&>(field));
    auto y = resultField(equation, "equ.apply");
    auto* values = impl::fieldData(y);
    const auto* xv = impl::fieldData(field);
    const auto& coefficients = equation_state.coefficients;
    const auto& geometry = impl::meshData(mesh);
    for (Index i : geometry.owned_cells)
        values[i] = coefficients.diagonal[i] * xv[i];
    // upper couples owner <- neighbour; lower couples neighbour <- owner.
    // Ghost values are readable, but only owned rows are written.
    for (Index f : geometry.owned_faces) {
        const Index o = geometry.face_owner[f], n = geometry.face_neighbour[f];
        if (n == invalid_index)
            continue;
        if (impl::isOwned(mesh, o))
            values[o] += coefficients.upper[f] * xv[n];
        if (impl::isOwned(mesh, n))
            values[n] += coefficients.lower[f] * xv[o];
    }
    return y;
}
} // namespace

template <class T>
Field<T> apply(const Equation<T>& equation, const Field<T>& field) {
    auto y = applyOwned(equation, field);
    backend().synchronizeUpdated(y);
    return y;
}
template <class T>
Field<T> residual(const Equation<T>& equation, const Field<T>& field) {
    auto r = applyOwned(equation, field);
    const auto& coefficients = equationState(equation).coefficients;
    auto* values = impl::fieldData(r);
    for (Index i : impl::meshData(field.mesh()).owned_cells)
        values[i] = coefficients.source[i] - values[i];
    backend().synchronizeUpdated(r);
    return r;
}

ScalarField offDiagonalSum(const Equation<Vec3>& equation) {
    const auto& state = equationState(equation);
    const auto& mesh = state.unknown->mesh();
    const auto& coefficients = state.coefficients;
    ScalarField result(mesh, FieldLocation::Cell, "equ.offDiagonalSum");
    auto* values = impl::fieldData(result);
    for (Index face = 0; face < mesh.faceCount(); ++face) {
        const Index owner = mesh.owner(face), neighbour = mesh.neighbour(face);
        if (neighbour == invalid_index)
            continue;
        if (impl::isOwned(mesh, owner))
            values[owner] += coefficients.upper[face];
        if (impl::isOwned(mesh, neighbour))
            values[neighbour] += coefficients.lower[face];
    }
    synchronize(result);
    return result;
}

template <class T>
void relax(Equation<T>& equation, const Field<T>& previous, double alpha) {
    if (!(alpha > 0 && alpha <= 1))
        throw std::invalid_argument("relaxation must be in (0,1]");
    auto& equation_state = equationState(equation);
    requireCellField(previous, equation_state.unknown->mesh());
    auto& coefficients = equation_state.coefficients;
    for (Index i : impl::meshData(previous.mesh()).owned_cells) {
        const double extra = coefficients.diagonal[i] * (1 / alpha - 1);
        coefficients.diagonal[i] += extra;
        coefficients.source[i] += extra * impl::fieldData(previous)[i];
    }
}
namespace {
// Reference constraints use global cell IDs; every rank participates.
Index applyReference(ScalarDiscreteEquation& coefficients, Index global_cell, double value) {
    requireFinite(value);
    const auto& mesh = *coefficients.mesh;
    Index localCell = invalid_index;
    for (Index i : impl::meshData(mesh).owned_cells) {
        if (impl::globalCellId(mesh, i) != global_cell)
            continue;
        if (!(coefficients.diagonal[i] > 0))
            throw std::invalid_argument("reference requires a positive diagonal");
        coefficients.source[i] += coefficients.diagonal[i] * value;
        coefficients.diagonal[i] *= 2;
        localCell = i;
    }
    if (backend().all(localCell == invalid_index))
        throw std::invalid_argument("reference global cell does not exist");
    return localCell;
}

bool hasConstantNullMode(const ScalarDiscreteEquation& coefficients) {
    const auto& mesh = *coefficients.mesh;
    auto rowSum = coefficients.diagonal;
    auto magnitude = coefficients.diagonal;
    for (auto& v : magnitude)
        v = std::abs(v);
    for (Index f = 0; f < mesh.faceCount(); ++f) {
        const Index o = mesh.owner(f), n = mesh.neighbour(f);
        if (n == invalid_index)
            continue;
        if (impl::isOwned(mesh, o)) {
            rowSum[o] += coefficients.upper[f];
            magnitude[o] += std::abs(coefficients.upper[f]);
        }
        if (impl::isOwned(mesh, n)) {
            rowSum[n] += coefficients.lower[f];
            magnitude[n] += std::abs(coefficients.lower[f]);
        }
    }
    bool constantMode = true;
    for (Index i : impl::meshData(mesh).owned_cells)
        constantMode =
            constantMode && std::isfinite(rowSum[i]) &&
            std::abs(rowSum[i]) <= 64 * std::numeric_limits<double>::epsilon() * magnitude[i];
    return backend().all(constantMode);
}

ScalarField frozenDiffusionFlux(const ScalarField& coefficient, const ScalarField& iterate,
                                const ScalarField& correction, const VectorField& unusedGradient) {
    ScalarField flux(iterate.mesh(), FieldLocation::Face, "equ.diffusionFluxTerm");
    diffusionFlux(coefficient, iterate, unusedGradient, flux, DiffusionMethod::Orthogonal);
    flux += correction;
    return flux;
}
} // namespace

void reference(Equation<double>& equation, Index global_cell, double value) {
    applyReference(equationState(equation).coefficients, global_cell, value);
}

void reference(Equation<double>& equation, double value) {
    requireFinite(value);
    auto& coefficients = equationState(equation).coefficients;
    if (hasConstantNullMode(coefficients))
        applyReference(coefficients, 0, value);
}

ScalarField faceFlux(const Equation<double>& equation, const ScalarField& solution) {
    const auto& equation_state = equationState(equation);
    requireCellField(solution, equation_state.unknown->mesh());
    ScalarField result(solution.mesh(), FieldLocation::Face, "equ.diffusionFlux");
    // Orthogonal flux ignores this gradient; the explicit part is already frozen.
    const VectorField unusedGradient(solution.mesh(), FieldLocation::Cell);
    for (const auto& contribution : equation_state.diffusion) {
        const auto& snapshot = *contribution.snapshot;
        auto value = snapshot.iterate;
        value = solution; // keep the assembled boundary constraints
        synchronize(value);
        // Only consumers of a face flux need this second evaluation. The
        // gradient, limiter inputs and boundary conditions are assembly-time
        // snapshots, independent of subsequent solves or field changes.
        const auto correction = diffusionCorrection(
            snapshot.coefficient, snapshot.iterate,
            snapshot.gradient ? &*snapshot.gradient : nullptr, snapshot.method);
        const auto increment =
            frozenDiffusionFlux(snapshot.coefficient, value, correction, unusedGradient);
        result.addScaled(contribution.factor, increment);
    }
    synchronize(result);
    return result;
}

ScalarDiffusion::ScalarDiffusion(ScalarField& unknown, const EquationControl& control,
                                 const std::string& term)
    : equation_(unknown, control), options_(equation_.options(term)) {}

void ScalarDiffusion::assemble(const ScalarField& faceCoefficient,
                               std::optional<double> reference) {
    equation_.reset();
    correction_.reset();
    auto& state = equationState(equation_);
    auto& matrix = state.coefficients;
    if (&faceCoefficient.mesh() != matrix.mesh || faceCoefficient.location() != FieldLocation::Face)
        throw std::invalid_argument("scalar diffusion requires a compatible face coefficient");
    coefficient_.emplace(faceCoefficient);
    synchronize(*coefficient_);
    iterate_.emplace(*state.unknown);
    referenceCell_ = invalid_index;
    addDiffusion(matrix, *coefficient_, *iterate_, *options_.gradient, DiffusionMethod::Orthogonal);
    if (reference) {
        requireFinite(*reference);
        if (hasConstantNullMode(matrix)) {
            referenceCell_ = applyReference(matrix, 0, *reference);
            if (referenceCell_ != invalid_index)
                referenceRhs_ = 0.5 * matrix.diagonal[referenceCell_] * *reference;
        }
    }
    matrix.immutable_matrix = std::make_shared<const impl::MatrixIdentity>();
}

SolveResult ScalarDiffusion::solve(const ScalarField& source, const LinearSolverConfig& config,
                                   const VectorField* reconstructedGradient) {
    auto& state = equationState(equation_);
    if (!state.coefficients.immutable_matrix)
        throw std::logic_error("assemble scalar diffusion before solving");
    correction_.reset();
    *iterate_ = *state.unknown; // preserve the assembled boundary constraints
    synchronize(*iterate_);
    std::optional<VectorField> gradientWorkspace;
    const VectorField* correctionGradient = nullptr;
    if (*options_.diffusion != DiffusionMethod::Orthogonal) {
        if (reconstructedGradient) {
            requireCellField(*reconstructedGradient, iterate_->mesh());
            synchronize(const_cast<VectorField&>(*reconstructedGradient));
            correctionGradient = reconstructedGradient;
        } else {
            gradientWorkspace.emplace(iterate_->mesh(), FieldLocation::Cell, "equ.frozenGradient");
            impl::operators::gradientForExchange(*iterate_, *gradientWorkspace, *options_.gradient);
            backend().synchronizeUpdated(*gradientWorkspace);
            correctionGradient = &*gradientWorkspace;
        }
    }
    refreshDiffusionRhs(state.coefficients, *coefficient_, *iterate_, correctionGradient,
                        *options_.diffusion);
    equ::source(equation_, source);
    if (referenceCell_ != invalid_index)
        state.coefficients.source[referenceCell_] += referenceRhs_;
    auto correction =
        diffusionCorrection(*coefficient_, *iterate_, correctionGradient, *options_.diffusion);
    const auto result = equ::solve(equation_, config);
    correction_.emplace(std::move(correction));
    return result;
}

ScalarField ScalarDiffusion::faceFlux() {
    if (!correction_)
        throw std::logic_error("solve scalar diffusion before recovering its face flux");
    *iterate_ = *equationState(equation_).unknown;
    synchronize(*iterate_);
    ScalarField result(iterate_->mesh(), FieldLocation::Face, "equ.diffusionFlux");
    const VectorField unusedGradient(iterate_->mesh(), FieldLocation::Cell);
    result.addScaled(-1.0,
                     frozenDiffusionFlux(*coefficient_, *iterate_, *correction_, unusedGradient));
    backend().synchronizeUpdated(result);
    return result;
}

template <class T>
SolveResult solveConfigured(const Equation<T>& equation, Field<T>& field,
                            const LinearSolverConfig* config) {
    const auto& equation_state = equationState(equation);
    if (equation_state.unknown != &field)
        throw std::invalid_argument("solve target differs from matrix unknown");
    if constexpr (std::is_same_v<T, double>)
        return backend().solve(equation_state.coefficients, field,
                               config ? *config : equation_state.control.linear);
    else {
        const auto components = backend().solve(equation_state.coefficients, field,
                                                config ? *config : equation_state.control.linear);
        SolveResult result;
        result.status = SolveStatus::Converged;
        for (const auto& c : components) {
            if (!c.healthy())
                result.status = SolveStatus::NumericalFailure;
            else if (!c.converged() && result.status != SolveStatus::NumericalFailure)
                result.status = SolveStatus::MaxIterations;
            result.iterations += c.iterations;
            result.initial_residual = std::hypot(result.initial_residual, c.initial_residual);
            result.final_residual = std::hypot(result.final_residual, c.final_residual);
            result.relative_residual = std::max(result.relative_residual, c.relative_residual);
            result.performance += c.performance;
        }
        return result;
    }
}
template <class T>
SolveResult solve(const Equation<T>& equation, Field<T>& field) {
    return solveConfigured(equation, field, nullptr);
}
template <class T>
SolveResult solve(const Equation<T>& equation, Field<T>& field, const LinearSolverConfig& config) {
    config.validate();
    return solveConfigured(equation, field, &config);
}
template <class T>
SolveResult solve(const Equation<T>& equation) {
    return solveConfigured(equation, *EquationAccess::get(equation).unknown, nullptr);
}
template <class T>
SolveResult solve(const Equation<T>& equation, const LinearSolverConfig& config) {
    config.validate();
    return solveConfigured(equation, *EquationAccess::get(equation).unknown, &config);
}
template <class T>
void Equation<T>::reset() {
    clear(*this);
}
template <class T>
Equation<T> Equation<T>::copy() const {
    return equ::copy(*this);
}
template <class T>
ScalarField Equation<T>::diagonal() const {
    return equ::diagonal(*this);
}
template <class T>
Field<T> Equation<T>::rhs() const {
    return equ::rhs(*this);
}
template <class T>
void Equation<T>::reference(Index cellIndex, double value) {
    if constexpr (std::is_same_v<T, double>)
        equ::reference(*this, cellIndex, value);
    else
        throw std::invalid_argument("reference requires a scalar equation");
}
template <class T>
void Equation<T>::referenceIfUnanchored(double value) {
    if constexpr (std::is_same_v<T, double>)
        equ::reference(*this, value);
    else
        throw std::invalid_argument("reference requires a scalar equation");
}
// Keep template definitions private; the public API supports scalar and vector equations.
#define INSTANTIATE(T)                                                                             \
    template class Equation<T>;                                                                    \
    template void clear(Equation<T>&);                                                             \
    template Equation<T> copy(const Equation<T>&);                                                 \
    template void div(Equation<T>&, const ScalarField&, double, const std::string&,                \
                      const OperatorOptions&);                                                     \
    template void laplacian(Equation<T>&, double, double, const std::string&,                      \
                            const OperatorOptions&);                                               \
    template void laplacian(Equation<T>&, const ScalarField&, double, const std::string&,          \
                            const OperatorOptions&);                                               \
    template void reaction(Equation<T>&, double);                                                  \
    template void reaction(Equation<T>&, const ScalarField&, double);                              \
    template void source(Equation<T>&, T);                                                         \
    template void source(Equation<T>&, const Field<T>&, double);                                   \
    template void ddt(Equation<T>&, double, const Field<T>&, double, TimeMethod, const Field<T>*,  \
                      double);                                                                     \
    template void ddt(Equation<T>&, const ScalarField&, const Field<T>&, double, TimeMethod,       \
                      const Field<T>*, double);                                                    \
    template void ddt(Equation<T>&, double, const time::History<T>&);                              \
    template void ddt(Equation<T>&, const ScalarField&, const time::History<T>&);                  \
    template void scale(Equation<T>&, double);                                                     \
    template void add(Equation<T>&, const Equation<T>&, double);                                   \
    template void addDiagonal(Equation<T>&, const ScalarField&);                                   \
    template void addRhs(Equation<T>&, const Field<T>&);                                           \
    template ScalarField diagonal(const Equation<T>&);                                             \
    template Field<T> rhs(const Equation<T>&);                                                     \
    template Field<T> apply(const Equation<T>&, const Field<T>&);                                  \
    template Field<T> residual(const Equation<T>&, const Field<T>&);                               \
    template void relax(Equation<T>&, const Field<T>&, double);                                    \
    template SolveResult solve(const Equation<T>&, Field<T>&);                                     \
    template SolveResult solve(const Equation<T>&, Field<T>&, const LinearSolverConfig&);          \
    template SolveResult solve(const Equation<T>&);                                                \
    template SolveResult solve(const Equation<T>&, const LinearSolverConfig&);
INSTANTIATE(double)
INSTANTIATE(Vec3)
#undef INSTANTIATE
} // namespace babelsim::equ

namespace babelsim::diagnostics {
namespace {
double square(double value) { return value * value; }
double square(Vec3 value) { return squaredNorm(value); }
}

template <class T>
double relativeResidual(const equ::Equation<T>& equation, const Field<T>& field) {
    const auto ax = equ::applyOwned(equation, field);
    const auto* values = detail::fieldData(ax);
    const auto& source = equ::equationState(equation).coefficients.source;
    double local[3]{}, global[3]{};
    for (Index cell : detail::meshData(field.mesh()).owned_cells) {
        local[0] += square(source[cell] - values[cell]);
        local[1] += square(values[cell]);
        local[2] += square(source[cell]);
    }
    equ::backend().sum(local, global, 3);
    return std::sqrt(global[0]) / std::max(std::sqrt(global[1]) + std::sqrt(global[2]), 1e-30);
}
template double relativeResidual(const equ::Equation<double>&, const ScalarField&);
template double relativeResidual(const equ::Equation<Vec3>&, const VectorField&);
}
