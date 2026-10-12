#include "babelsim/case.h"
#include "babelsim/equ.h"
#include "babelsim/history.h"
#include "babelsim/runtime.h"
#include "babelsim/time_discretization.h"
#include "babelsim/time_step.h"
#include "internal/boundary_evaluation.h"
#include "internal/field_access.h"
#include "test_util.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>

using namespace babelsim;

namespace {
template<class Function>
void rejects(Function function, const std::string& message) {
    bool rejected = false;
    try { function(); } catch (const std::exception&) { rejected = true; }
    require(rejected, message);
}

void requireValues(const ScalarField& field, double value, const std::string& message) {
    for (std::size_t index = 0; index < field.size(); ++index)
        require(near(detail::fieldData(field)[index], value), message);
}

void requireValues(const VectorField& field, Vec3 value, const std::string& message) {
    for (std::size_t index = 0; index < field.size(); ++index)
        require(near(detail::fieldData(field)[index], value), message);
}

void checkWeightsAndPrediction() {
    const auto first = time::weights(TimeMethod::BDF2, 0.2);
    require(first.method == TimeMethod::Euler && near(first.current, 5.0)
        && near(first.previous, 5.0) && near(first.older, 0.0), "BDF2 startup is not Euler");
    const auto uniform = time::weights(TimeMethod::BDF2, 0.2, 0.2, 2);
    require(uniform.method == TimeMethod::BDF2 && near(uniform.current, 7.5)
        && near(uniform.previous, 10.0) && near(uniform.older, -2.5), "uniform BDF2 weights");
    const auto variable = time::weights(TimeMethod::BDF2, 0.5, 1.0, 2);
    // Independent quadratic derivative: q(0)=0, q(1)=1, q(1.5)=2.25.
    require(near(variable.current * 2.25 - variable.previous * 1.0 - variable.older * 0.0, 3.0),
        "variable-step BDF2 does not differentiate a quadratic exactly");
    require(near(variable.current - variable.previous - variable.older, 0.0),
        "BDF2 derivative of a constant is nonzero");
    for (double bad : {0.0, -0.1, std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::quiet_NaN()})
        rejects([&] { (void)time::weights(TimeMethod::Euler, bad); }, "invalid dt accepted");
    rejects([&] { (void)time::weights(TimeMethod::Euler, std::numeric_limits<double>::denorm_min()); },
        "nonfinite Euler weights accepted");
    rejects([&] { (void)time::weights(TimeMethod::BDF2, 0.1, -0.1, 2); },
        "invalid previous dt accepted");

    const auto mesh = makeHexBox({2, 1, 1}, {0, 0, 0}, {4, 1, 1});
    RuntimeControl settings;
    settings.methods.time = TimeMethod::BDF2;
    auto runtime = RunTime::forMesh(mesh, settings);
    ScalarField x(mesh, FieldLocation::Cell, "x", 1.0);
    x.boundary("minus_x") = fixedValue(11.0);
    auto history = time::history(x);
    rejects([&] { (void)time::weights(TimeMethod::BDF2, history); }, "unsaved history accepted");
    history.save(x, 0.2);
    const auto compareAssembly = [&] {
        auto equation = testEquation(x);
        equ::ddt(equation, 3.0, history);
        const auto w = time::weights(TimeMethod::BDF2, history);
        const auto diagonal = equation.diagonal();
        const auto rhs = equation.rhs();
        for (Index cell = 0; cell < mesh.cellCount(); ++cell) {
            const double capacityVolume = 3.0 * mesh.cellVolume(cell);
            require(near(detail::fieldData(diagonal)[cell], capacityVolume * w.current),
                "public time weights disagree with production time diagonal");
            const double expected = capacityVolume * (w.previous * detail::fieldData(history.previous())[cell]
                + w.older * detail::fieldData(history.older())[cell]);
            require(near(detail::fieldData(rhs)[cell], expected),
                "public time weights disagree with production time RHS");
        }
    };
    compareAssembly();
    x.fill(1.4); // q(t)=1+2t at t=0.2.
    history.save(x, 0.125);
    compareAssembly();
    time::predict(x, history, TimeMethod::BDF2);
    requireValues(x, 1.65, "variable-step linear prediction");
    require(x.boundary("minus_x").type == BoundaryType::FixedValue
        && near(x.boundary("minus_x").value, 11.0), "prediction changed physical boundary constraints");
    time::predict(x, history, TimeMethod::Euler);
    requireValues(x, 1.4, "Euler prediction must use the previous value");

    VectorField velocity(mesh, FieldLocation::Cell, "velocity", {1, 2, 3});
    auto vectorHistory = time::history(velocity);
    vectorHistory.save(velocity, 0.2);
    velocity.fill({3, 6, 9});
    vectorHistory.save(velocity, 0.1);
    time::predict(velocity, vectorHistory, TimeMethod::BDF2);
    requireValues(velocity, {4, 8, 12}, "vector time prediction");
}

void checkHistoryGroups() {
    const auto mesh = makeHexBox({2, 1, 1}, {0, 0, 0}, {2, 1, 1});
    ScalarField scalar(mesh, FieldLocation::Cell, "scalar", 1.0);
    VectorField velocity(mesh, FieldLocation::Cell, "velocity", {1, 2, 3});
    ScalarField flux(mesh, FieldLocation::Face, "flux", 2.0);
    auto scalarHistory = time::history(scalar);
    auto velocityHistory = time::history(velocity);
    auto fluxHistory = time::history(flux);
    time::HistoryGroup group;
    group.add(scalar, scalarHistory);
    group.add(velocity, velocityHistory);
    group.add(flux, fluxHistory);
    rejects([&] { group.add(scalar, scalarHistory); }, "duplicate history registration accepted");
    const auto otherMesh = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
    ScalarField other(otherMesh, FieldLocation::Cell, "other");
    auto otherHistory = time::history(other);
    rejects([&] { group.add(other, otherHistory); }, "group accepted fields from another mesh");
    rejects([&] { group.track(scalar); }, "history field registered again for tracking");
    group.save(0.1);
    group.validate();
    require(scalarHistory.levels() == 1 && velocityHistory.levels() == 1 && fluxHistory.levels() == 1,
        "group did not advance all histories");
    rejects([&] { group.save(-0.1); }, "group accepted invalid dt");
    require(scalarHistory.levels() == 1 && velocityHistory.levels() == 1 && fluxHistory.levels() == 1,
        "invalid group save partially advanced history");
    scalarHistory.save(scalar, 0.1);
    rejects([&] { group.validate(); }, "group failed to detect independently advanced history");
    rejects([&] { group.save(0.1); }, "group silently repaired mismatched history levels");
    require(velocityHistory.levels() == 1 && fluxHistory.levels() == 1,
        "mismatched save partially advanced other histories");

    auto oneHistory = time::history(other);
    time::HistoryGroup one;
    one.add(other, oneHistory);
    one.save(0.1);
    oneHistory.save(other, 0.2);
    rejects([&] { one.validate(); }, "one-member group failed to detect external history drift");

    auto a = time::history(scalar);
    auto b = time::history(flux);
    a.save(scalar, 0.1);
    b.save(flux, 0.2);
    time::HistoryGroup inconsistent;
    inconsistent.add(scalar, a);
    rejects([&] { inconsistent.add(flux, b); inconsistent.validate(); },
        "same-level group accepted mismatched dt");
}

struct ModelState {
    int updates = 0;
    double coefficient = 1.0;
    std::string phase = "initial";
};

struct CaseFixture {
    std::filesystem::path root = std::filesystem::temp_directory_path()
        / ("babelsim-time-tools-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    CaseFixture() {
        std::filesystem::copy("cases/heat", root, std::filesystem::copy_options::recursive);
        std::ofstream(root / "control.bs") << "startTime 0\nendTime 0.25\ndeltaT 0.1\n";
        std::ofstream(root / "numerics/methods.bs") << "time bdf2\n"
            "equation.temperature.interpolation linear\n"
            "equation.temperature.gradient greenGauss\n"
            "equation.temperature.convection upwind\n"
            "equation.temperature.diffusion orthogonal\n";
    }
    ~CaseFixture() { std::filesystem::remove_all(root); }
};

void checkStepTransactions() {
    using Step = decltype(time::beginStep(std::declval<TimeStepper&>(), std::declval<time::HistoryGroup&>()));
    static_assert(!std::is_copy_constructible_v<Step> && std::is_move_constructible_v<Step>,
        "step trials must have unique movable ownership");
    CaseFixture fixture;
    Case problem(fixture.root);
    auto& temperature = problem.scalarField("T");
    auto& velocity = problem.createVectorField("velocity", {1, 2, 3});
    auto& flux = problem.createFaceScalarField("flux");
    auto& derived = problem.createScalarField("derived", 4.0);
    temperature.fill(1.0);
    temperature.boundary(0) = fixedValue(8.0);
    flux.fill(2.0);
    velocity.boundary(0) = BoundaryCondition<Vec3>::inletOutlet({5, 6, 7});
    velocity.setBoundaryFlux(flux);
    derived.useCalculatedBoundary();
    derived.fill(4.0);
    const Index boundaryFace = *problem.mesh().patchFaces(0).begin();
    (void)problem.physics().positive("density");
    (void)problem.physics().positive("heatCapacity");
    (void)problem.physics().nonnegative("conductivity");
    (void)problem.physics().number("source");
    auto equation = equ::createEquation(problem, "temperature", temperature, {"diffusion"});
    auto clock = time::start(problem);
    auto temperatureHistory = time::history(temperature);
    auto velocityHistory = time::history(velocity);
    auto fluxHistory = time::history(flux);
    ModelState model;
    time::HistoryGroup group;
    group.add(temperature, temperatureHistory);
    group.add(velocity, velocityHistory);
    group.add(flux, fluxHistory);
    group.track(derived);
    group.watch(model);
    rejects([&] { group.track(derived); }, "duplicate tracked field accepted");
    rejects([&] { group.watch(model); }, "duplicate watched model state accepted");
    rejects([&] { group.save(clock); }, "history saved before advancing clock");
    rejects([&] { auto invalid = time::beginStep(clock, group, -0.1); }, "invalid trial dt accepted");
    {
        const auto otherMesh = makeHexBox({1, 1, 1}, {0, 0, 0}, {1, 1, 1});
        ScalarField other(otherMesh, FieldLocation::Cell, "other");
        auto otherHistory = time::history(other);
        time::HistoryGroup wrongGroup;
        wrongGroup.add(other, otherHistory);
        rejects([&] { auto invalid = time::beginStep(clock, wrongGroup); },
            "trial accepted histories from a different mesh");
        require(clock.step() == 0 && near(clock.value(), 0.0) && otherHistory.levels() == 0,
            "failed trial construction partially advanced clock/history");
    }
    {
        auto trial = time::beginStep(clock, group);
        require(near(clock.value(), 0.1) && clock.step() == 1 && temperatureHistory.levels() == 1,
            "trial did not advance clock and grouped histories together");
        rejects([&] { auto nested = time::beginStep(clock, group); }, "nested trial accepted");
        rejects([&] { clock.advance(); }, "clock advanced during pending trial");
        rejects([&] { group.save(clock); }, "group saved again during pending trial");
        temperature.fill(9.0);
        temperature.boundary(0) = fixedValue(99.0);
        temperature.useCalculatedBoundary();
        velocity.fill({9, 9, 9});
        flux.fill(-12.0);
        velocity.setBoundaryFlux(flux);
        derived.fill(19.0);
        model = {7, 9.0, "failed"};
        auto moved = std::move(trial);
        moved.reject();
        rejects([&] { trial.accept(); }, "rejected trial accepted twice");
        rejects([&] { moved.reject(); }, "rejected trial rejected twice");
    }
    require(near(clock.value(), 0.0) && clock.step() == 0 && near(clock.dt(), 0.1)
        && near(problem.time(), 0.0) && problem.step() == 0, "rejection did not restore clock/Case");
    requireValues(temperature, 1.0, "rejection did not restore unknown");
    requireValues(velocity, {1, 2, 3}, "rejection did not restore vector unknown");
    requireValues(flux, 2.0, "rejection did not restore face field");
    requireValues(derived, 4.0, "rejection did not restore tracked derived field");
    require(!temperature.calculatedBoundary(), "rejection did not restore physical boundary mode");
    require(near(temperature.boundary(0).value, 8.0), "rejection did not restore boundary constraints");
    require(derived.calculatedBoundary() && near(boundaryFaceValue(derived, boundaryFace), 4.0),
        "rejection did not restore calculated boundary traces");
    require(near(boundaryFaceValue(velocity, boundaryFace), Vec3{1, 2, 3}),
        "rejection did not restore inletOutlet flux context");
    require(temperatureHistory.levels() == 0 && velocityHistory.levels() == 0 && fluxHistory.levels() == 0
        && near(temperatureHistory.dt(), 0.0) && near(temperatureHistory.previousDt(), 0.0),
        "rejection did not restore history metadata");
    require(model.updates == 0 && near(model.coefficient, 1.0) && model.phase == "initial",
        "rejection did not restore model state");
    group.validate();

    {
        auto retry = time::beginStep(clock, group, 0.05);
        equation.reset();
        equ::ddt(equation, 1.0, temperatureHistory);
        equ::source(equation, 1.0);
        require(equ::solve(equation).converged(), "retry solve failed");
        requireValues(temperature, 1.05, "retry did not use restored Euler startup history");
        model.updates = 1;
        retry.accept();
        rejects([&] { retry.accept(); }, "accepted trial accepted twice");
        rejects([&] { retry.reject(); }, "accepted trial rejected afterwards");
    }
    require(clock.step() == 1 && near(clock.value(), 0.05) && temperatureHistory.levels() == 1,
        "accepted retry lost time/history");
    {
        auto trial = time::beginStep(clock, group, 0.03);
        temperatureHistory.save(temperature, 0.03);
        rejects([&] { trial.accept(); }, "accept ignored independent history advancement");
        trial.reject();
        require(clock.step() == 1 && temperatureHistory.levels() == 1,
            "reject after failed acceptance did not restore accepted state");
    }
    {
        auto trial = time::beginStep(clock, group); // Default dt follows accepted variable step.
        require(near(clock.value(), 0.15) && near(clock.dt(), 0.1), "default dt jumped over accepted retry");
        temperature.fill(80.0);
        model.updates = 99;
        // Unaccepted scope exit rolls back the complete trial.
    }
    require(near(clock.value(), 0.05) && near(clock.dt(), 0.05) && clock.step() == 1,
        "trial destructor did not restore clock");
    requireValues(temperature, 1.05, "trial destructor did not restore field");
    require(temperatureHistory.levels() == 1 && model.updates == 1, "trial destructor lost accepted state");
    {
        auto trial = time::beginStep(clock, group, 0.1);
        const auto coefficients = time::weights(TimeMethod::BDF2, temperatureHistory);
        require(near(coefficients.current, 50.0 / 3.0) && near(coefficients.previous, 30.0)
            && near(coefficients.older, -40.0 / 3.0), "retry corrupted variable-step BDF2 weights");
        equation.reset();
        equ::ddt(equation, 1.0, temperatureHistory);
        equ::source(equation, 1.0);
        require(equ::solve(equation).converged(), "second accepted solve failed");
        requireValues(temperature, 1.15, "BDF2 retry differs from direct accepted trajectory");
        trial.accept();
    }
    {
        auto last = time::beginStep(clock, group, 1.0);
        require(near(clock.value(), 0.25) && near(clock.dt(), 0.1), "requested dt did not respect endTime");
        last.accept();
    }
    require(clock.finished() && clock.step() == 3 && temperatureHistory.levels() == 3,
        "accepted shortened interval did not finish correctly");
    rejects([&] { auto extra = time::beginStep(clock, group); }, "trial beyond endTime accepted");
}

void checkClockSaveGuards() {
    CaseFixture fixture;
    Case problem(fixture.root);
    auto& field = problem.scalarField("T");
    (void)problem.physics().positive("density");
    (void)problem.physics().positive("heatCapacity");
    (void)problem.physics().nonnegative("conductivity");
    (void)problem.physics().number("source");
    auto equation = equ::createEquation(problem, "temperature", field, {"diffusion"});
    auto clock = time::start(problem);
    auto otherClock = time::start(problem);
    auto history = time::history(field);
    time::HistoryGroup group;
    group.add(field, history);
    clock.advance();
    group.save(clock);
    rejects([&] { group.save(clock); }, "one physical step saved twice");
    rejects([&] { group.save(clock.dt()); }, "clock-bound group accepted an untracked dt save");
    otherClock.advance();
    rejects([&] { group.save(otherClock); }, "group rebound to a different clock");
    clock.advance();
    clock.advance();
    rejects([&] { group.save(clock); }, "history group skipped a physical step");
    require(history.levels() == 1, "clock-save guard mutated history");
}
} // namespace

int main() {
    checkWeightsAndPrediction();
    checkHistoryGroups();
    checkStepTransactions();
    checkClockSaveGuards();
    std::cout << "time_tools_test: public weights/production assembly, prediction, grouped histories, "
        "rejected-step restoration and variable-step retry passed\n";
}
