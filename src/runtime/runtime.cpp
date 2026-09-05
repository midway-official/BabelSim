#include "babelsim/runtime.h"
#include "babelsim/mpi_support.h"
#include "babelsim/parallel.h"
#include "internal/compute_backend.h"
#include "internal/fvm_execution.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace babelsim {
namespace {
thread_local RunTime* active_run_time = nullptr;
}  // 匿名命名空间

struct RunTime::Implementation {
    Implementation(const Mesh& mesh_value, RuntimeControl settings, ParallelContext parallel_value)
        : mesh(&mesh_value), control(std::move(settings)), primary_rank(parallel_value.rank == 0),
          fvm(mesh_value, control.methods,
              detail::makeComputeBackend(
                  mesh_value, control.scalar_solver, control.vector_solver,
                  std::move(parallel_value)),
              control.time.delta_t) {}
    const Mesh* mesh;
    RuntimeControl control;
    bool primary_rank;
    detail::FvmExecution fvm;
    double current_time = 0.0;
    double current_delta_t = 0.0;
    int current_step = 0;
};

void RuntimeControl::validate() const {
    time.validate();
    const double steps = (time.end_time - time.start_time) / time.delta_t;
    if (!std::isfinite(steps) || steps > std::numeric_limits<int>::max())
        throw std::invalid_argument("time interval contains too many steps");
    // 当前 BDF2 是等步长离散，不能把缩短的末步冒充等步长 BDF2。
    if (methods.time == TimeMethod::BDF2 &&
        std::abs(steps - std::round(steps)) >
            64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, steps))
        throw std::invalid_argument("BDF2 requires an integral number of uniform time steps");
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
    control.validate();
    m_implementation = std::make_unique<Implementation>(mesh, std::move(control), parallel);
    m_implementation->current_time = m_implementation->control.time.start_time;
    m_implementation->current_delta_t = m_implementation->control.time.delta_t;
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
        throw std::logic_error("solve/math/diagnostics require an active RunTime");
    }
    return *active_run_time;
}

const Mesh& RunTime::mesh() const { return *m_implementation->mesh; }
const Methods& RunTime::methods() const { return m_implementation->control.methods; }
double RunTime::time() const { return m_implementation->current_time; }
double RunTime::deltaT() const { return m_implementation->current_delta_t; }
int RunTime::step() const { return m_implementation->current_step; }
bool RunTime::primary() const { return m_implementation->primary_rank; }


bool RunTime::loop() {
    Implementation& state = *m_implementation;
    const TimeControl& control = state.control.time;
    const double tolerance = 32.0 * std::numeric_limits<double>::epsilon() *
        std::max({std::abs(control.start_time), std::abs(control.end_time), control.delta_t});
    const double remaining = control.end_time - state.current_time;
    if (remaining <= tolerance) return false;

    const double next = static_cast<double>(static_cast<long double>(control.start_time) +
        static_cast<long double>(state.current_step + 1) * control.delta_t);
    state.current_delta_t = remaining < control.delta_t - tolerance ? remaining : control.delta_t;
    state.current_time = next >= control.end_time - tolerance ? control.end_time : next;
    state.fvm.beginStep(state.current_delta_t);
    ++state.current_step;
    return true;
}


detail::FvmExecution& detail::execution() {
    return RunTime::current().m_implementation->fvm;
}
}  // babelsim 命名空间
