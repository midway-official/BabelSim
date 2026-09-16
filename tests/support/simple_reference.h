// Test-only numerical reference retained for low-level mesh/operator regressions.
// Production steady/transient solvers are independent programs in their main.cpp.
#pragma once

#include "babelsim/config.h"
#include "babelsim/solver.h"
#include "babelsim/geometry.h"
#include "physics/RANS/api.h"

#include <cmath>
#include <stdexcept>

namespace babelsim {

class Case;

// 稳态与瞬态 SIMPLE 共用的私有物理量和算法控制，不属于 BabelSim Solver SDK。
struct FluidProperties {
    double density = 1.0;
    double dynamic_viscosity = 1e-3;

    void validate() const {
        if (!(density > 0.0) || !std::isfinite(density) ||
            !(dynamic_viscosity > 0.0) || !std::isfinite(dynamic_viscosity)) {
            throw std::invalid_argument("density and dynamic viscosity must be positive");
        }
    }
};

struct SimpleControl {
    int max_iterations = 1000;
    int non_orthogonal_corrections = 1;
    double velocity_relaxation = 0.7;
    double pressure_relaxation = 0.3;
    double continuity_tolerance = 1e-8;
    double velocity_tolerance = 1e-7;
    double momentum_tolerance = 1e-6;
    // p' 是本轮用于更新 p 的未松弛压力修正。仅检查速度变化会让不同
    // 分区在压力仍变化时过早停止，因此它必须有独立的外迭代门槛。
    double pressure_correction_tolerance = 1e-6;

    void validate() const {
        if (max_iterations <= 0 || non_orthogonal_corrections < 0 ||
            non_orthogonal_corrections > 20 ||
            !(velocity_relaxation > 0.0 && velocity_relaxation <= 1.0) ||
            !(pressure_relaxation > 0.0 && pressure_relaxation <= 1.0) ||
            !(continuity_tolerance > 0.0) || !std::isfinite(continuity_tolerance) ||
            !(velocity_tolerance > 0.0) || !std::isfinite(velocity_tolerance) ||
            !(pressure_correction_tolerance > 0.0) ||
            !std::isfinite(pressure_correction_tolerance) || !(momentum_tolerance > 0.0) ||
            !std::isfinite(momentum_tolerance)) {
            throw std::invalid_argument("SIMPLE controls are invalid");
        }
    }
};

inline SimpleControl readSimpleControl(const Parameters& settings) {
    SimpleControl result;
    result.max_iterations = settings.integer("maxIterations", result.max_iterations);
    result.non_orthogonal_corrections = settings.integer(
        "nonOrthogonalCorrections", result.non_orthogonal_corrections);
    result.velocity_relaxation = settings.number("velocityRelaxation", result.velocity_relaxation);
    result.pressure_relaxation = settings.number("pressureRelaxation", result.pressure_relaxation);
    result.continuity_tolerance = settings.number("continuityTolerance", result.continuity_tolerance);
    result.velocity_tolerance = settings.number("velocityTolerance", result.velocity_tolerance);
    result.momentum_tolerance = settings.number("momentumTolerance", result.momentum_tolerance);
    result.pressure_correction_tolerance = settings.number(
        "pressureCorrectionTolerance", result.pressure_correction_tolerance);
    result.validate();
    return result;
}

struct IncompressibleFields {
    explicit IncompressibleFields(const Mesh& mesh)
        : velocity(mesh, FieldLocation::Cell, "U"),
          pressure(mesh, FieldLocation::Cell, "p"),
          face_flux(mesh, FieldLocation::Face, "phi")
    {}

    VectorField velocity;
    ScalarField pressure;
    ScalarField face_flux;
};

struct SimpleIterationResult {
    SolveResult velocity;
    SolveResult pressure;
    rans::TransportResult turbulence;
    FluxBalance continuity;
    double relative_velocity_change = 0.0;
    double relative_momentum_residual = 0.0;
    double relative_pressure_correction = 0.0;
    double relative_turbulence_change = 0.0;
    double relative_turbulence_residual = 0.0;
    bool turbulence_active = false;
    bool healthy = false;
    bool linear_converged = false;
    bool converged = false;
};

// 分量 Laplacian 隐式处理 muEff*grad(U)，其余偏应力作为显式通用张量散度。
// 仅由启用涡黏性闭合的动量路径调用；层流保留原 NS 动量离散。
inline void evaluateStressCorrection(
    VectorField& velocity, const ScalarField& phi, const ScalarField& viscosity,
    TensorField& gradient, TensorField& stress, VectorField& divergence)
{
    velocity.setBoundaryFlux(phi);
    gradient.useCalculatedBoundary();
    stress.useCalculatedBoundary();
    gradient = math::grad(velocity);
    stress.evaluate(gradient, [](const Tensor3& g) {
        Tensor3 value = transpose(g);
        const double isotropic = (2.0 / 3.0) * trace(g);
        for (int i = 0; i < 3; ++i) value[i][i] -= isotropic;
        return value;
    });
    stress.assignProduct(viscosity, stress);
    divergence = math::div(stress);
}

}  // babelsim 命名空间

namespace babelsim {
SimpleIterationResult solveIncompressible(IncompressibleFields&, FluidProperties,
    const SimpleControl&, int* iterations = nullptr);
// Ordinary numerical procedure: all nonlinear/correction loops are in solver.cpp.
SimpleIterationResult solveIncompressible(
    VectorField& U, ScalarField& p, ScalarField& phi,
    double density, double viscosity, const SimpleControl& control,
    ScalarField& effectiveViscosity, rans::Model* turbulence = nullptr,
    const VectorField* previous = nullptr, const VectorField* older = nullptr,
    double dt = 0.0, double previousDt = 0.0, TimeMethod timeMethod = TimeMethod::Steady,
    bool log = false, int* iterations = nullptr);
}

#include <iostream>
#include "babelsim/equ.h"
#include <iomanip>
#include <sstream>

namespace babelsim {
SimpleIterationResult solveIncompressible(IncompressibleFields& fields, FluidProperties fluid,
    const SimpleControl& control, int* iterations) {
    ScalarField viscosity(fields.velocity.mesh(),FieldLocation::Cell,"muEffective",fluid.dynamic_viscosity);
    fields.face_flux = math::flux(fields.velocity);
    return solveIncompressible(fields.velocity,fields.pressure,fields.face_flux,
        fluid.density,fluid.dynamic_viscosity,control,viscosity,nullptr,nullptr,nullptr,
        0,0,TimeMethod::Steady,false,iterations);
}
namespace {
void assembleMomentum(equ::Equation<Vec3>& A, VectorField& U,
    const ScalarField& p, const ScalarField& phi, double density, double viscosity,
    const ScalarField& effectiveViscosity, bool turbulent,
    const VectorField* previous, const VectorField* older,
    double dt, double previousDt, TimeMethod timeMethod)
{
    A.reset();
    if (timeMethod != TimeMethod::Steady) {
        if (!previous) throw std::invalid_argument("transient momentum needs explicit history");
        equ::ddt(A, density, *previous, dt, timeMethod, older, previousDt);
    }
    equ::div(A, phi, density);
    if (turbulent) {
        equ::laplacian(A, effectiveViscosity, -1);
        TensorField gradient(U.mesh(), FieldLocation::Cell, "stressGradU");
        TensorField stress(U.mesh(), FieldLocation::Cell, "stressCorrection");
        VectorField divergence(U.mesh(), FieldLocation::Cell, "stressDivergence");
        evaluateStressCorrection(U, phi, effectiveViscosity, gradient, stress, divergence);
        equ::source(A, divergence);
    } else equ::laplacian(A, viscosity, -1);
    VectorField gradP(U.mesh(), FieldLocation::Cell, "gradP");
    gradP = math::grad(p);
    equ::source(A, gradP, -1.0);
}
}

SimpleIterationResult solveIncompressible(
    VectorField& U, ScalarField& p, ScalarField& phi,
    double density, double viscosity, const SimpleControl& control,
    ScalarField& effectiveViscosity, rans::Model* turbulence,
    const VectorField* previous, const VectorField* older,
    double dt, double previousDt, TimeMethod timeMethod, bool log, int* iterations)
{
    FluidProperties{density, viscosity}.validate(); control.validate();
    const auto& mesh=U.mesh();
    if (&p.mesh()!=&mesh || &phi.mesh()!=&mesh ||
        U.location()!=FieldLocation::Cell || p.location()!=FieldLocation::Cell ||
        phi.location()!=FieldLocation::Face)
        throw std::invalid_argument("incompressible field layouts do not match");
    ScalarField pPrime(mesh,FieldLocation::Cell,"pPrime");
    ScalarField rAU(mesh,FieldLocation::Cell,"rAU");
    ScalarField phiHbyA(mesh,FieldLocation::Face,"phiHbyA");
    ScalarField divPhi(mesh,FieldLocation::Cell,"divPhiHbyA");
    VectorField previousIteration(U);
    VectorField gradP(mesh,FieldLocation::Cell,"gradP");
    VectorField rAUgradP(mesh,FieldLocation::Cell,"rAUGradP");
    VectorField rAUgradPf(mesh,FieldLocation::Face,"rAUGradPFace");
    ScalarField rAUf(mesh,FieldLocation::Face,"rAUFace");
    const auto V = geometry::cellVolumes(mesh);
    const bool fixedPressure=setHomogeneousCorrectionBoundaries(pPrime,p);
    auto A=equ::createEquation(U);
    auto P=equ::createEquation(pPrime);
    const auto& methods=numericalMethods();
    const int corrections=methods.diffusionFor(pPrime.name())==DiffusionMethod::Orthogonal
        ? 1 : control.non_orthogonal_corrections+1;
    SimpleIterationResult result;
    result.turbulence_active=turbulence!=nullptr;
    for(int iter=0; iter<control.max_iterations; ++iter) {
        previousIteration=U;
        assembleMomentum(A,U,p,phi,density,viscosity,effectiveViscosity,turbulence,
            previous,older,dt,previousDt,timeMethod);
        equ::relax(A,previousIteration,control.velocity_relaxation);
        // Preserve the existing SIMPLE row normalization explicitly. This scales
        // both sides; V/aP therefore describes precisely the solved matrix.
        equ::scale(A,control.velocity_relaxation);
        rAU=V / A.diagonal();
        result.velocity=equ::solve(A,U);
        if (!diagnostics::all(result.velocity.healthy())) { result.healthy=false; result.converged=false; return result; }

        // Rhie-Chow momentum interpolation, preserving physical boundary fluxes.
        gradP = math::grad(p);
        phiHbyA = math::flux(U);
        rAUgradP=rAU*gradP;
        rAUgradPf = math::interpolate(rAUgradP);
        rAUf = math::interpolate(rAU);
        math::add(math::flux(rAUgradPf),phiHbyA,math::FaceRegion::Interior);
        math::subtract(math::flux(rAUf, p, gradP),phiHbyA,math::FaceRegion::Interior);
        divPhi = math::div(phiHbyA);
        pPrime.fill(0);
        bool pressureHealthy=true, pressureConverged=true;
        for(int correction=0; correction<corrections; ++correction) {
            P.reset();
            equ::laplacian(P,rAU, -1);
            equ::source(P,divPhi,-1.0);
            if(!fixedPressure) P.reference(0, 0.0);
            result.pressure=equ::solve(P,pPrime);
            pressureHealthy=pressureHealthy && result.pressure.healthy();
            pressureConverged=pressureConverged && result.pressure.converged();
        }
        if (!diagnostics::all(pressureHealthy)) { result.healthy=false; result.converged=false; return result; }
        p+=control.pressure_relaxation*pPrime;
        // Preserve the current transient pressure/velocity relaxation coupling.
        auto pressureFlux=equ::faceFlux(P,pPrime);
        if(timeMethod!=TimeMethod::Steady) {
            pPrime*=control.pressure_relaxation;
            pressureFlux*=control.pressure_relaxation;
        }
        math::subtract(rAU,math::grad(pPrime),U);
        phi=phiHbyA+pressureFlux;

        if(turbulence) {
            result.turbulence=turbulence->solveTransport();
            result.relative_turbulence_change=result.turbulence.relativeChange();
            result.relative_turbulence_residual=result.turbulence.initialResidual();
        }
        result.relative_velocity_change=diagnostics::relativeChange(U,previousIteration);
        result.relative_pressure_correction=diagnostics::relativeMagnitude(pPrime,p);
        result.continuity=diagnostics::fluxBalance(phi);
        // Nonlinear residual uses freshly assembled physical coefficients,
        // separate from the relaxed linear system solved above.
        assembleMomentum(A,U,p,phi,density,viscosity,effectiveViscosity,turbulence,
            previous,older,dt,previousDt,timeMethod);
        result.relative_momentum_residual=diagnostics::relativeResidual(A,U);
        const bool turbulenceHealthy=!turbulence || (result.turbulence.healthy() &&
            std::isfinite(result.relative_turbulence_change) && std::isfinite(result.relative_turbulence_residual));
        result.healthy=diagnostics::all(pressureHealthy && turbulenceHealthy && result.velocity.healthy() &&
            std::isfinite(result.relative_velocity_change) && std::isfinite(result.relative_momentum_residual) &&
            std::isfinite(result.relative_pressure_correction) && std::isfinite(result.continuity.relative));
        result.linear_converged=diagnostics::all(pressureConverged && result.velocity.converged() &&
            (!turbulence || result.turbulence.linearConverged()));
        result.converged=result.healthy && result.linear_converged && diagnostics::all(
            result.continuity.relative<=control.continuity_tolerance &&
            result.relative_velocity_change<=control.velocity_tolerance &&
            result.relative_momentum_residual<=control.momentum_tolerance &&
            result.relative_pressure_correction<=control.pressure_correction_tolerance &&
            (!turbulence || (result.relative_turbulence_change<=turbulence->tolerance() &&
                result.relative_turbulence_residual<=turbulence->tolerance())));
        if(iterations) *iterations=iter+1;
        if(log && (iter==0 || (iter+1)%100==0 || result.converged || !result.healthy || iter+1==control.max_iterations)) {
            std::ostringstream message;
            message << (timeMethod==TimeMethod::Steady ? "SIMPLE " : "Transient SIMPLE ") << iter+1
                << std::scientific << std::setprecision(6)
                << " mass=" << result.continuity.relative << " dU=" << result.relative_velocity_change
                << " rU=" << result.relative_momentum_residual << " dP=" << result.relative_pressure_correction
                << " linP=" << result.pressure.relative_residual << " dTurb=" << result.relative_turbulence_change
                << " rTurb=" << result.relative_turbulence_residual
                << " linear=" << (result.linear_converged ? "ok" : "inexact")
                << " converged=" << (result.converged ? "true" : "false");
            if (primaryProcess()) std::cout << message.str() << '\n';
        }
        if(result.converged || !result.healthy) break;
    }
    return result;
}
}
