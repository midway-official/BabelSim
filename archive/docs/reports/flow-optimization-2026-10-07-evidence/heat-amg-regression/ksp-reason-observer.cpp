#include <petscksp.h>
#include <dlfcn.h>

// Diagnostic-only public API observer. Never loaded by formal measurements.
extern "C" PetscErrorCode KSPSolve(KSP ksp, Vec rhs, Vec solution) {
    using Solve = PetscErrorCode (*)(KSP, Vec, Vec);
    static const auto actual = reinterpret_cast<Solve>(dlsym(RTLD_NEXT, "KSPSolve"));
    if (!actual) return PETSC_ERR_LIB;
    const auto error = actual(ksp, rhs, solution);
    if (error) return error;
    KSPConvergedReason reason;
    PetscReal recurrence = 0, truth = 0;
    PetscInt iterations = 0;
    Mat matrix = nullptr;
    Vec residual = nullptr;
    PetscCall(KSPGetConvergedReason(ksp, &reason));
    PetscCall(KSPGetResidualNorm(ksp, &recurrence));
    PetscCall(KSPGetIterationNumber(ksp, &iterations));
    PetscCall(KSPGetOperators(ksp, &matrix, nullptr));
    PetscCall(VecDuplicate(rhs, &residual));
    PetscCall(MatMult(matrix, solution, residual));
    PetscCall(VecAXPY(residual, -1.0, rhs));
    PetscCall(VecNorm(residual, NORM_2, &truth));
    PetscCall(PetscPrintf(PetscObjectComm(reinterpret_cast<PetscObject>(ksp)),
        "DIAGNOSTIC KSP reason=%d iterations=%" PetscInt_FMT " recurrence=%g true=%g\n",
        static_cast<int>(reason), iterations, static_cast<double>(recurrence),
        static_cast<double>(truth)));
    PetscCall(VecDestroy(&residual));
    return 0;
}
