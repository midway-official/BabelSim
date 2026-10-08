#include <petscksp.h>
#include <dlfcn.h>
#include <cstdlib>

// Diagnostic-only LD_PRELOAD observer. It calls the unchanged real KSPSetUp
// first and records a public PCView for the first collective Hypre setup.
// Formal timing runs do
// not load this library. No PETSc/Hypre private data is accessed or modified.
extern "C" PetscErrorCode KSPSetUp(KSP ksp) {
    using Setup = PetscErrorCode (*)(KSP);
    static const auto actual = reinterpret_cast<Setup>(dlsym(RTLD_NEXT, "KSPSetUp"));
    if (!actual) return PETSC_ERR_LIB;
    const auto error = actual(ksp);
    if (error) return error;
    const char* path = std::getenv("BABELSIM_PROFILE_PC_VIEW");
    if (!path) return 0;
    PC pc = nullptr;
    PetscBool hypre = PETSC_FALSE;
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PetscObjectTypeCompare(reinterpret_cast<PetscObject>(pc), PCHYPRE, &hypre));
    if (!hypre) return 0;
    // Local pointer reuse can differ between ranks. It must never decide
    // participation in collective PCView/viewer operations.
    static bool viewed = false;
    if (viewed) return 0;
    viewed = true;
    PetscViewer viewer = nullptr;
    PetscCall(PetscViewerASCIIOpen(PetscObjectComm(reinterpret_cast<PetscObject>(pc)), path, &viewer));
    const auto view_error = PCView(pc, viewer);
    PetscViewerDestroy(&viewer);
    return view_error;
}
