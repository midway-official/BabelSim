#define _GNU_SOURCE
#include <petscksp.h>
#include <dlfcn.h>
#include <stdio.h>
#include <string.h>

/* Diagnostic overlay only. Production source is unchanged during screening.
 * Enable the public PETSc CG recurrence for every distributed CG system.
 * Do not alter the norm, tolerances, preconditioner, or convergence callback.
 */
PetscErrorCode KSPSetType(KSP ksp, KSPType type) {
    typedef PetscErrorCode (*SetType)(KSP, KSPType);
    static SetType original = NULL;
    if (!original) original = (SetType)dlsym(RTLD_NEXT, "KSPSetType");
    if (!original) return PETSC_ERR_LIB;
    PetscErrorCode error = original(ksp, type);
    if (error || !type || strcmp(type, KSPCG)) return error;
    MPI_Comm comm = PetscObjectComm((PetscObject)ksp);
    int size = 1, rank = 0;
    if (MPI_Comm_size(comm, &size) != MPI_SUCCESS) return PETSC_ERR_MPI;
    if (size == 1) return PETSC_SUCCESS;
    error = KSPCGUseSingleReduction(ksp, PETSC_TRUE);
    if (!error && MPI_Comm_rank(comm, &rank) == MPI_SUCCESS && rank == 0)
        fprintf(stderr, "DIAGNOSTIC distributed CG single reduction enabled\n");
    return error;
}
