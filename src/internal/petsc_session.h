#pragma once

namespace babelsim::detail {

// PETSc is process-global. The numerical backend initializes it lazily, while
// applications call finalizePetscSession() after destroying the Case/RunTime
// and before finalizing MPI.
void ensurePetscSession();
void finalizePetscSession();

}  // namespace babelsim::detail
