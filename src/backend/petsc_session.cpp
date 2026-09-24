#include "internal/petsc_session.h"

#include <petscsys.h>

#include <mutex>
#include <stdexcept>
#include <string>

namespace babelsim::detail {
namespace {

std::mutex session_mutex;
bool initialized_by_babelsim = false;
bool return_handler_pushed = false;

void checkPetsc(PetscErrorCode code, const char* operation) {
    if (code == 0) return;
    const char* message = nullptr;
    PetscErrorMessage(code, &message, nullptr);
    throw std::runtime_error(std::string(operation) + " failed with PETSc error " +
                             std::to_string(code) +
                             (message ? std::string(": ") + message : std::string{}));
}

}  // namespace

void ensurePetscSession() {
    std::lock_guard<std::mutex> lock(session_mutex);
    PetscBool initialized = PETSC_FALSE;
    checkPetsc(PetscInitialized(&initialized), "PetscInitialized");
    if (!initialized) {
        // Do not let PETSc consume BabelSim's command line; all solver settings
        // are parsed and validated from the case configuration.
        checkPetsc(PetscInitializeNoArguments(), "PetscInitializeNoArguments");
        initialized_by_babelsim = true;
    }
    if (!return_handler_pushed) {
        checkPetsc(PetscPushErrorHandler(PetscReturnErrorHandler, nullptr),
                   "PetscPushErrorHandler");
        return_handler_pushed = true;
    }
}

void finalizePetscSession() {
    std::lock_guard<std::mutex> lock(session_mutex);
    PetscBool initialized = PETSC_FALSE;
    PetscBool finalized = PETSC_FALSE;
    checkPetsc(PetscInitialized(&initialized), "PetscInitialized");
    checkPetsc(PetscFinalized(&finalized), "PetscFinalized");
    if (initialized && !finalized && initialized_by_babelsim) {
        checkPetsc(PetscFinalize(), "PetscFinalize");
        initialized_by_babelsim = false;
        return_handler_pushed = false;
    }
}

}  // namespace babelsim::detail
