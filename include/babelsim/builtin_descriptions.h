#pragma once
#include "babelsim/model_description.h"

namespace babelsim::builtin {
const ModuleDescription& heat();
const ModuleDescription& transport();
const ModuleDescription& simple();
const ModuleDescription& transientSimple();
const ModuleDescription& piso();
const ModuleDescription& spalartAllmaras();
const ModuleDescription& kOmega();
const ModuleDescription& kEpsilon();
} // namespace babelsim::builtin
