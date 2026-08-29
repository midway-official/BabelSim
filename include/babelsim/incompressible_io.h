#pragma once

#include "babelsim/case.h"
#include "babelsim/incompressible.h"

namespace babelsim {

struct IncompressibleCaseControl {
    FluidProperties fluid;
    Methods methods;
    SimpleControl simple;
};

IncompressibleCaseControl readIncompressibleCase(const CaseDefinition& definition);

}  // babelsim 命名空间
