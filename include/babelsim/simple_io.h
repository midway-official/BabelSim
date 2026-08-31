#pragma once

#include "babelsim/case.h"
#include "babelsim/simple.h"

namespace babelsim {

struct SimpleCaseControl {
    FluidProperties fluid;
    RuntimeControl runtime;
    SimpleControl simple;
};

SimpleCaseControl readSimpleCase(const CaseDefinition& definition);

}  // babelsim 命名空间
