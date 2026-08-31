#pragma once

#include "babelsim/case.h"
#include "babelsim/runtime.h"

namespace babelsim {

struct TransportCaseControl {
    double storage = 1.0;
    double diffusivity = 0.0;
    double source = 0.0;
    RuntimeControl runtime;
};

TransportCaseControl readTransportCase(const CaseDefinition& definition);

}  // babelsim 命名空间
