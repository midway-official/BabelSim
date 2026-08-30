#pragma once

#include "babelsim/case.h"
#include "babelsim/thermal.h"

namespace babelsim {

struct ThermalCaseControl {
    ThermalProperties material;
    double volumetric_source = 0.0;
    RuntimeControl runtime;
};

ThermalCaseControl readThermalCase(const CaseDefinition& definition);

}  // babelsim 命名空间
