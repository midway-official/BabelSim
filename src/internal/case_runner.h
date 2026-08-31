#pragma once

#include "babelsim/case.h"
#include "babelsim/parallel.h"

#include <string>

namespace babelsim {

// 各 Physics 目录拥有自己的“Case -> Field -> RunTime -> Output”适配器；通用启动器
// 只负责参数、MPI 生命周期和显式分派，不需要了解具体 Field 或数值控制。
int runHeatCase(const CaseDefinition& definition, const ParallelContext& parallel,
                const std::string& requested_time);
int runSimpleCase(const CaseDefinition& definition, const ParallelContext& parallel,
                  const std::string& requested_time);
int runTransportCase(const CaseDefinition& definition, const ParallelContext& parallel,
                     const std::string& requested_time);

}  // babelsim 命名空间
