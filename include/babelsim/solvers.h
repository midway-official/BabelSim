#pragma once

namespace babelsim {
class Case;

// 显式分派表的声明。新增求解器只需提供同样的函数签名，无基类、注册器或工厂。
int runSolver(Case& problem);
}  // babelsim 命名空间
