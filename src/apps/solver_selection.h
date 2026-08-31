#pragma once

// 应用层的显式分派声明，不属于普通 Solver 或数值后端的接口。
namespace babelsim {
class Case;

// 显式分派表的声明。新增求解器只需提供同样的函数签名，无基类、注册器或工厂。
int runSolver(Case& problem);
}  // babelsim 命名空间
