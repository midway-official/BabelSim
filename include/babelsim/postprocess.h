#pragma once

namespace babelsim {
// 独立后处理程序入口；文件选择、原生结果读取和格式转换都属于 IO 模块。
int runPostprocess(int argc, char* argv[]);
}  // babelsim 命名空间
