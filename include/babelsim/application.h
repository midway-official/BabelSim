#pragma once

namespace babelsim {
class Case;

// 每个 Solver 在自己的源文件命名空间作用域注册一次，名称使用字符串字面量。
// 注册只连接描述项：无堆分配、异常或 MPI 调用，校验与分派在 runApplication 内完成。
// 对象及名称必须在应用运行期间保持有效；不支持运行中的并发注册或动态插件卸载。
class SolverRegistration {
public:
    SolverRegistration(const char* name, int (*run)(Case&)) noexcept;
    ~SolverRegistration() noexcept;
    SolverRegistration(const SolverRegistration&) = delete;
    SolverRegistration& operator=(const SolverRegistration&) = delete;

private:
    const char* m_name;
    int (*m_run)(Case&);
    mutable const SolverRegistration* m_next;
    static const SolverRegistration*& first() noexcept;
    friend int runApplication(int argc, char* argv[]);
};

// 运行已链接源文件注册的 Solver；统一管理参数、MPI 生命周期、失败退出与最终输出。
int runApplication(int argc, char* argv[]);

}  // babelsim 命名空间
