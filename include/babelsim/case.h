#pragma once

#include "babelsim/config.h"
#include "babelsim/field.h"
#include "babelsim/time.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"

#include <memory>

#include <filesystem>
#include <string>

namespace babelsim {

// 与求解器无关的案例入口。物理专属字典分离存放，使启动器能先选择求解器再读取它们。
struct CaseDefinition {
    std::filesystem::path root;
    std::string solver;
    std::filesystem::path mesh_file;
    std::filesystem::path fields_directory;
    std::filesystem::path physics_file;
    std::filesystem::path methods_file;
    std::filesystem::path solution_file;
    std::filesystem::path control_file;
    std::filesystem::path output_file;
    Index ghost_layers = 3;
};

struct OutputControl {
    std::filesystem::path directory = "results";
    std::string time_name = "final";
    int write_interval = 1;
    // Optional names. Empty write_fields means the historical default: all
    // file-loaded cell fields, except exclude_fields.
    std::vector<std::string> write_fields;
    std::vector<std::string> exclude_fields;
};

CaseDefinition readCase(const std::filesystem::path& case_directory);
OutputControl readOutputControl(const CaseDefinition& definition);
// 算例拥有网格和已读取的场，并管理本次运行的时间/结果文件。场引用在 Case 生存期内
// 稳定；同名场只读一次。算法显式管理 time::History，不需要管理执行对象或并行 I/O。
class Case {
public:
    explicit Case(const std::filesystem::path& directory, const std::string& run_name = {});
    ~Case();
    Case(const Case&) = delete;
    Case& operator=(const Case&) = delete;

    const std::string& solver() const;
    const Mesh& mesh() const;
    // 读取 case.bs 的 physics 条目所指向的字典；solution() 对应 solution 条目。
    const Parameters& physics() const;
    const Parameters& solution() const;
    // 文件加载的单元场默认选入输出列表；真正写出仍需显式 write()。
    // These overloads load the named cell field from fields/initial. Repeated
    // loads return the same object; a program-created field cannot be loaded later.
    ScalarField& scalarField(const std::string& name);
    VectorField& vectorField(const std::string& name);
    // Access an already declared field without loading a file or creating a
    // second object. This is useful for program-owned derived quantities.
    ScalarField& existingScalarField(const std::string& name);
    VectorField& existingVectorField(const std::string& name);
    // Program-owned fields never read files. They are output only when listed in
    // output.bs or explicitly selected with output(field).
    ScalarField& createScalarField(const std::string& name, double initial = 0.0);
    VectorField& createVectorField(const std::string& name, Vec3 initial = {});
    TensorField& tensorField(const std::string& name);
    TensorField& existingTensorField(const std::string& name);
    TensorField& createTensorField(const std::string& name, Tensor3 initial = {});
    // Explicitly names the value type and location.  createFaceField is kept
    // as a source-compatible alias for older solvers.
    ScalarField& createFaceScalarField(const std::string& name);
    ScalarField& createFaceField(const std::string& name);
    VectorField& createFaceVectorField(const std::string& name);
    TensorField& createFaceTensorField(const std::string& name);
    ScalarField& existingFaceField(const std::string& name);
    VectorField& existingFaceVectorField(const std::string& name);
    TensorField& existingFaceTensorField(const std::string& name);

    // 选择 write() 将写出的 Case 自有单元场；可添加派生量或关闭某个输入场。
    // 只改变后续写出选择，不立即写文件。当前结果格式只保存单元场。
    template <typename T>
    void output(const Field<T>& field, bool enabled = true) {
        selectOutput(field.name(), &field, enabled);
    }

    // 下一次 loop() 前保存已完成时间步；自然退出时保证最终时刻写出。
    // 若求解失败，请提前返回，不调用 finish()，以免把失败步标成完整结果。
    bool loop();
    const TimeControl& timeControl() const;
    // Numerical methods are parsed once while the Case/runtime is constructed.
    // This accessor is a read-only view; it never reloads the methods file.
    const Methods& methods() const;
    LinearSolverConfig linearControl(bool vector) const;
    // Backward-compatible name for methods(); retained as an idempotent query.
    const Methods& loadMethods();
    const OutputControl& outputControl() const;
    // Procedural lifecycle: these calls never advance field histories.
    void setTime(double value, int step, double dt);
    void write();
    double time() const;
    int step() const;
    // validate 只校验；只有最外层 start/loop 关闭声明阶段，算法构造不改变 Case 状态。
    void validate() const;
    void start();
    void finish();
    // Data only. Solver code decides when and how to print these counters.
    PerformanceCounters performance() const;

private:
    void selectOutput(const std::string& name, const void* field, bool enabled);
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};

// Read-only configuration values and explicit I/O; no algorithm lifecycle.
struct TimeOptions { double start, end, dt; };
const Methods& loadMethods(Case&);
TimeOptions readTimeControl(const Case&);
LinearSolverConfig readLinearControl(const Case&, const ScalarField&);
LinearSolverConfig readLinearControl(const Case&, const VectorField&);
int readWriteInterval(const Case&);
void setTime(Case&, double time);
void write(Case&, double time, int step);
void write(Case&, const TimeStepper&);

}  // babelsim 命名空间
