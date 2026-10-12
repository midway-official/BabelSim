#pragma once

#include "babelsim/config.h"
#include "babelsim/equation_control.h"
#include "babelsim/field.h"
#include "babelsim/time.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"
#include "babelsim/model_description.h"

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
// 算例拥有网格和所有命名场，并管理本次运行的时间/结果文件。场引用在 Case 生存期内
// 稳定；同名场只读一次。算法显式管理 time::History，不需要管理执行对象或并行 I/O。
class Case {
public:
    // All named fields belong to this Case. Bind every load/create/fields result
    // with auto&; ordinary Field value copies remain independent snapshots.
    class Load {
        friend class Case;
        explicit Load(Case& owner) : owner_(owner) {}
        Case& owner_;
    public:
        Load(const Load&) = delete;
        Load& operator=(const Load&) = delete;
        // Load values and physical boundaries once; repeated loads reuse the field.
        ScalarField& scalarField(const std::string& name);
        VectorField& vectorField(const std::string& name);
        TensorField& tensorField(const std::string& name);
    };
    class Create {
        friend class Case;
        explicit Create(Case& owner) : owner_(owner) {}
        Case& owner_;
    public:
        Create(const Create&) = delete;
        Create& operator=(const Create&) = delete;
        // Create once, without file I/O. Initial values default to zero.
        // New cell fields start with zeroGradient boundaries.
        ScalarField& scalarField(const std::string& name, double initial = 0.0);
        VectorField& vectorField(const std::string& name, Vec3 initial = {});
        TensorField& tensorField(const std::string& name, Tensor3 initial = {});
        ScalarField& faceScalarField(const std::string& name, double initial = 0.0);
        VectorField& faceVectorField(const std::string& name, Vec3 initial = {});
        TensorField& faceTensorField(const std::string& name, Tensor3 initial = {});
        // Zero cell field, with the existing homogeneous correction boundary
        // mapping. Registered and owned exactly like every other created field.
        ScalarField& homogeneousCorrection(const ScalarField& original, const std::string& name);
    };
    class Fields {
        friend class Case;
        explicit Fields(Case& owner) : owner_(owner) {}
        Case& owner_;
    public:
        Fields(const Fields&) = delete;
        Fields& operator=(const Fields&) = delete;
        // Retrieve an already declared field, regardless of its initialization
        // source. Never read a file, create an object, or reset its values.
        ScalarField& scalarField(const std::string& name);
        VectorField& vectorField(const std::string& name);
        TensorField& tensorField(const std::string& name);
        ScalarField& faceScalarField(const std::string& name);
        VectorField& faceVectorField(const std::string& name);
        TensorField& faceTensorField(const std::string& name);
    };
    class Read {
        friend class Case;
        explicit Read(const Case& owner) : owner_(owner) {}
        const Case& owner_;
    public:
        Read(const Read&) = delete;
        Read& operator=(const Read&) = delete;
        // Configuration only: no field creation or time advancement. Dictionaries
        // and methods are loaded once; parameters applies the module description.
        const Parameters& physics() const;
        const Parameters& solution() const;
        ModuleParameters parameters(const ModuleDescription&) const;
        template<class Config> Config config(const ConfigurationDescription<Config>& description) const {
            return description.configuration(parameters(description));
        }
        const Methods& methods() const;
        const TimeControl& timeControl() const;
        const OutputControl& outputControl() const;
    };

    Load load;
    Create create;
    Fields fields;
    Read read;

    explicit Case(const std::filesystem::path& directory, const std::string& run_name = {});
    ~Case();
    Case(const Case&) = delete;
    Case& operator=(const Case&) = delete;

    const std::string& solver() const;
    const CaseDefinition& definition() const;
    const Mesh& mesh() const;
    // 选择 write() 将写出的 Case 自有单元场；可添加派生量或关闭某个输入场。
    // 只改变后续写出选择，不立即写文件。当前结果格式只保存单元场。
    template <typename T>
    void output(const Field<T>& field, bool enabled = true) {
        selectOutput(field.name(), &field, enabled);
    }

    // Procedural lifecycle: these calls never advance field histories.
    void setTime(double value, int step, double dt);
    void write();
    double time() const;
    int step() const;
    // Validate declared modules, fields, equation bindings and consumed config;
    // then close declarations. Does not advance time, save history or solve.
    void finalizeInitialization();
    void finish();
    // Data only. Solver code decides when and how to print these counters.
    // performance() is the critical-path (maximum rank) view; localPerformance()
    // is the rank-local view used by benchmark aggregation.
    PerformanceCounters performance() const;
    PerformanceCounters localPerformance() const;

private:
    template<class T, class Terms> EquationControl bindEquationControl(const std::string&,
        const Field<T>&, const Terms&, const std::string& expectedUnknown = {}) const;
    void bindEquationIdentity(const std::string&, const void*, bool, const std::vector<std::string>&) const;
    friend EquationControl readEquationControl(const Case&, const std::string&, const ScalarField&,
        std::initializer_list<std::string>);
    friend EquationControl readEquationControl(const Case&, const std::string&, const VectorField&,
        std::initializer_list<std::string>);
    friend EquationControl readEquationControl(const Case&, const EquationDescription&, const ScalarField&);
    friend EquationControl readEquationControl(const Case&, const EquationDescription&, const VectorField&);
    void selectOutput(const std::string& name, const void* field, bool enabled);
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};

// Read-only configuration values and explicit I/O; no algorithm lifecycle.
struct TimeOptions { double start, end, dt; };
TimeOptions readTimeControl(const Case&);
int readWriteInterval(const Case&);
void write(Case&, double time, int step);
void write(Case&, const TimeStepper&);

}  // babelsim 命名空间
