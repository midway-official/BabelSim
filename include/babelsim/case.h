#pragma once

#include "babelsim/config.h"
#include "babelsim/field.h"

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
};

struct OutputControl {
    std::filesystem::path directory = "results";
    std::string time_name = "final";
    int write_interval = 1;
};

CaseDefinition readCase(const std::filesystem::path& case_directory);
OutputControl readOutputControl(const CaseDefinition& definition);
// 算例拥有网格和已读取的场，并管理本次运行的时间/结果文件。场引用在 Case 生存期内
// 稳定；同名场只读一次。方程或算法不需要创建执行对象、管理历史场或调用并行 I/O。
class Case {
public:
    explicit Case(const std::filesystem::path& directory, const std::string& run_name = {});
    ~Case();
    Case(const Case&) = delete;
    Case& operator=(const Case&) = delete;

    const std::string& solver() const;
    const Mesh& mesh() const;
    const Parameters& properties() const;
    const Parameters& solution() const;
    ScalarField& scalarField(const std::string& name);
    VectorField& vectorField(const std::string& name);
    // 算法需要的中间数学场：不读取文件，不自动输出，默认零梯度边界。
    // 与输入场共用相同 Field/生命周期机制，无需为新算法引入工作区类。
    ScalarField& scalarField(const std::string& name, double initial);
    VectorField& vectorField(const std::string& name, Vec3 initial);
    TensorField& tensorField(const std::string& name);
    ScalarField& faceField(const std::string& name);
    ScalarField& faceFlux(const std::string& name, const VectorField& velocity);

    // 下一次 loop() 前保存已完成时间步；自然退出时保证最终时刻写出。
    // 若求解失败，请提前返回，不调用 finish()，以免把失败步标成完整结果。
    bool loop();
    double time() const;
    int step() const;
    void validate() const;
    void finish();

private:
    struct Implementation;
    std::unique_ptr<Implementation> m_implementation;
};

}  // babelsim 命名空间
