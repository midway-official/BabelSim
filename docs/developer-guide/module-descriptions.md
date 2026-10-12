# 模块描述、强类型配置与算例模板

[开发指南](README.md) · [新增求解器](extending-solvers.md) · [时间推进](time-integration.md)

模块描述集中声明字段、参数、方程、离散项、时间格式、边界能力和组合限制。`ConfigurationDescription<Config>` 在同一份描述上绑定普通 C++ 配置成员，成员类型决定参数类型；解析、范围校验、Markdown、JSON Schema 和算例模板共同使用这些声明。物理公式仍由 C++ 显式编写，已有 equ/math 的积分和符号语义不变。

## 查看描述与生成算例

```bash
build/babelsim-solve -list
build/babelsim-solve -describe heat
build/babelsim-solve -describe piso -format json
build/babelsim-solve -describe kOmega
build/babelsim-solve -template heat -output work/describedHeat
build/babelsim-solve -case work/describedHeat -time first
```

`-describe` 不构造 Case 或推进时间。heat、transport 的通用模板包括单单元立方体网格、输入场、physics、methods、solution、control、output、README 和 module.schema.json。生成器拒绝覆盖已有目录。必需参数的 `template_value` 是示范值，不会将参数变为可省略；methods/solution 生成底层要求的完整配置。

模板默认均匀源和零梯度边界；transport 的初始 U 为零，可检查 `C=S*t/storage`。它们用于验证 SDK 和时间推进流程，不代表空间精度或工程物理验证。流动/RANS 提供描述，通用模板不替代具体流动算例的设计。

## 一份声明同时绑定参数和成员

配置结构体只声明类型，不重复存放默认值、范围或解析代码：

```cpp
struct ScalarConfig {
    struct Material { double diffusivity; } material;
    int corrections;
    std::optional<int> optionalLimit;
};

ConfigurationDescription<ScalarConfig> spec;
spec.name = "myScalar";
spec.explanation = "Unit-capacity scalar diffusion.";
spec.fields = {{"C", DescribedFieldType::Scalar, DescribedFieldLocation::Cell,
    FieldRole::Input, "1", "Concentration.", "0"}};

ParameterDescription diffusion;
diffusion.name = "diffusivity";
diffusion.unit = "m^2/s";
diffusion.explanation = "Nonnegative scalar diffusion coefficient.";
diffusion.required = false;
diffusion.default_value = 0.01;
diffusion.minimum = 0.0;
spec.bind(diffusion, &ScalarConfig::material, &ScalarConfig::Material::diffusivity);

ParameterDescription corrections;
corrections.name = "corrections";
corrections.dictionary = ParameterDictionary::Solution;
corrections.required = false;
corrections.default_value = 1;
corrections.minimum = 1;
corrections.explanation = "Number of explicit correction passes.";
spec.bind(corrections, &ScalarConfig::corrections);

auto limit = corrections;
limit.name = "optionalLimit";
limit.default_value.reset();
limit.explanation = "Optional solver-controlled limit.";
spec.bind(limit, &ScalarConfig::optionalLimit);

spec.equations = {{"concentration", "C", {"diffusion"}, "Scalar equation."}};
spec.time_methods = {TimeMethod::Euler, TimeMethod::BDF2};
spec.boundary_types = {"fixedValue", "fixedGradient", "zeroGradient", "symmetry"};
spec.supports_case_template = true;
```

在求解器中：

```cpp
const auto config = problem.read.config(spec);
auto& C = problem.load.scalarField("C");
auto equation = equ::createEquation(problem, spec.equation("concentration"), C);
time::HistoryGroup histories;
auto& oldC = histories.addHistory(C);
auto clock = time::makeStepper(problem);
problem.finalizeInitialization();
// 循环直接使用 config.material.diffusivity / config.corrections。
```

double/int/bool/string 成员分别决定 number/integer/boolean/word 参数类型。嵌套成员绑定提供分组补全，不改变 `.bs` 配置文件结构。无默认值的可选参数必须绑定 `std::optional<T>`；缺省值保持 unset，由具体算法决定。Config 不持有 Field、Equation 或指向自身成员的引用，复制/移动不会改变方程绑定。

运行时仅在初始化解析和绑定一次。循环读取普通成员，不进行参数键查找或动态解释。所有配置参数应通过 `bind()` 登记；直接向 typed description 的 parameters 添加未绑定项，会在读取强类型配置时被拒绝。底层动态工具仍可使用 `ModuleDescription` 与 `problem.read.parameters(spec)`。

## 复用内置描述和外部 SDK

内置接口返回强类型描述，转换为 `const ModuleDescription&` 后仍能用于注册、模板和前端工具。复用描述时保留其类型：

```cpp
const ConfigurationDescription<builtin::HeatConfig>& description() {
    static const auto result = [] {
        auto spec = builtin::heat();
        spec.name = "customHeat";
        spec.explanation = "My heat-conduction solver.";
        return spec;
    }();
    return result;
}
// SolverRegistration("customHeat", runCustomHeat, &description());
```

注册引用的描述必须比注册对象存活更久；函数局部静态对象满足这一契约。[external_heat.cpp](../examples/external_heat.cpp) 是完整可编译示例，演示强类型参数、描述方程、组持有的历史和可选步骤事务。[coupled_scalar.cpp](../examples/coupled_scalar.cpp) 演示从头声明双标量模块。外部程序使用同一个 runApplication，即可使用 describe/template 命令。

## 字段角色与初始化来源

所有命名场由 Case 持有，不存在需要求解器作者选择的局部/全局命名场。角色和初值来源是描述信息，不能改变数学 Field 的计算语义。

| 属性 | 取值与含义 |
| --- | --- |
| role | Input：模块输入；Output：模型或求解器输出；Work：工作/修正场 |
| initialization | File：从场文件加载；Provided：宿主已经准备；Constant：代码给定常量；Computed：代码计算；HomogeneousCorrection：从物理场建立齐次修正边界 |
| source_field | HomogeneousCorrection 的原始物理场名称 |
| type / location | scalar/vector/tensor 与 cell/face，初始化检查要求匹配 |

例如 RANS 的 U、phi 和 muEffective 由宿主提供；mut 是模型输出。SIMPLE 的 pPrime 是 Case 持有的 Work，来源为 HomogeneousCorrection，source_field 为 p。文档与 JSON 同时输出角色和来源；模板仅为 File 来源生成初值文件。场的输出选择仍由 `problem.output(field)` 或 output.bs 控制。

常量数值、初值公式和物理有效性仍由显式 C++ 初始化与模型检查负责。描述不执行表达式，也不能仅凭“有场值”判断物理初值正确。

## 默认值、依赖与组合约束

| 描述属性 | 语义 |
| --- | --- |
| required | 必需参数必须由字典提供，不能同时有运行时默认值 |
| default_value | 类型正确的缺省值；integer 用 int，number 用 double |
| template_value | 教程示范值，不改变必需性 |
| minimum / maximum / exclusive_minimum / exclusive_maximum | 数值范围及端点 |
| choices | word 的允许值 |
| default_factory / default_explanation | 按声明顺序从同一字典中先前的有效值计算默认值，并给出公式说明 |

SA 的 saCw1 仍从 saCb1/saCb2/saKappa/saSigma 计算；PISO 的 maxCorrectors 和初始非正交修正次数同样复用依赖默认值。turbulenceCorrectors 保持 optional，由宿主按模型与耦合策略计算，不重复声明描述默认值。

`parameter_constraints` 声明两个数值参数的 <=、>= 或 == 关系；双方有效值都存在时检查。`time_methods` 与 `compatible_solvers` 在前端和实际模块读取时检查。复杂的条件耦合规则仍由宿主显式检查，limitations 保存开发者需要理解的限制。复用 RANS 到自定义宿主时，明确扩展 compatible_solvers 并满足场与时间契约。

## 前端检查与初始化完成

注册描述的前端 `validateModuleConfiguration` 检查参数、场文件头/边界及完整数值配置，不加载字段，不消费实际方程的数值选项，不关闭声明阶段。

实际 `problem.read.config/parameters` 登记本次使用的模块。所有字段、方程和物理初值准备好后，调用 `problem.finalizeInitialization()` 检查必需场及类型/位置、来源、修正场来源、方程未知量与离散项、边界能力及未消费配置，然后关闭资源声明。它不保存历史、不推进时间、不求解，也不替代模型的物理初始化检查。没有描述的低层自定义算法仍须显式完成初始化，保留配置消费与资源阶段检查。

验收命令：`make test-descriptions`、`make test-external` 和文档的 `verify_examples.py`。它们分别检查描述/模板/错误诊断、SDK 链接和实际外部程序使用。
