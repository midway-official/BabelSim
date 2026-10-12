# 模块描述、参数与算例模板

[开发指南](README.md) · [新增求解器](extending-solvers.md) · [时间推进](time-integration.md)

`ModuleDescription` 是一份公开的 C++ 数据描述。它声明字段、单位、参数类型和范围、方程及离散项、时间格式、边界能力和组合限制。参数解析、前端校验、Markdown 文档、JSON Schema 和小型算例模板共同读取这份描述。物理算法仍由求解器显式编写，`equ::ddt/laplacian/source` 等离散接口保持原有语义。

## 直接查看内置模块

```bash
build/babelsim-solve -list
build/babelsim-solve -describe heat
build/babelsim-solve -describe piso -format json
build/babelsim-solve -describe SA
build/babelsim-solve -describe kOmega
build/babelsim-solve -describe kEpsilon
```

`-describe` 只输出描述，不构造 Case 或推进时间。JSON 输出适合由编辑器、配置工具和后续 GUI 消费；当前实现没有增加 GUI 页面。

heat、transport 提供完整的单单元立方体教程模板：

```bash
build/babelsim-solve -template heat -output work/describedHeat
build/babelsim-solve -case work/describedHeat -time first
build/babelsim-solve -template transport -output work/describedTransport
build/babelsim-solve -case work/describedTransport -time first
```

生成目录包含合法网格、输入场、physics、methods、solution、control、output、README 和 `module.schema.json`。描述中的必需参数用 `template_value` 写入教程；这不会把它们变成运行时可省略参数。methods 和 solution 始终生成现有底层要求的完整逐方程配置，包括纯热方程的 convection 选项。生成器拒绝覆盖已有目录。

教程默认零梯度边界与均匀源项；transport 的教程速度为零，使结果可用 `C=S*t/storage` 核对。增加对流时，应同时配置速度和入口/出口边界。教程适合检查 SDK 和时间推进流程，不代表空间精度或工程物理验证。流动与 RANS 模型的边界/耦合需要更具体的设计，因此目前只提供描述，不通过通用模板命令生成流动验证算例。

## 给新求解器添加描述

描述由求解器注册对象引用，必须在注册对象存活期间保持有效。内置描述通过 `builtin_descriptions.h` 提供，可复用后修改模块名称：

```cpp
#include "babelsim/application.h"
#include "babelsim/builtin_descriptions.h"
#include "babelsim/case.h"

const babelsim::ModuleDescription& description() {
    static const babelsim::ModuleDescription result = [] {
        auto value = babelsim::builtin::heat();
        value.name = "customHeat";
        value.explanation = "My heat-conduction solver.";
        return value;
    }();
    return result;
}

babelsim::SolverResult runCustomHeat(babelsim::Case& problem) {
    const auto parameters = problem.parameters(description());
    const double rho = parameters.physics.number("density");
    const double cp = parameters.physics.number("heatCapacity");
    const double k = parameters.physics.number("conductivity");
    const double Q = parameters.physics.number("source");
    // Declare fields/history/equations, validate actual numerical use,
    // and explicitly assemble rho*cp*ddt(T) - laplacian(k,T) = Q.
    // See external_heat.cpp for the complete algorithm.
    (void)rho; (void)cp; (void)k; (void)Q;
    return babelsim::SolverResult::completed();
}

const babelsim::SolverRegistration registration(
    "customHeat", runCustomHeat, &description());
```

上述短片段说明接口连接；[external_heat.cpp](../examples/external_heat.cpp)提供完整可编译求解算法。外部程序沿用同一 `runApplication`，即可使用 `-describe customHeat` 和 `-template customHeat -output ...`。新模块统一提供描述；不需要维护单独的模板和参数校验实现。

从头声明一个新模块时，可在一个函数局部静态对象中集中填写：

```cpp
babelsim::ModuleDescription value;
value.name = "myTransport";
value.explanation = "Scalar diffusion with constant capacity.";
value.fields = {{"C", babelsim::DescribedFieldType::Scalar,
    babelsim::DescribedFieldLocation::Cell, babelsim::FieldRole::Input,
    "kg/m^3", "Concentration.", "0"}};

babelsim::ParameterDescription diffusivity;
diffusivity.name = "diffusivity";
diffusivity.type = babelsim::ParameterType::Number;
diffusivity.unit = "m^2/s";
diffusivity.explanation = "Scalar diffusion coefficient.";
diffusivity.minimum = 0.0;
diffusivity.required = false;
diffusivity.default_value = 0.01;
value.parameters = {diffusivity};
value.equations = {{"concentration", "C", {"diffusion"},
    "Unit-capacity transient concentration equation."}};
value.time_methods = {babelsim::TimeMethod::Euler, babelsim::TimeMethod::BDF2};
value.boundary_types = {"fixedValue", "fixedGradient", "zeroGradient", "symmetry"};
value.supports_case_template = true;
```

字符串和列表采用标准 C++ 类型；description 不持有 Case、Field、Equation 或后端对象。字段 `Input` 表示需要用户提供场文件，`Derived` 表示程序创建或宿主提供。当前 Case 场文件初始化路径支持 cell 输入，face 通量通常声明为 Derived。

## 参数类型、默认值与组合约束

`problem.parameters(description)` 返回 `ModuleParameters`，其中 `physics` 和 `solution` 各有 `number/integer/boolean/word/contains`。解析通过现有 Parameters 完成，继承文件与行号诊断。字段类型和参数类型独立；整数参数必须使用 `int` 默认值，普通数值参数使用 `double` 默认值。

| 描述属性 | 语义 |
| --- | --- |
| required=true | 字典必须提供该参数；不允许同时声明运行时默认值 |
| default_value | 参数缺省时使用并校验的类型化值 |
| template_value | 写入教程算例的示范值；不改变 required 契约 |
| minimum / maximum | 数值与整数的范围 |
| exclusive_minimum / exclusive_maximum | 是否排除范围端点 |
| choices | word 参数的允许值列表 |
| default_factory | 从同一字典中先声明的已解析参数计算默认值 |
| default_explanation | 依赖默认值的公式说明，写入文档与 JSON 元数据 |

依赖默认值按参数声明顺序计算。SA 的 `saCw1` 从已解析的 `saCb1/saCb2/saKappa/saSigma` 计算，所以用户改变这些常数时，缺省 cw1 仍保持原来的依赖语义。PISO 的 maxCorrectors 和初始非正交修正次数也采用此机制。依赖模型选择和算法状态的 turbulenceCorrectors 则保持可选且无描述默认值，由宿主显式计算。

跨参数数值约束采用结构化声明：

```cpp
value.parameter_constraints.push_back({
    babelsim::ParameterDictionary::Solution,
    "maxCorrectors", babelsim::ParameterRelation::GreaterEqual,
    "nCorrectors", "The maximum must cover the requested minimum."});
```

两个有效值都存在时，参数解析执行该规则。`time_methods` 和 `compatible_solvers` 在前端校验中执行；`limitations` 保存需要开发者理解的说明。复杂的条件耦合规则继续由物理算法显式检查。复用 RANS 描述到自定义宿主时，应明确扩展 `compatible_solvers`，并满足该模型的场、时间和耦合契约。

## 前端校验与现有严格配置检查

带描述的 SolverRegistration 在调用求解函数前执行 `validateModuleConfiguration`：

1. 解析声明参数，检查类型、范围、choices 和跨参数数值关系。
2. 检查时间格式与宿主兼容性。
3. 只读检查 Input 场文件的名称、类型、位置及声明边界能力；使用统一边界名称。
4. 检查方程的完整 methods profile 与必需线性求解配置是否存在。

此过程不加载或创建字段，不自动声明派生场，不改变时间，也不把数值 profile 标记为 used。求解器随后仍需创建自身字段、历史和方程，并调用 `problem.validate()`。参数描述会消费其声明的配置项；描述之外的未知键与实际未使用的数值 profile 仍由原检查拒绝。

`moduleJson()` 对归一化 physics/solution 对象生成标准 JSON Schema；原生 `.bs` 仍通过现有行式解析器解析。Schema 保留 `additionalProperties=true`，因为宿主与多个模型共享字典；`.bs` 的严格未知/未用检查以 Case 为准。字段、方程、时间和组合信息保存在 `x-*` 元数据中，适合生成表单与文档。

组合模型可用 metadata-only 注册进行工具发现：

```cpp
const babelsim::ModuleRegistration modelDescription(description());
```

它只发布描述，不建立新的物理模型工厂或动态插件 ABI。RANS 三模型已通过这一路径公开描述。

## 验证

```bash
make test-descriptions
python3 tests/module_description_workflow_test.py \
    --solver build/babelsim-solve --external build/examples/custom-heat
```

C++ 测试覆盖类型、依赖默认值、参数范围、时间/宿主/字段/边界检查、数值配置不被提前消费、未知键拒绝和注册生命周期。工作流测试解析全部内置 JSON/Markdown 描述，生成并运行 heat/transport 的 Euler/BDF2 教程，检查解析结果与拒绝覆盖；提供 external 参数时也检查外部求解器同一路径。
