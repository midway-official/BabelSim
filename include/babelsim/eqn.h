#pragma once

#include "babelsim/field.h"
#include "babelsim/math.h"

#include <utility>
#include <vector>

namespace babelsim {
namespace detail { class FvmExecution; }

// 这些轻量描述符是连续/离散数学之间的边界：它们只引用 Field 和常数，绝不分配
// LDU、CSR、Eigen 向量或 MPI 缓冲区。真正组装只发生在 solve() 内部。
enum class EquationTermKind {
    TimeDerivative,
    Convection,
    Laplacian,
    Source,
    Gradient,
};

struct ScalarEquationTerm {
    EquationTermKind kind = EquationTermKind::Source;
    int sign = 1;
    double coefficient = 1.0;
    const ScalarField* field = nullptr;
    const ScalarField* coefficient_field = nullptr;
    const ScalarField* flux = nullptr;
};

struct VectorEquationTerm {
    EquationTermKind kind = EquationTermKind::Source;
    int sign = 1;
    double coefficient = 1.0;
    const VectorField* vector_field = nullptr;
    const ScalarField* scalar_field = nullptr;
    const ScalarField* coefficient_field = nullptr;
    const ScalarField* flux = nullptr;
    Vec3 vector_source{};
};

class ScalarExpression {
public:
    ScalarExpression() = default;
    explicit ScalarExpression(ScalarEquationTerm term);

private:
    std::vector<ScalarEquationTerm> m_terms;

    ScalarExpression& add(const ScalarExpression& other, int sign);
    ScalarExpression& negate();

    friend ScalarExpression operator+(
        ScalarExpression left, const ScalarExpression& right);
    friend ScalarExpression operator-(
        ScalarExpression left, const ScalarExpression& right);
    friend ScalarExpression operator-(ScalarExpression expression);
    friend class detail::FvmExecution;
    friend struct ScalarEquationDefinition;
};

class VectorExpression {
public:
    VectorExpression() = default;
    explicit VectorExpression(VectorEquationTerm term);
    explicit VectorExpression(math::ScalarGradient gradient);

private:
    std::vector<VectorEquationTerm> m_terms;

    VectorExpression& add(const VectorExpression& other, int sign);
    VectorExpression& negate();

    friend VectorExpression operator+(
        VectorExpression left, const VectorExpression& right);
    friend VectorExpression operator-(
        VectorExpression left, const VectorExpression& right);
    friend VectorExpression operator-(VectorExpression expression);
    friend class detail::FvmExecution;
    friend struct VectorEquationDefinition;
};

inline ScalarExpression operator+(
    ScalarExpression left, const ScalarExpression& right)
{
    return left.add(right, 1);
}

inline ScalarExpression operator-(
    ScalarExpression left, const ScalarExpression& right)
{
    return left.add(right, -1);
}

inline ScalarExpression operator-(ScalarExpression expression) {
    return expression.negate();
}

inline VectorExpression operator+(
    VectorExpression left, const VectorExpression& right)
{
    return left.add(right, 1);
}

inline VectorExpression operator-(
    VectorExpression left, const VectorExpression& right)
{
    return left.add(right, -1);
}

inline VectorExpression operator-(VectorExpression expression) {
    return expression.negate();
}

struct ScalarEquationDefinition {
    ScalarExpression lhs;
    ScalarExpression rhs;
};

struct VectorEquationDefinition {
    VectorExpression lhs;
    VectorExpression rhs;
};

inline ScalarEquationDefinition operator==(
    ScalarExpression lhs, ScalarExpression rhs)
{
    return {std::move(lhs), std::move(rhs)};
}

inline VectorEquationDefinition operator==(
    VectorExpression lhs, VectorExpression rhs)
{
    return {std::move(lhs), std::move(rhs)};
}

namespace eqn {

// 构造方程贡献，不立即装配。输运项隐式进入方程，source 则是已知的显式体源。

inline ScalarExpression ddt(double density, const ScalarField& field) {
    return ScalarExpression({EquationTermKind::TimeDerivative, 1, density, &field});
}

inline ScalarExpression ddt(
    const ScalarField& volumetric_capacity,
    const ScalarField& field)
{
    return ScalarExpression(
        {EquationTermKind::TimeDerivative, 1, 1.0, &field, &volumetric_capacity});
}

inline ScalarExpression ddt(const ScalarField& field) { return ddt(1.0, field); }

inline ScalarExpression div(
    double flux_scale,
    const ScalarField& flux,
    const ScalarField& field)
{
    return ScalarExpression(
        {EquationTermKind::Convection, 1, flux_scale, &field, nullptr, &flux});
}

inline ScalarExpression div(const ScalarField& flux, const ScalarField& field) {
    return div(1.0, flux, field);
}

inline ScalarExpression laplacian(double diffusivity, const ScalarField& field) {
    return ScalarExpression({EquationTermKind::Laplacian, 1, diffusivity, &field});
}

inline ScalarExpression laplacian(
    const ScalarField& diffusivity, const ScalarField& field)
{
    return ScalarExpression(
        {EquationTermKind::Laplacian, 1, 1.0, &field, &diffusivity});
}

inline ScalarExpression source(double value) {
    return ScalarExpression({EquationTermKind::Source, 1, value});
}

inline ScalarExpression source(const ScalarField& field) {
    return ScalarExpression({EquationTermKind::Source, 1, 1.0, &field});
}

// 显式体源 a*field；直接复用源项系数，不创建缩放后的临时 Field。
inline ScalarExpression source(double coefficient, const ScalarField& field) {
    return ScalarExpression({EquationTermKind::Source, 1, coefficient, &field});
}

inline VectorExpression ddt(double density, const VectorField& field) {
    VectorEquationTerm term;
    term.kind = EquationTermKind::TimeDerivative;
    term.coefficient = density;
    term.vector_field = &field;
    return VectorExpression(term);
}

inline VectorExpression ddt(const ScalarField& density, const VectorField& field) {
    VectorEquationTerm term;
    term.kind = EquationTermKind::TimeDerivative;
    term.vector_field = &field;
    term.coefficient_field = &density;
    return VectorExpression(term);
}

inline VectorExpression ddt(const VectorField& field) { return ddt(1.0, field); }

inline VectorExpression div(
    double flux_scale,
    const ScalarField& flux,
    const VectorField& field)
{
    VectorEquationTerm term;
    term.kind = EquationTermKind::Convection;
    term.coefficient = flux_scale;
    term.vector_field = &field;
    term.flux = &flux;
    return VectorExpression(term);
}

inline VectorExpression div(const ScalarField& flux, const VectorField& field) {
    return div(1.0, flux, field);
}

inline VectorExpression laplacian(double diffusivity, const VectorField& field) {
    VectorEquationTerm term;
    term.kind = EquationTermKind::Laplacian;
    term.coefficient = diffusivity;
    term.vector_field = &field;
    return VectorExpression(term);
}

inline VectorExpression laplacian(
    const ScalarField& diffusivity, const VectorField& field)
{
    VectorEquationTerm term;
    term.kind = EquationTermKind::Laplacian;
    term.vector_field = &field;
    term.coefficient_field = &diffusivity;
    return VectorExpression(term);
}

inline VectorExpression source(const Vec3& value) {
    VectorEquationTerm term;
    term.kind = EquationTermKind::Source;
    term.vector_source = value;
    return VectorExpression(term);
}

inline VectorExpression source(double coefficient, const VectorField& field) {
    VectorEquationTerm term;
    term.kind = EquationTermKind::Source;
    term.coefficient = coefficient;
    term.vector_field = &field;
    return VectorExpression(term);
}

inline VectorExpression source(const VectorField& field) { return source(1.0, field); }

}  // eqn 命名空间

inline ScalarEquationDefinition operator==(ScalarExpression lhs, double rhs) {
    return std::move(lhs) == eqn::source(rhs);
}

inline VectorEquationDefinition operator==(VectorExpression lhs, Vec3 rhs) {
    return std::move(lhs) == eqn::source(rhs);
}

// 均匀体源可以直接写成方程右端的数值；这仍表示体源项，不会创建临时 Field。
inline ScalarExpression operator+(ScalarExpression left, double source) {
    return left + eqn::source(source);
}

inline ScalarExpression operator+(double source, ScalarExpression right) {
    return eqn::source(source) + right;
}

inline ScalarExpression operator-(ScalarExpression left, double source) {
    return left - eqn::source(source);
}

namespace math {

inline VectorExpression asExpression(ScalarGradient gradient) {
    return VectorExpression(gradient);
}

inline VectorExpression operator-(ScalarGradient gradient) {
    return -asExpression(gradient);
}

}  // math 命名空间

}  // babelsim 命名空间
