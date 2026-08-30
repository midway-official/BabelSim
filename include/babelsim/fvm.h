#pragma once

#include "babelsim/field.h"
#include "babelsim/fvc.h"

#include <utility>
#include <vector>

namespace babelsim {

// 这些轻量描述符是连续/离散数学之间的边界：它们只引用 Field 和常数，绝不分配
// LDU、CSR、Eigen 向量或 MPI 缓冲区。真正组装只发生在 solve() 内部。
enum class FvmTermKind {
    TimeDerivative,
    Convection,
    Laplacian,
    Source,
    Gradient,
};

struct ScalarFvmTerm {
    FvmTermKind kind = FvmTermKind::Source;
    int sign = 1;
    double coefficient = 1.0;
    const ScalarField* field = nullptr;
    const ScalarField* coefficient_field = nullptr;
    const ScalarField* flux = nullptr;
};

struct VectorFvmTerm {
    FvmTermKind kind = FvmTermKind::Source;
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
    explicit ScalarExpression(ScalarFvmTerm term);

private:
    std::vector<ScalarFvmTerm> m_terms;

    ScalarExpression& add(const ScalarExpression& other, int sign);
    ScalarExpression& negate();

    friend ScalarExpression operator+(
        ScalarExpression left, const ScalarExpression& right);
    friend ScalarExpression operator-(
        ScalarExpression left, const ScalarExpression& right);
    friend ScalarExpression operator-(ScalarExpression expression);
    friend class RunTime;
    friend struct ScalarEquationDefinition;
};

class VectorExpression {
public:
    VectorExpression() = default;
    explicit VectorExpression(VectorFvmTerm term);
    explicit VectorExpression(fvc::ScalarGradient gradient);

private:
    std::vector<VectorFvmTerm> m_terms;

    VectorExpression& add(const VectorExpression& other, int sign);
    VectorExpression& negate();

    friend VectorExpression operator+(
        VectorExpression left, const VectorExpression& right);
    friend VectorExpression operator-(
        VectorExpression left, const VectorExpression& right);
    friend VectorExpression operator-(VectorExpression expression);
    friend class RunTime;
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

namespace fvm {

inline ScalarExpression ddt(double density, const ScalarField& field) {
    return ScalarExpression({FvmTermKind::TimeDerivative, 1, density, &field});
}

inline ScalarExpression ddt(
    const ScalarField& volumetric_capacity,
    const ScalarField& field)
{
    return ScalarExpression(
        {FvmTermKind::TimeDerivative, 1, 1.0, &field, &volumetric_capacity});
}

inline ScalarExpression ddt(const ScalarField& field) { return ddt(1.0, field); }

inline ScalarExpression div(
    double flux_scale,
    const ScalarField& flux,
    const ScalarField& field)
{
    return ScalarExpression(
        {FvmTermKind::Convection, 1, flux_scale, &field, nullptr, &flux});
}

inline ScalarExpression div(const ScalarField& flux, const ScalarField& field) {
    return div(1.0, flux, field);
}

inline ScalarExpression laplacian(double diffusivity, const ScalarField& field) {
    return ScalarExpression({FvmTermKind::Laplacian, 1, diffusivity, &field});
}

inline ScalarExpression laplacian(
    const ScalarField& diffusivity, const ScalarField& field)
{
    return ScalarExpression(
        {FvmTermKind::Laplacian, 1, 1.0, &field, &diffusivity});
}

inline ScalarExpression source(double value) {
    return ScalarExpression({FvmTermKind::Source, 1, value});
}

inline ScalarExpression source(const ScalarField& field) {
    return ScalarExpression({FvmTermKind::Source, 1, 1.0, &field});
}

inline VectorExpression ddt(double density, const VectorField& field) {
    VectorFvmTerm term;
    term.kind = FvmTermKind::TimeDerivative;
    term.coefficient = density;
    term.vector_field = &field;
    return VectorExpression(term);
}

inline VectorExpression ddt(const ScalarField& density, const VectorField& field) {
    VectorFvmTerm term;
    term.kind = FvmTermKind::TimeDerivative;
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
    VectorFvmTerm term;
    term.kind = FvmTermKind::Convection;
    term.coefficient = flux_scale;
    term.vector_field = &field;
    term.flux = &flux;
    return VectorExpression(term);
}

inline VectorExpression div(const ScalarField& flux, const VectorField& field) {
    return div(1.0, flux, field);
}

inline VectorExpression laplacian(double diffusivity, const VectorField& field) {
    VectorFvmTerm term;
    term.kind = FvmTermKind::Laplacian;
    term.coefficient = diffusivity;
    term.vector_field = &field;
    return VectorExpression(term);
}

inline VectorExpression laplacian(
    const ScalarField& diffusivity, const VectorField& field)
{
    VectorFvmTerm term;
    term.kind = FvmTermKind::Laplacian;
    term.vector_field = &field;
    term.coefficient_field = &diffusivity;
    return VectorExpression(term);
}

inline VectorExpression source(const Vec3& value) {
    VectorFvmTerm term;
    term.kind = FvmTermKind::Source;
    term.vector_source = value;
    return VectorExpression(term);
}

}  // fvm 命名空间

// 均匀体源可以直接写成方程右端的数值；这仍表示体源项，不会创建临时 Field。
inline ScalarExpression operator+(ScalarExpression left, double source) {
    return left + fvm::source(source);
}

inline ScalarExpression operator+(double source, ScalarExpression right) {
    return fvm::source(source) + right;
}

inline ScalarExpression operator-(ScalarExpression left, double source) {
    return left - fvm::source(source);
}

namespace fvc {

inline VectorExpression asExpression(ScalarGradient gradient) {
    return VectorExpression(gradient);
}

inline VectorExpression operator-(ScalarGradient gradient) {
    return -asExpression(gradient);
}

}  // fvc 命名空间

}  // babelsim 命名空间
