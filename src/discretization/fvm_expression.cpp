#include "babelsim/fvm.h"

#include <utility>

namespace babelsim {

namespace {

template <typename Term>
void append(std::vector<Term>& destination, const std::vector<Term>& source, int sign) {
    destination.reserve(destination.size() + source.size());
    for (Term term : source) {
        term.sign *= sign;
        destination.push_back(term);
    }
}

}  // 匿名命名空间

ScalarExpression::ScalarExpression(ScalarFvmTerm term) {
    m_terms.push_back(term);
}

ScalarExpression& ScalarExpression::add(const ScalarExpression& other, int sign) {
    append(m_terms, other.m_terms, sign);
    return *this;
}

ScalarExpression& ScalarExpression::negate() {
    for (ScalarFvmTerm& term : m_terms) term.sign = -term.sign;
    return *this;
}

VectorExpression::VectorExpression(VectorFvmTerm term) {
    m_terms.push_back(term);
}

VectorExpression::VectorExpression(fvc::ScalarGradient gradient) {
    VectorFvmTerm term;
    term.kind = FvmTermKind::Gradient;
    term.scalar_field = &gradient.field;
    m_terms.push_back(term);
}

VectorExpression& VectorExpression::add(const VectorExpression& other, int sign) {
    append(m_terms, other.m_terms, sign);
    return *this;
}

VectorExpression& VectorExpression::negate() {
    for (VectorFvmTerm& term : m_terms) term.sign = -term.sign;
    return *this;
}

}  // babelsim 命名空间
