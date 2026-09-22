#pragma once

#include "babelsim/field.h"
#include "babelsim/methods.h"
#include "babelsim/solver_control.h"
#include <initializer_list>
#include <map>

namespace babelsim {
class Case;

// Resolved numerical data only: no dictionary, execution context or backend.
struct EquationControl {
    std::string name;
    OperatorOptions spatial;
    std::map<std::string, OperatorOptions> terms;
    LinearSolverConfig linear;

    void requireSpatial(const std::string& context) const { spatial.requireComplete(context); }

    OperatorOptions options(const std::string& term = {}) const {
        if (term.empty()) return spatial;
        const auto found = terms.find(term);
        if (found == terms.end()) throw std::invalid_argument("undeclared equation term " + name + "." + term);
        return found->second;
    }
};

EquationControl readEquationControl(const Case&, const std::string&, const ScalarField&,
    std::initializer_list<std::string> terms = {});
EquationControl readEquationControl(const Case&, const std::string&, const VectorField&,
    std::initializer_list<std::string> terms = {});
} // namespace babelsim
