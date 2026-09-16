#pragma once

#include <initializer_list>
#include <string>
#include <string_view>
#include <variant>

namespace babelsim::monitor {

// Semantic labels and observations, consumed synchronously (never retained).
// No field, equation, tolerance, or solver-specific concept belongs here.
struct Metric {
    std::string_view name;
    std::variant<double, int, bool, std::string_view> value;
    Metric(std::string_view label, double v) : name(label), value(v) {}
    Metric(std::string_view label, int v) : name(label), value(v) {}
    Metric(std::string_view label, bool v) : name(label), value(v) {}
    Metric(std::string_view label, const char* v) : name(label), value(std::string_view(v)) {}
    Metric(std::string_view label, std::string_view v) : name(label), value(v) {}
};

class Reporter {
public:
    explicit Reporter(const char* name, int iterationInterval = 100);
    void record(std::initializer_list<Metric>) const;
    // First, periodic and final iterations are reported. urgent is decided by
    // the caller, e.g. for convergence or a warning; it never drives control flow.
    void iteration(int index, int limit, bool urgent, std::initializer_list<Metric>) const;
private:
    std::string name_;
    int interval_;
};
} // namespace babelsim::monitor
