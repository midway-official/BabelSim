#include "babelsim/monitor.h"
#include "babelsim/solver.h"
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace babelsim::monitor {
Reporter::Reporter(const char* name, int interval) : name_(name ? name : ""), interval_(interval) {
    if (!name || interval <= 0) throw std::invalid_argument("invalid reporter configuration");
}
namespace {
void printMetrics(std::initializer_list<Metric> metrics) {
    for (const auto& metric : metrics) {
        std::cout << ' ' << metric.name << '=';
        std::visit([](auto value) {
            if constexpr (std::is_same_v<decltype(value), bool>)
                std::cout << (value ? "true" : "false");
            else std::cout << value;
        }, metric.value);
    }
    std::cout << '\n';
}
}
void Reporter::record(std::initializer_list<Metric> metrics) const {
    if (!primaryProcess()) return;
    std::cout << name_;
    printMetrics(metrics);
}
void Reporter::iteration(int index, int limit, bool urgent, std::initializer_list<Metric> metrics) const {
    if (!primaryProcess() || !(index == 1 || index % interval_ == 0 || urgent || index == limit)) return;
    std::cout << name_ << ' ' << index;
    printMetrics(metrics);
}
} // namespace babelsim::monitor
