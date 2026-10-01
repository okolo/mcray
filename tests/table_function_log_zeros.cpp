// Standalone regression test using the production TableFunction implementation.
#include "TableFunction.h"
#include <cmath>
#include <memory>
#include <limits>
#include <stdexcept>

template<typename X> void check(bool value, const char* message) {
    if(!value)
        throw std::runtime_error(message);
}

template<typename X> void testLogZeros() {
    using namespace Utils;
    const std::vector<X> x = {1, 10, 100, 1000, 10000, 100000};
    const std::vector<X> y = {0, 0, 4, 16, 0, 0};
    LinearFuncX<X> function(x, y, new LogScaleX<X>(), new LogScaleX<X>(),
                          X(7), X(9), true);
    check<X>(function.f(X(1)) == 0 && function.f(X(10)) == 0,
             "zero endpoints must be exactly zero");
    check<X>(function.f(X(3)) == 0 && function.f(X(30000)) == 0,
             "two zero nodes must define a zero interval");
    check<X>(std::abs(function.f(std::sqrt(X(1000))) - X(2)) < X(1e-12),
             "zero-to-positive interval must be linear in Y versus log(X)");
    check<X>(std::abs(function.f(std::sqrt(X(10000000))) - X(8)) < X(1e-12),
             "positive-to-zero interval must be linear in Y versus log(X)");
    check<X>(std::abs(function.f(std::sqrt(X(100000))) - X(8)) < X(1e-12),
             "positive interval must retain geometric log-log interpolation");
    check<X>(function.f(X(10000)) == 0 && function.f(X(100000)) == 0,
             "trailing zero endpoints must be exactly zero");
    check<X>(function.f(X(0.5)) == 7 && function.f(X(200000)) == 9,
             "out-of-table values must retain the configured limits");
    std::unique_ptr<FunctionX<X>> clone(function.Clone());
    check<X>(clone->f(X(10)) == 0 &&
             std::abs(clone->f(std::sqrt(X(1000))) - X(2)) < X(1e-12),
             "Clone must preserve zero handling and original Y values");

    // Legacy callers still use the original transformed-zero behavior.
    LinearFuncX<X> legacy(x, y, new LogScaleX<X>(), new LogScaleX<X>());
    check<X>(legacy.f(std::sqrt(X(1000))) < X(1e-100),
             "the new parameter must not change legacy interpolation");
    LinearFuncX<X> ordinary(x, y);
    check<X>(std::abs(ordinary.f(X(55)) - X(2)) <=
             X(8) * std::numeric_limits<X>::epsilon() * X(2),
             "ordinary linear interpolation must remain unchanged");
    bool rejected = false;
    try {
        LinearFuncX<X> invalid(x, y, new LogScaleX<X>(), nullptr,
                              X(0), X(0), true);
    } catch(Exception* error) {
        delete error;
        rejected = true;
    }
    check<X>(rejected, "zero-aware log-log mode requires both logarithmic scales");
}

int main() {
    testLogZeros<double>();
    testLogZeros<long double>();
    return 0;
}
