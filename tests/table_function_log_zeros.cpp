#include "TableFunction.h"
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <memory>

// Exercise the production interpolator for both supported floating-point types.
template<typename X>
class LogInterpolationTest : public ::testing::Test {
protected:
    const std::vector<X> x = {1, 10, 100, 1000, 10000, 100000};
    const std::vector<X> y = {0, 0, 4, 16, 0, 0};
    Utils::LinearFuncX<X> function{x, y, new Utils::LogScaleX<X>(),
        new Utils::LogScaleX<X>(), X(7), X(9), true};
};

using FloatingTypes = ::testing::Types<double, long double>;
TYPED_TEST_SUITE(LogInterpolationTest, FloatingTypes);

TYPED_TEST(LogInterpolationTest, ZeroNodesAndIntervalsAreExactlyZero) {
    for (TypeParam value : {TypeParam(1), TypeParam(3), TypeParam(10),
                           TypeParam(10000), TypeParam(30000), TypeParam(100000)})
        EXPECT_EQ(this->function.f(value), TypeParam(0)) << "X=" << value;
}

TYPED_TEST(LogInterpolationTest, ZeroTransitionsAreLinearInSigmaAgainstLogK) {
    EXPECT_LT(std::abs(this->function.f(std::sqrt(TypeParam(1000))) - TypeParam(2)),
              TypeParam(1e-12));
    EXPECT_LT(std::abs(this->function.f(std::sqrt(TypeParam(10000000))) - TypeParam(8)),
              TypeParam(1e-12));
}

TYPED_TEST(LogInterpolationTest, PositiveIntervalsRetainLogLogInterpolation) {
    EXPECT_LT(std::abs(this->function.f(std::sqrt(TypeParam(100000))) - TypeParam(8)),
              TypeParam(1e-12));
}

TYPED_TEST(LogInterpolationTest, ExternalValuesRetainConfiguredLimits) {
    EXPECT_EQ(this->function.f(TypeParam(0.5)), TypeParam(7));
    EXPECT_EQ(this->function.f(TypeParam(200000)), TypeParam(9));
}

TYPED_TEST(LogInterpolationTest, ClonePreservesZeroHandling) {
    std::unique_ptr<Utils::FunctionX<TypeParam>> clone(this->function.Clone());
    EXPECT_EQ(clone->f(TypeParam(10)), TypeParam(0));
    EXPECT_LT(std::abs(clone->f(std::sqrt(TypeParam(1000))) - TypeParam(2)),
              TypeParam(1e-12));
}

TYPED_TEST(LogInterpolationTest, LegacyCallsRetainTransformedZeroBehavior) {
    Utils::LinearFuncX<TypeParam> legacy(this->x, this->y,
        new Utils::LogScaleX<TypeParam>(), new Utils::LogScaleX<TypeParam>());
    EXPECT_LT(legacy.f(std::sqrt(TypeParam(1000))), TypeParam(1e-100));
}

TYPED_TEST(LogInterpolationTest, OrdinaryLinearInterpolationIsUnchanged) {
    Utils::LinearFuncX<TypeParam> ordinary(this->x, this->y);
    EXPECT_LE(std::abs(ordinary.f(TypeParam(55)) - TypeParam(2)),
              TypeParam(8) * std::numeric_limits<TypeParam>::epsilon() * TypeParam(2));
}

TYPED_TEST(LogInterpolationTest, ZeroHandlingRejectsIncompatibleScales) {
    bool rejected = false;
    try {
        Utils::LinearFuncX<TypeParam> invalid(this->x, this->y,
            new Utils::LogScaleX<TypeParam>(), nullptr, TypeParam(0), TypeParam(0), true);
    } catch (Utils::Exception* error) {
        delete error;
        rejected = true;
    }
    EXPECT_TRUE(rejected);
}
