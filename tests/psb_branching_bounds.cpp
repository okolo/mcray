#include "Nucleus.h"
#include <gtest/gtest.h>
#include <climits>

namespace {
using namespace Interactions;

TIsotope makeIsotope(const TBranching& branching) {
    return {56, 26, {{10., 20., 1., 5., 5.}, {10., 20., 1., 5., 5.}},
            1., &branching};
}

class PsbBranchingTest : public ::testing::TestWithParam<int> {};

TEST_P(PsbBranchingTest, InvalidDeltaAReturnsZero) {
    const double coefficients[] = {0.1, 0.3, 0.6};
    const TBranching branching = {GetParam(), coefficients, -1., -1.};
    const TIsotope isotope = makeIsotope(branching);
    CPhotoDisintegrationCanal channel(isotope, EMultiple);
    for (int delta : {INT_MIN, -1, 0, GetParam() + 1, INT_MAX})
        EXPECT_EQ(channel.getRate(delta), 0.) << "deltaA=" << delta;
}

TEST_P(PsbBranchingTest, ValidCoefficientsAreUnchanged) {
    const double coefficients[] = {0.1, 0.3, 0.6};
    const TBranching branching = {GetParam(), coefficients, -1., -1.};
    const TIsotope isotope = makeIsotope(branching);
    CPhotoDisintegrationCanal channel(isotope, EMultiple);
    for (int delta = 1; delta <= GetParam(); ++delta)
        EXPECT_EQ(channel.getRate(delta), coefficients[delta - 1]) << "deltaA=" << delta;
}

INSTANTIATE_TEST_SUITE_P(ModeCounts, PsbBranchingTest, ::testing::Values(1, 2, 3));

TEST(PsbLegacyChannels, SingleNucleonChannelIsUnchanged) {
    const double coefficients[] = {1.};
    const TBranching branching = {1, coefficients, -1., -1.};
    const TIsotope isotope = makeIsotope(branching);
    CPhotoDisintegrationCanal channel(isotope, ESingle);
    for (int delta = -1; delta <= 4; ++delta)
        EXPECT_EQ(channel.getRate(delta), delta == 1 ? 1. : 0.);
}

TEST(PsbLegacyChannels, DoubleNucleonChannelIsUnchanged) {
    const double coefficients[] = {1.};
    const TBranching branching = {1, coefficients, -1., -1.};
    const TIsotope isotope = makeIsotope(branching);
    CPhotoDisintegrationCanal channel(isotope, EDouble);
    for (int delta = -1; delta <= 4; ++delta)
        EXPECT_EQ(channel.getRate(delta), delta == 2 ? 1. : 0.);
}
} // namespace
