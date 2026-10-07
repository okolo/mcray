#include <gtest/gtest.h>
#include "gzk_spectral_statistics.h"

TEST(GzkSpectralStatistics, BinsIncludeInjectionAndSeparateUnsupportedTails) {
    EXPECT_EQ(GzkSpectral::bin(1e13), -1);
    EXPECT_EQ(GzkSpectral::bin(1e14), 0);
    EXPECT_EQ(GzkSpectral::bin(1e17), 1);
    EXPECT_EQ(GzkSpectral::bin(1e18), 2);
    EXPECT_EQ(GzkSpectral::bin(1e19), 3);
    EXPECT_EQ(GzkSpectral::bin(1e20), 3);
    EXPECT_EQ(GzkSpectral::bin(1.01e20), 4);
    EXPECT_THROW(GzkSpectral::bin(0.), std::invalid_argument);
}

TEST(GzkSpectralStatistics, BatchErrorPreservesCorrelatedHistoryProducts) {
    // Two siblings in each occupied history must remain in the same batch.
    const auto result = GzkSpectral::moments({0., 2., 0., 2.});
    EXPECT_DOUBLE_EQ(result.mean, 1.);
    EXPECT_NEAR(result.se, 1./std::sqrt(3.), 1e-14);
    const auto scaled = GzkSpectral::moments({0., 6., 0., 6.});
    EXPECT_NEAR(scaled.se, 3.*result.se, 1e-14);
    EXPECT_THROW(GzkSpectral::moments({1.}), std::invalid_argument);
}

TEST(GzkSpectralStatistics, BonferroniQuantileMatchesAnalyticCauchyTail) {
    // Student t with one degree of freedom is a standard Cauchy distribution.
    const double expected = 1./std::tan(std::acos(-1.)*.01/(2.*28));
    EXPECT_NEAR(GzkSpectral::critical(2), expected, 1e-8*expected);
    EXPECT_GT(GzkSpectral::critical(20), GzkSpectral::critical(20, 1));
}

TEST(GzkSpectralStatistics, EquivalenceRequiresEntireIntervalInsideBudget) {
    using GzkSpectral::Outcome;
    EXPECT_EQ(GzkSpectral::classify(1.02, 1., .01, .005, true, true), Outcome::Agreement);
    EXPECT_EQ(GzkSpectral::classify(1.045, 1., .01, .005, true, true), Outcome::Inconclusive);
    EXPECT_EQ(GzkSpectral::classify(1.08, 1., .01, .005, true, true), Outcome::Disagreement);
    EXPECT_EQ(GzkSpectral::classify(1., 1., .08, 0., true, true), Outcome::Inconclusive);
}

TEST(GzkSpectralStatistics, MissingConvergenceAndUnresolvedZerosCannotPass) {
    using GzkSpectral::Outcome;
    EXPECT_EQ(GzkSpectral::classify(0., 0., 0., 0., true, true), Outcome::Inconclusive);
    EXPECT_EQ(GzkSpectral::classify(1., 1., 0., 0., false, true), Outcome::Inconclusive);
    EXPECT_EQ(GzkSpectral::classify(1., 1., 0., 0., true, false), Outcome::Inconclusive);
    EXPECT_THROW(GzkSpectral::classify(1., 1., -.01, 0., true, true), std::invalid_argument);
}


TEST(GzkSpectralStatistics, CrossingCellIsSharedByLogEnergyOverlap) {
    const auto half = GzkSpectral::cellFractions(1e17, .02);
    EXPECT_NEAR(half[1], .5, 1e-12);
    EXPECT_NEAR(half[2], .5, 1e-12);
    const auto quarter = GzkSpectral::cellFractions(std::pow(10., 17.005), .02);
    EXPECT_NEAR(quarter[1], .25, 1e-12);
    EXPECT_NEAR(quarter[2], .75, 1e-12);
}

TEST(GzkSpectralStatistics, ReferenceReconstructionConservesCountsIncludingTails) {
    for (double energy : {1e13, 1e14, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21}) {
        const auto fractions = GzkSpectral::cellFractions(energy, .01);
        double count = 0.;
        for (double f : fractions) {
            EXPECT_GE(f, 0.);
            EXPECT_LE(f, 1.);
            count += 7.*f;
        }
        EXPECT_NEAR(count, 7., 1e-12);
    }
    EXPECT_THROW(GzkSpectral::cellFractions(1e17, 0.), std::invalid_argument);
}

TEST(GzkSpectralStatistics, GridWidthToleratesPrintedRoundingButRejectsMissingCells) {
    EXPECT_NEAR(GzkSpectral::gridLogWidth({1e17, 1.02329e17, 1.04713e17}), .01, 1e-12);
    EXPECT_THROW(GzkSpectral::gridLogWidth({1e17, 1.02329e17, 1.07152e17}), std::invalid_argument);
    EXPECT_THROW(GzkSpectral::gridLogWidth({1e17}), std::invalid_argument);
}
