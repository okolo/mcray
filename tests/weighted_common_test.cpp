#include "crpropa/weighted_common.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(WeightedChannels, ExactTwoEventLikelihoodAndSplittingConservation) {
    const std::vector<double> p{.99999,.00001},guide{1.,100000.};
    const double total=p[0]*guide[0]+p[1]*guide[1];
    std::vector<double> q;for(size_t i=0;i<p.size();++i) q.push_back(.2*p[i]+.8*p[i]*guide[i]/total);
    double expected=0.;
    for(unsigned i=0;i<2;++i) for(unsigned j=0;j<2;++j) {
        const auto a=weighted::select(p,guide,.8,(i?q[0]:0.)+q[i]/2);
        const auto b=weighted::select(p,guide,.8,(j?q[0]:0.)+q[j]/2);
        expected+=q[i]*q[j]*a.likelihood*b.likelihood*(i+j);
    }
    EXPECT_NEAR(expected,2*p[1],1e-18);
    EXPECT_DOUBLE_EQ(4*weighted::reweight(.25,2.,4),.5);
    EXPECT_THROW(weighted::reweight(1e-300,1e-300),std::runtime_error);
}
TEST(WeightedChannels, OfficialNearestNodeSamplerRetainsResidualLastBranch) {
    const auto under=weighted::officialProbabilities({.2,.3,.1});
    EXPECT_DOUBLE_EQ(under[0],.2);EXPECT_DOUBLE_EQ(under[1],.3);EXPECT_DOUBLE_EQ(under[2],.5);
    const auto over=weighted::officialProbabilities({.8,.4,.1});
    EXPECT_DOUBLE_EQ(over[0],.8);EXPECT_NEAR(over[1],.2,1e-15);EXPECT_DOUBLE_EQ(over[2],0.);
    EXPECT_THROW(weighted::select({.9,.1},{1.,2.},1.,.5),std::runtime_error);
}
TEST(WeightedChannels, TopologyProposalResolvesRareBranchWithoutChangingExpectation) {
    const std::vector<double> p{.9999965,.0000035,0.},g{.1,1.,100.};
    const double q1=.2*p[1]+.8/1.1,q0=1.-q1;
    const auto a=weighted::selectDensity(p,g,.8,q0/2);
    const auto b=weighted::selectDensity(p,g,.8,q0+q1/2);
    EXPECT_EQ(a.index,0u);EXPECT_EQ(b.index,1u);
    EXPECT_NEAR(q1*b.likelihood,p[1],1e-20);
    EXPECT_NEAR(q0*a.likelihood+q1*b.likelihood,1.,1e-14);
}
