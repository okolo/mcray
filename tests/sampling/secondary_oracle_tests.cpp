#include "secondary_oracle.h"
#include <gtest/gtest.h>
using namespace secondary;
TEST(SecondaryOracle, PairChargeAndMultiplicityAreCorrelated) {
    Quadrature q(32);
    std::vector<Region> obs{{"all","inclusive","bin",0,1},
        {"e","electron","bin",0,.5L},{"p","positron","bin",0,.5L},
        {"low","low","bin",0,.5L},{"forbidden","low","tail",0,1e-5L}};
    const auto mass=conditional(true,100,0,obs,q);
    EXPECT_NEAR(double(mass[0][2]),1,1e-14);
    EXPECT_NEAR(double(mass[1][1]),.5,1e-14);
    EXPECT_NEAR(double(mass[2][1]),.5,1e-14);
    EXPECT_NEAR(double(mass[3][1]),1,1e-14);
    EXPECT_EQ(mass[4][1],0);
}
TEST(SecondaryOracle, ThomsonLimitHasKnownMoments) {
    // In the Thomson limit photon energy / its maximum has density
    // proportional to 1 - 2z + 2z^2. Its mean is exactly 1/2.
    Quadrature q(32);const Real s=1.000001L,a=1/s,b=1-a;
    const auto norm=q.logarithmic([&](Real y){return comptonDensity(y,s);},a,1);
    const auto moment=q.logarithmic([&](Real y){return (1-y)/b*comptonDensity(y,s);},a,1);
    EXPECT_NEAR(double(moment/norm),.5,2e-7);
}
TEST(SecondaryOracle, ExtremeComptonTailIsIntegratedDirectly) {
    Quadrature q(32);
    std::vector<Region> obs{{"tail","recoil","tail",1e-12L,1e-10L},
        {"all","recoil","bin",0,1}};
    const auto mass=conditional(false,1e14L,1e-12L,obs,q);
    EXPECT_GT(mass[0][1],0);
    EXPECT_NEAR(double(mass[1][1]),1,1e-14);
    EXPECT_NEAR(double(mass[0][0]+mass[0][1]),1,1e-14);
}
