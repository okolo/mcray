#include "PhotoDisintegration.h"
#include "GZK.h"
#include <gtest/gtest.h>
#include <cmath>

namespace {
using namespace mcray;
using namespace Interactions;

TEST(PsbConvolution, FrozenCmbAndKSAreConsistent) {
    if (!cosmology.IsInitialized()) cosmology.Init(0.1);
    const double t = 2.7255 * units.K;
    PlankBackground cmb(t, t * 1e-3, t * 1e3, 0., 0.02, true);
    PBackgroundIntegral integral = new ContinuousBackgroundIntegral(
        cmb, 0.02, std::pow(10., 0.01), 0.02, 1e-5);
    PhotoDisintegration byK(integral, true), byS(integral, false);
    for (double gamma : {1e9, 1e10}) {
        Particle p(A56, 0.);
        p.Energy = p.Mass() * gamma;
        const double rate = byK.Rate(p);
        ASSERT_GT(rate, 0.);
        EXPECT_NEAR(byS.Rate(p) / rate, 1., 1e-4);
        p.Time.setZ(0.01);
        EXPECT_NEAR(byK.Rate(p) / rate, 1., 1e-12);
        // Planck I(q)=-(T/pi^2) log(1-exp(-q/T)); beta=1 at these gamma.
        // Integrate independently in K, channel by channel.
        CPhotoDisintegrationMap map;
        double expected = 0.;
        for (auto* channel : *map.getOutcome(56)) {
            const auto* sigma = channel->Sigma();
            constexpr int steps = 20000;
            const double dx = (sigma->Xmax() - sigma->Xmin()) / steps;
            long double sum = 0.;
            for (int i = 0; i <= steps; ++i) {
                const double k = sigma->Xmin() + i * dx;
                const double kernel = -std::log(-std::expm1(-k / (2 * gamma * t)));
                sum += (i == 0 || i == steps ? 1 : i % 2 ? 4 : 2) * k * sigma->f(k) * kernel;
            }
            expected += static_cast<double>(sum * dx / 3.) * t / (2 * M_PI * M_PI * gamma * gamma);
        }
        EXPECT_NEAR(rate / expected, 1., 1e-3);
    }
}

TEST(PsbConvolution, DaughterSamplingCompletesWithoutDiscardingHistories) {
    if (!cosmology.IsInitialized()) cosmology.Init(0.1);
    const double t = 2.7255 * units.K;
    PlankBackground cmb(t, t * 1e-3, t * 1e3, 0., 0.02, true);
    PhotoDisintegration interaction(new ContinuousBackgroundIntegral(
        cmb, 0.02, std::pow(10., 0.01), 0.02, 1e-4));
    Particle p(A55, 0.);
    p.Energy = 55 * 2.238721138568347e18 * units.eV;
    p.Time.setZ(0.0044801);
    Randomizer random(1);
    std::vector<Particle> products;
    try {
        ASSERT_TRUE(interaction.GetSecondaries(p, products, random));
    } catch (Utils::Exception* error) {
        const std::string message = error->Message();
        delete error;
        FAIL() << message;
    }
}

TEST(NucleonConvolution, GzkRateConvergesOnFrozenCmb) {
    if (!cosmology.IsInitialized()) cosmology.Init(0.1);
    const double t = 2.7255 * units.K;
    PlankBackground cmb(t, t * 1e-3, t * 1e3, 0., 0.02, true);
    Interactions::GZK precise(new ContinuousBackgroundIntegral(
        cmb, 0.02, std::pow(10., 0.01), 0.02, 1e-4));
    Interactions::GZK coarse(new ContinuousBackgroundIntegral(
        cmb, 0.02, std::pow(10., 0.01), 0.02, 1e-3));
    Particle p(Proton, 0.);
    p.Energy = 1e20 * units.eV;
    try {
        const double rate = precise.Rate(p);
        ASSERT_GT(rate, 0.);
        EXPECT_NEAR(coarse.Rate(p) / rate, 1., 2e-3);
    } catch (Utils::Exception* error) {
        const std::string message = error->Message();
        delete error;
        FAIL() << message;
    }
}

} // namespace
