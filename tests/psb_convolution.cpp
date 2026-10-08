#include "PhotoDisintegration.h"
#include "GZK.h"
#include "frozen_inoue.h"
#include <gtest/gtest.h>
#include <cmath>

namespace {
using namespace mcray;
using namespace Interactions;

TEST(BackgroundConvolution, MonochromaticKUsesTheUpperCrossSectionBound) {
    class BoundedSigma : public Utils::Function {
    public:
        double f(double) const override { return 2.; }
        double Xmin() const override { return units.MeV; }
        double Xmax() const override { return 3*units.MeV; }
    } sigma;
    Utils::ConstFunction density(3.), energy(units.MeV);
    MonochromaticBackgroundIntegral background(density,energy,1.1,1e-6);
    Particle p(Proton,0.); p.Energy=10*p.Mass();
    // Integral K*sigma dK = sigma*(Kmax^2-Kmin^2)/2.
    const double expected=0.5*2.*(9.-1.)/2.*3./(100.*p.beta());
    EXPECT_NEAR(background.GetRateK(sigma,p)/expected,1.,1e-12);
    Randomizer random(20261008); double k=0.;
    EXPECT_NEAR(background.GetRateAndSampleK(sigma,p,random,k)/expected,1.,1e-5);
    EXPECT_GE(k,units.MeV); EXPECT_LE(k,3*units.MeV);
}

TEST(PsbConvolution, ProductionInoueExtendsToZeroInsideCompoundBackground) {
    Backgrounds::Inoue12BaselineIROSpectrum direct;
    ASSERT_DOUBLE_EQ(direct.MinZ(), 0.);
    CompoundBackground compound;
    compound.AddComponent(new Backgrounds::Inoue12BaselineIROSpectrum);
    for (double energy : {0.01, 0.1, 1.}) {
        const double e = energy * units.eV;
        const double firstPlane = direct.n(e, 0.01);
        ASSERT_GT(firstPlane, 0.);
        for (double z : {0., 0.005, 0.01}) {
            EXPECT_DOUBLE_EQ(direct.n(e, z), firstPlane);
            EXPECT_DOUBLE_EQ(compound.n(e, z), firstPlane);
        }
    }
}

TEST(PsbConvolution, FrozenInoueRatesDoNotEvolveAlongThePath) {
    if (!cosmology.IsInitialized()) cosmology.Init(0.1);
    FrozenInoueBaseline frozen;
    Backgrounds::Inoue12BaselineIROSpectrum evolving;
    for (double energy : {0.01, 0.1, 1.}) {
        const double e = energy * units.eV;
        ASSERT_GT(frozen.n(e, 0.), 0.);
        EXPECT_DOUBLE_EQ(frozen.n(e, 0.01), evolving.n(e, 0.));
        // The table's first redshift plane can cover the entire short path;
        // probe later evolution as well to detect accidental use of z.
        EXPECT_NE(frozen.n(e, 1.), evolving.n(e, 1.));
    }
    CompoundBackground photons;
    photons.AddComponent(new FrozenInoueBaseline);
    const double t = 2.7255 * units.K;
    photons.AddComponent(new PlankBackground(t, t * 1e-3, t * 1e3, 0., 0.02, true));
    for (double z : {0., 0.005, 0.01})
        EXPECT_DOUBLE_EQ(photons.n(1. * units.eV, z), frozen.n(1. * units.eV, 0.));
    PhotoDisintegration interaction(new ContinuousBackgroundIntegral(
        photons, 0.02, std::pow(10., 0.01), 0.02, 1e-4));
    Particle p(56u, 0.);
    p.Energy = 56 * 2.238721138568347e18 * units.eV;
    const double rate = interaction.Rate(p);
    ASSERT_GT(rate, 0.);
    for (double z : {0.005, 0.01}) {
        p.Time.setZ(z);
        EXPECT_NEAR(interaction.Rate(p) / rate, 1., 1e-12);
    }
}

TEST(PsbConvolution, FrozenCmbAndKSAreConsistent) {
    if (!cosmology.IsInitialized()) cosmology.Init(0.1);
    const double t = 2.7255 * units.K;
    PlankBackground cmb(t, t * 1e-3, t * 1e3, 0., 0.02, true);
    PBackgroundIntegral integral = new ContinuousBackgroundIntegral(
        cmb, 0.02, std::pow(10., 0.01), 0.02, 1e-5);
    PhotoDisintegration byK(integral, true), byS(integral, false);
    for (double gamma : {1e9, 1e10}) {
        Particle p(56u, 0.);
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
    Particle p(55u, 0.);
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

TEST(PsbConvolution, InoueDaughterRatesConvergeAndAgreeInKAndS) {
    if (!cosmology.IsInitialized()) cosmology.Init(0.1);
    const double t = 2.7255 * units.K;
    CompoundBackground photons;
    photons.AddComponent(new PlankBackground(t, t * 1e-3, t * 1e3, 0., 0.02, true));
    photons.AddComponent(new FrozenInoueBaseline);
    PBackgroundIntegral integral = new ContinuousBackgroundIntegral(
        photons, 0.02, std::pow(10., 0.01), 0.02, 1e-4);
    PhotoDisintegration byK(integral, true), byS(integral, false);
    PhotoDisintegration refined(new ContinuousBackgroundIntegral(
        photons, 0.02, std::pow(10., 0.005), 0.02, 1e-4));
    PhotoDisintegration stricter(new ContinuousBackgroundIntegral(
        photons, 0.02, std::pow(10., 0.01), 0.02, 1e-5));
    // Include the A30 daughter that previously raised GSL_EROUND and lost
    // complete Fe56 histories, as well as every other PSB daughter mass.
    for (double energy : {1e18, std::pow(10., 18.5)}) {
        for (int a = 2; a <= 56; ++a) {
            SCOPED_TRACE(::testing::Message() << "A=" << a << " E/A=" << energy);
            Particle p(static_cast<unsigned int>(a), 0.01);
            p.Energy = a * energy * units.eV;
            try {
                const double rate = byK.Rate(p);
                ASSERT_TRUE(std::isfinite(rate));
                ASSERT_GT(rate, 0.);
                EXPECT_NEAR(byS.Rate(p) / rate, 1., 5e-4);
                EXPECT_NEAR(refined.Rate(p) / rate, 1., 5e-4);
                if (a == 30 || a == 56)
                    EXPECT_NEAR(stricter.Rate(p) / rate, 1., 5e-4);
            } catch (Utils::Exception* error) {
                const std::string message = error->Message();
                delete error;
                FAIL() << message;
            }
        }
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
