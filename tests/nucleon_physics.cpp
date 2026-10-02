#include "GZK.h"
#include "PPP.h"
#include "ProtonPP.h"
#include <gtest/gtest.h>
#include <memory>
#include <gsl/gsl_integration.h>

namespace {
using namespace mcray;
using namespace Interactions;

class FixedRestFramePhoton : public BackgroundIntegral {
public:
    explicit FixedRestFramePhoton(double k) : k(k) {}
    double GetRateS(const Utils::Function& sigma, const Particle& p, double) override {
        return sigma.f(p.Mass() * (p.Mass() + 2 * k));
    }
    double GetRateK(const Utils::Function& sigma, const Particle&, double) override {
        return sigma.f(k);
    }
    double GetRateAndSampleS(const Utils::Function& sigma, const Particle& p,
                            Randomizer&, double& s, double error) override {
        s = p.Mass() * (p.Mass() + 2 * k);
        return GetRateS(sigma, p, error);
    }
    double GetRateAndSampleK(const Utils::Function& sigma, const Particle& p,
                            Randomizer&, double& sampled, double error) override {
        sampled = k;
        return GetRateK(sigma, p, error);
    }
    BackgroundIntegral* Clone() const override { return new FixedRestFramePhoton(k); }
private:
    double k;
};

TEST(NucleonGzk, SupportedRatesAndClonePreserveCrossSections) {
    GZK interaction(new FixedRestFramePhoton(0.5 * units.GeV));
    std::unique_ptr<RandomInteraction> clone(interaction.Clone());
    for (auto type : {Proton, Neutron}) {
        Particle p(type, 0.);
        p.Energy = 1e20 * units.eV;
        ASSERT_GT(interaction.Rate(p), 0.);
        EXPECT_DOUBLE_EQ(clone->Rate(p), interaction.Rate(p));
    }
    EXPECT_DOUBLE_EQ(interaction.Rate(Particle(A56, 0.)), 0.);
}

TEST(NucleonGzk, SophiaEventsPreserveNucleonNumberWeightAndEnergy) {
    GZK interaction(new FixedRestFramePhoton(0.5 * units.GeV), 20261071);
    Randomizer random(20261071);
    for (auto type : {Proton, Neutron}) {
        Particle p(type, 0.);
        p.Energy = 1e20 * units.eV;
        p.Weight = 2.5;
        for (int event = 0; event < 20; ++event) {
            std::vector<Particle> products;
            ASSERT_TRUE(interaction.GetSecondaries(p, products, random));
            int nucleons = 0;
            double energy = 0.;
            for (const auto& product : products) {
                EXPECT_DOUBLE_EQ(product.Weight, p.Weight);
                EXPECT_GT(product.Energy, 0.);
                energy += product.Energy;
                if (product.Type == Proton || product.Type == Neutron) ++nucleons;
            }
            EXPECT_EQ(nucleons, 1);
            EXPECT_NEAR(energy / p.Energy, 1., 1e-6);
        }
    }
}

TEST(ProtonPairs, ScalingAndClonePreserveWeightedProductionRate) {
    PPP ordinary(new FixedRestFramePhoton(10. * units.MeV));
    PPP scaled(new FixedRestFramePhoton(10. * units.MeV), 3.);
    std::unique_ptr<RandomInteraction> clone(scaled.Clone());
    Particle p(Proton, 0.);
    p.Energy = 1e19 * units.eV;
    p.Weight = 2.5;
    ASSERT_GT(ordinary.Rate(p), 0.);
    EXPECT_DOUBLE_EQ(clone->Rate(p), scaled.Rate(p));
    EXPECT_NEAR(3. * scaled.Rate(p) / ordinary.Rate(p), 1., 1e-14);
    EXPECT_FALSE(scaled.KillsPrimary());
    Randomizer random(20261081);
    std::vector<Particle> products;
    try {
        ASSERT_TRUE(scaled.GetSecondaries(p, products, random));
    } catch (Utils::Exception* error) {
        const std::string message = error->Message();
        delete error;
        FAIL() << message;
    }
    ASSERT_EQ(products.size(), 1u);
    EXPECT_TRUE(products[0].Type == Electron || products[0].Type == Positron);
    EXPECT_DOUBLE_EQ(products[0].Weight, 6. * p.Weight);
    EXPECT_GT(products[0].Energy, 0.);
    EXPECT_LT(products[0].Energy, p.Energy);
    EXPECT_DOUBLE_EQ(p.Energy, 1e19 * units.eV);
}

TEST(ProtonPairs, CelCloneAndUnsupportedSpeciesAreExplicit) {
    ProtonPPcel cel(new FixedRestFramePhoton(10. * units.MeV));
    std::unique_ptr<CELInteraction> clone(cel.Clone());
    Particle p(Proton, 0.);
    p.Energy = 1e19 * units.eV;
    ASSERT_GT(cel.Rate(p), 0.);
    EXPECT_DOUBLE_EQ(clone->Rate(p), cel.Rate(p));
    PPP random(new FixedRestFramePhoton(10. * units.MeV));
    for (auto type : {Neutron, A56}) {
        Particle unsupported(type, 0.);
        EXPECT_DOUBLE_EQ(cel.Rate(unsupported), 0.);
        EXPECT_DOUBLE_EQ(random.Rate(unsupported), 0.);
    }
}

TEST(ProtonPairs, NearThresholdSamplingRemainsFinite) {
    const double m = Particle::Mass(Proton)/Particle::Mass(Electron);
    const double samplingThreshold = 2.*m/(m-1.);
    for (double k : {2., 2.0005, samplingThreshold-1e-8,
                     samplingThreshold, samplingThreshold+1e-6, 2.01}) {
        SCOPED_TRACE(k);
        PPP interaction(new FixedRestFramePhoton(k*Particle::Mass(Electron)));
        Particle p(Proton, 0.);
        const double energy = 1e19*units.eV;
        p.Energy = energy;
        const double rate = interaction.Rate(p);
        ASSERT_TRUE(std::isfinite(rate));
        ASSERT_GE(rate, 0.);
        if (rate == 0.) continue;
        Randomizer random(20266701);
        std::vector<Particle> products;
        try {
            const bool generated = interaction.GetSecondaries(p, products, random);
            EXPECT_DOUBLE_EQ(p.Energy, energy);
            EXPECT_DOUBLE_EQ(p.Weight, 1.);
            if (!generated) {
                EXPECT_TRUE(products.empty());
                continue;
            }
            ASSERT_EQ(products.size(), 1u);
            EXPECT_TRUE(std::isfinite(products[0].Energy));
            EXPECT_GT(products[0].Energy, 0.);
            EXPECT_LE(products[0].Energy, p.Energy);
        } catch (Utils::Exception* error) {
            const std::string message = error->Message();
            delete error;
            FAIL() << message;
        }
    }
}

class ProtonPppSampling : public ::testing::TestWithParam<double> {};

TEST_P(ProtonPppSampling, MatchesIndependentRestFrameQuadrature) {
    const double k = GetParam(); // Rest-frame photon energy in electron masses.
    const double massRatio = Particle::Mass(Proton)/Particle::Mass(Electron);
    // Integrate independently over rest-frame electron energy and angle.
    // The sampler integrates a transformed one-dimensional lab-energy kernel.
    auto integrate = [&](int nodes) {
        std::unique_ptr<gsl_integration_glfixed_table, decltype(&gsl_integration_glfixed_table_free)>
            table(gsl_integration_glfixed_table_alloc(nodes), gsl_integration_glfixed_table_free);
        std::array<long double, 3> sums{};
        for (int i = 0; i < nodes; ++i) {
            double energy, energyWeight;
            gsl_integration_glfixed_point(1., k-1., i, &energy, &energyWeight, table.get());
            for (int j = 0; j < nodes; ++j) {
                double angle, angleWeight;
                gsl_integration_glfixed_point(-1., 1., j, &angle, &angleWeight, table.get());
                const double density = PPP::DiffSigma(k, energy, angle)*energyWeight*angleWeight;
                const double r = (energy-angle*std::sqrt(energy*energy-1.))/massRatio;
                sums[0] += density;
                sums[1] += density*r;
                sums[2] += density*r*r;
            }
        }
        return std::array<double, 2>{static_cast<double>(sums[1]/sums[0]), static_cast<double>(sums[2]/sums[0])};
    };
    const auto coarse = integrate(256), expected = integrate(512);
    for (int moment = 0; moment < 2; ++moment) {
        std::cout << "k=" << k << " moment=" << moment+1 << " quadrature_relative_change="
                  << std::abs(coarse[moment]/expected[moment]-1.) << '\n';
        ASSERT_NEAR(coarse[moment]/expected[moment], 1., 1e-3);
    }
    PPP interaction(new FixedRestFramePhoton(k*Particle::Mass(Electron)));
    Particle p(Proton, 0.);
    p.Energy = 1e19 * units.eV;
    Randomizer random(20261601 + static_cast<unsigned int>(100*k));
    constexpr int samples = 10000;
    std::array<double, 2> sums{}, squares{};
    try {
        for (int i = 0; i < samples; ++i) {
            std::vector<Particle> products;
            ASSERT_TRUE(interaction.GetSecondaries(p, products, random));
            ASSERT_EQ(products.size(), 1u);
            const double r = products[0].Energy/p.Energy;
            for (int moment = 0; moment < 2; ++moment) {
                const double value = moment == 0 ? r : r*r;
                sums[moment] += value;
                squares[moment] += value*value;
            }
        }
    } catch (Utils::Exception* error) {
        const std::string message = error->Message();
        delete error;
        FAIL() << message;
    }
    for (int moment = 0; moment < 2; ++moment) {
        const double mean = sums[moment]/samples;
        const double se = std::sqrt((squares[moment]/samples-mean*mean)/(samples-1));
        std::cout << "k=" << k << " moment=" << moment+1 << " quadrature="
                  << expected[moment] << " MC=" << mean << " SE=" << se << '\n';
        // The kernel neglects recoil; the sampler retains recoil in its bounds.
        EXPECT_NEAR(mean, expected[moment], 4*se+.01*expected[moment]);
    }
}
INSTANTIATE_TEST_SUITE_P(RestFrameEnergies, ProtonPppSampling,
                        ::testing::Values(2.2, 3., 5., 10., 20., 50.));
} // namespace
