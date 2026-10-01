#include "GZK.h"
#include "PPP.h"
#include "ProtonPP.h"
#include <gtest/gtest.h>
#include <memory>

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
} // namespace
