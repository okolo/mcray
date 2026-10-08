#include "PhotoDisintegration.h"
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <memory>

namespace {
using namespace mcray;
using namespace Interactions;

// Probe routing and sample the production aggregate sigma at a fixed K.
class ProbeIntegral : public BackgroundIntegral {
public:
    explicit ProbeIntegral(double k = 31. / units.Eunit) : k(k) {}
    double GetRateK(const Utils::Function& sigma, const Particle&, double) override {
        return sigma.f(k);
    }
    double GetRateS(const Utils::Function& sigma, const Particle& p, double) override {
        return sigma.f(p.Mass() * (p.Mass() + 2 * k));
    }
    double GetRateAndSampleK(const Utils::Function& sigma, const Particle& p,
                            Randomizer&, double& sampled, double error) override {
        sampled = k;
        return GetRateK(sigma, p, error);
    }
    double GetRateAndSampleS(const Utils::Function& sigma, const Particle& p,
                            Randomizer&, double& sampled, double error) override {
        sampled = p.Mass() * (p.Mass() + 2 * k);
        return GetRateS(sigma, p, error);
    }
    BackgroundIntegral* Clone() const override { return new ProbeIntegral(k); }
    double k;
};

class RoutingIntegral : public ProbeIntegral {
public:
    double GetRateK(const Utils::Function&, const Particle&, double) override { return 1.; }
    double GetRateS(const Utils::Function&, const Particle&, double) override { return 2.; }
    BackgroundIntegral* Clone() const override { return new RoutingIntegral; }
};

Particle nucleus(int a) {
    Particle p(static_cast<unsigned int>(a), 0.);
    p.Energy = 1e20 * units.eV;
    p.Weight = 2.5;
    return p;
}

double integrate(const Utils::Function& sigma) {
    // Independent composite Simpson quadrature; the production normalization
    // uses an error function instead.
    constexpr int steps = 20000;
    const double dx = (sigma.Xmax() - sigma.Xmin()) / steps;
    long double sum = sigma.f(sigma.Xmin()) + sigma.f(sigma.Xmax());
    for (int i = 1; i < steps; ++i)
        sum += (i % 2 ? 4 : 2) * sigma.f(sigma.Xmin() + i * dx);
    return static_cast<double>(sum * dx / 3.);
}

TEST(PsbStorage, DefaultAndExplicitIndexOffsetsAreDefined) {
    Utils::AutoDeletePtrArray<int> direct;
    direct.add(new int(7));
    EXPECT_EQ(*direct[0], 7);
    Utils::AutoDeletePtrArray<int> shifted;
    shifted.setIndexShift(-2);
    shifted.add(new int(9));
    EXPECT_EQ(*shifted[2], 9);
}

TEST(PsbCrossSections, GaussianAndPlateauIntegralsMatchConfiguredStrengths) {
    const TGaussCrossSection cs{10., 20., 0.7, 5., 8.};
    for (bool compatibility : {false, true}) {
        CPhotoDisintegrationCanal::TGaussianSigma sigma(cs, 3., compatibility);
        EXPECT_NEAR(integrate(sigma) * units.Eunit, 2.1, 1e-10);
        EXPECT_DOUBLE_EQ(sigma.Xmin() * units.Eunit, 10.);
        EXPECT_GE(sigma.Xmax() * units.Eunit, 30.);
    }
    CPhotoDisintegrationCanal::TLinearSigma plateau(0.7, 3.);
    EXPECT_NEAR(integrate(plateau) * units.Eunit, 2.1, 1e-10);
    EXPECT_DOUBLE_EQ(plateau.Xmin() * units.Eunit, 30.);
    EXPECT_DOUBLE_EQ(plateau.Xmax() * units.Eunit, 150.);
}

TEST(PsbNetwork, IdentityMassAndExpectedBaryonYieldAcrossTheChain) {
    CPhotoDisintegrationMap map;
    for (int a = 2; a <= 56; ++a) {
        Particle p = nucleus(a);
        ASSERT_EQ(p.MassNumber(), a);
        EXPECT_GE(p.AtomicNumber(), 0);
        EXPECT_LE(p.AtomicNumber(), a);
        EXPECT_DOUBLE_EQ(p.Mass(), a * 0.5 * (Particle::Mass(Proton) + Particle::Mass(Neutron)));
        const auto* channels = map.getOutcome(a);
        ASSERT_NE(channels, nullptr);
        for (auto* channel : *channels) {
            SCOPED_TRACE(::testing::Message() << "A=" << a << " group=" << channel->type());
            ASSERT_LT(channel->getMaxDeltaA(), a);
            double baryons = 0.;
            // Income entries independently encode the expected product yields
            // used in the deterministic transport matrix.
            for (int product = 0; product < a; ++product)
                for (auto* entry : *map.getIncome(product))
                    if (entry->iCanal == channel) {
                        EXPECT_GE(entry->iRate, 0.);
                        baryons += (product ? product : 1) * entry->iRate;
                    }
            EXPECT_NEAR(baryons, a, 1e-12);
            EXPECT_GT(integrate(*channel->Sigma()), 0.);
        }
    }
}

TEST(PsbSampling, HighEnergyChannelDoesNotIncludeGaussianBeyondItsUpperLimit) {
    CPhotoDisintegrationMap map;
    Particle p = nucleus(56);
    const double k = 31. / units.Eunit;
    const CPhotoDisintegrationCanal* multiple = nullptr;
    for (auto* channel : *map.getOutcome(56))
        if (channel->type() == EMultiple) multiple = channel;
    ASSERT_NE(multiple, nullptr);
    std::array<double, 57> expected{};
    for (int product = 0; product < 56; ++product)
        for (auto* entry : *map.getIncome(product))
            if (entry->iCanal->type() == EMultiple && entry->iCanal->isotope().A == 56)
                expected[product] += entry->iRate;
    PhotoDisintegration interaction(new ProbeIntegral(k));
    Randomizer random(20261001);
    for (int event = 0; event < 100; ++event) {
        std::vector<Particle> products;
        interaction.SampleSecondariesWithK(p, products, k, random);
        std::array<double, 57> observed{};
        double energy = 0.;
        for (const auto& product : products) {
            const int a = product.Type == Neutron ? 0 : product.MassNumber();
            observed[a] += product.Weight / p.Weight;
            EXPECT_NEAR(product.Energy / p.Energy, product.Mass() / p.Mass(), 1e-15);
            energy += product.Energy * product.Weight;
        }
        for (int a = 0; a <= 56; ++a)
            ASSERT_NEAR(observed[a], expected[a], 1e-12) << "event=" << event << " A=" << a;
        EXPECT_NEAR(energy / (p.Energy * p.Weight), 1., 1e-3);
    }
}

TEST(PsbSampling, LightNucleiHaveCorrectWeightedMultiplicities) {
    CPhotoDisintegrationMap map;
    Randomizer random(42);
    PhotoDisintegration interaction(new ProbeIntegral);
    for (int a : {3, 4, 9, 12}) {
        SCOPED_TRACE(a);
        Particle p = nucleus(a);
        const double k = (a == 3 ? 3. : 31.) / units.Eunit;
        const auto selected = a == 3 ? EDouble : EMultiple;
        std::vector<Particle> products;
        interaction.SampleSecondariesWithK(p, products, k, random);
        std::array<double, 57> expected{}, observed{};
        // Isolate A=3 double emission below the single-emission threshold;
        // isolate the high-energy plateau for the other light nuclei.
        for (int product = 0; product < a; ++product)
            for (auto* entry : *map.getIncome(product))
                if (entry->iCanal->isotope().A == a && entry->iCanal->type() == selected)
                    expected[product] += entry->iRate;
        double baryons = 0.;
        for (const auto& product : products) {
            const int output = product.Type == Neutron ? 0 : product.MassNumber();
            observed[output] += product.Weight / p.Weight;
            baryons += product.MassNumber() * product.Weight / p.Weight;
        }
        for (int i = 0; i <= 56; ++i) EXPECT_NEAR(observed[i], expected[i], 1e-12);
        EXPECT_NEAR(baryons, a, 1e-12);
    }
}

TEST(PsbSampling, KAndSProduceIdenticalWeightedProducts) {
    Particle p = nucleus(56);
    PhotoDisintegration interaction(new ProbeIntegral);
    for (double mev : {15., 31., 100.}) {
        Randomizer r1(13), r2(13);
        std::vector<Particle> byK, byS;
        const double k = mev / units.Eunit;
        interaction.SampleSecondariesWithK(p, byK, k, r1);
        interaction.SampleSecondariesWithS(p, byS, p.Mass() * (p.Mass() + 2*k), r2);
        ASSERT_EQ(byK.size(), byS.size());
        for (size_t i = 0; i < byK.size(); ++i) {
            EXPECT_EQ(byK[i].Type, byS[i].Type);
            EXPECT_DOUBLE_EQ(byK[i].Weight, byS[i].Weight);
            EXPECT_DOUBLE_EQ(byK[i].Energy, byS[i].Energy);
        }
    }
}

TEST(PsbClone, PreservesBothSamplingModes) {
    Particle p = nucleus(56);
    for (bool useK : {true, false}) {
        PhotoDisintegration original(new RoutingIntegral, useK);
        std::unique_ptr<RandomInteraction> clone(original.Clone());
        EXPECT_EQ(original.Rate(p), useK ? 1. : 2.);
        EXPECT_EQ(clone->Rate(p), original.Rate(p));
    }
}

TEST(PsbRate, ZeroSampledCrossSectionProducesNoEvent) {
    PhotoDisintegration interaction(new ProbeIntegral(0.));
    Particle p = nucleus(56);
    Randomizer random(7);
    std::vector<Particle> products;
    EXPECT_EQ(interaction.Rate(p), 0.);
    EXPECT_FALSE(interaction.GetSecondaries(p, products, random));
    EXPECT_TRUE(products.empty());
}

TEST(PsbRate, NonNuclearParticlesHaveZeroRate) {
    PhotoDisintegration interaction(new ProbeIntegral);
    for (ParticleType type : {Electron, Positron, Photon, Proton, Neutron}) {
        Particle p(type, 0.);
        p.Energy = 1e20 * units.eV;
        EXPECT_EQ(interaction.Rate(p), 0.);
    }
}
} // namespace
