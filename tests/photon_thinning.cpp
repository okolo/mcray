#include "PropagationEngine.h"
#include "Thinning.h"
#include <cmath>
#include <functional>
#include <gtest/gtest.h>
#include <limits>

namespace {
using namespace mcray;

void Fatal(const std::function<void()>& action, const std::string& text) {
    try {
        action();
        FAIL() << "Expected fatal error: " << text;
    } catch (Utils::Exception* error) {
        EXPECT_FALSE(error->Retryable());
        EXPECT_NE(error->Message().find(text), std::string::npos) << error->Message();
        delete error;
    }
}
Particle Make(ParticleType type, double energy, double weight = 3) {
    Particle particle(type, 0);
    particle.Energy = energy;
    particle.Weight = weight;
    return particle;
}

TEST(PhotonThinning, ConstructorRejectsInvalidSettingsInRelease) {
    for (double alpha : {-1., 10.01, std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::quiet_NaN()})
        Fatal([&] { PhotonThinning policy(alpha, 10); }, "alpha");
    for (double threshold : {0., -1., std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()})
        Fatal([&] { PhotonThinning policy(1, threshold); }, "threshold");
}
TEST(PhotonThinning, DisabledAndInapplicableParticlesDoNotConsumeRandomNumbers) {
    for (double alpha : {0., 1., 10.}) {
        PhotonThinning policy(alpha, 10);
        Randomizer random(9201), reference(9201);
        std::vector<Particle> particles{Make(Electron, 1), Make(Positron, 2),
                                       Make(Proton, 3), Make(Photon, 10), Make(Photon, 20)};
        if (alpha == 0)
            particles.push_back(Make(Photon, 1));
        const auto original = particles;
        policy.Run(particles, random);
        ASSERT_EQ(particles.size(), original.size());
        for (size_t i = 0; i < particles.size(); ++i) {
            EXPECT_EQ(particles[i].id, original[i].id);
            EXPECT_DOUBLE_EQ(particles[i].Weight, original[i].Weight);
        }
        EXPECT_DOUBLE_EQ(random.Rand(), reference.Rand());
    }
}
TEST(PhotonThinning, SinglePhotonUsesBernoulliSamplingAndInverseProbabilityWeight) {
    PhotonThinning policy(2, 4); // E=2 gives p=1/4, also for a one-element list.
    Randomizer random(527), reference(527);
    size_t survived = 0, removed = 0;
    for (int i = 0; i < 100; ++i) {
        std::vector<Particle> particles{Make(Photon, 2)};
        const bool keep = reference.Rand() <= .25;
        policy.Run(particles, random);
        ASSERT_EQ(particles.size(), keep ? 1u : 0u);
        if (keep) {
            EXPECT_DOUBLE_EQ(particles[0].Weight, 12);
            ++survived;
        } else ++removed;
    }
    EXPECT_GT(survived, 0u);
    EXPECT_GT(removed, 0u);
    EXPECT_DOUBLE_EQ(random.Rand(), reference.Rand());
}
TEST(PhotonThinning, CompactionPreservesOrderAndParticleMetadata) {
    Particle primary(Proton, 0);
    std::vector<Particle> particles;
    for (int i = 0; i < 50; ++i) {
        auto particle = Make(i % 2 ? Photon : Electron, 1);
        particle.SourceParticle = &primary;
        particle.Ninteractions = 7;
        particle.fPrevId = primary.id;
        particle.X[0] = i;
        particle.Pdir[0] = .6;
        particle.Pdir[2] = .8;
        particle.interactionData[0] = &primary;
        particles.push_back(particle);
    }
    Particle nucleus(56u, 0);
    nucleus.Energy = 1;
    particles.push_back(nucleus);
    const auto original = particles;
    Randomizer random(103), reference(103);
    PhotonThinning policy(1, 2);
    policy.Run(particles, random);
    size_t kept = 0;
    for (const auto& before : original) {
        const bool photon = before.Type == Photon;
        if (photon && reference.Rand() > .5) continue;
        ASSERT_LT(kept, particles.size());
        const auto& after = particles[kept++];
        EXPECT_EQ(after.id, before.id);
        EXPECT_EQ(after.Type, before.Type);
        EXPECT_EQ(after.SourceParticle, before.SourceParticle);
        EXPECT_EQ(after.fPrevId, before.fPrevId);
        EXPECT_EQ(after.Ninteractions, before.Ninteractions);
        EXPECT_EQ(after.interactionData[0], before.interactionData[0]);
        EXPECT_DOUBLE_EQ(after.Energy, before.Energy);
        EXPECT_DOUBLE_EQ(after.Weight, before.Weight * (photon ? 2 : 1));
        for (size_t j = 0; j < 3; ++j) {
            EXPECT_EQ(after.X[j], before.X[j]);
            EXPECT_EQ(after.Pdir[j], before.Pdir[j]);
        }
        if (after.Type == Nucleus) {
            EXPECT_EQ(after.MassNumber(), 56u);
            EXPECT_EQ(after.AtomicNumber(), 26u);
            EXPECT_EQ(after.Mass(), before.Mass());
        }
    }
    EXPECT_EQ(kept, particles.size());
    EXPECT_DOUBLE_EQ(random.Rand(), reference.Rand());
}
TEST(PhotonThinning, InvalidPhotonDataAndNumericalRangeFailBeforeDrawing) {
    for (const auto& particle : {Make(Photon, 0), Make(Photon, -1),
             Make(Photon, std::numeric_limits<double>::infinity()),
             Make(Photon, std::numeric_limits<double>::quiet_NaN()),
             Make(Photon, 1, -1), Make(Photon, 1, std::numeric_limits<double>::infinity()),
             Make(Photon, 1, std::numeric_limits<double>::quiet_NaN()),
             Make(Photon, 1e-300), Make(Photon, 1, std::numeric_limits<double>::max())}) {
        Randomizer random(45), reference(45);
        PhotonThinning policy(10, 2);
        std::vector<Particle> particles{particle};
        Fatal([&] { policy.Run(particles, random); }, "Photon thinning");
        EXPECT_DOUBLE_EQ(random.Rand(), reference.Rand());
    }
}
TEST(PhotonThinning, ExtremeUnitsDoNotOverflowSeparatePowersOrUnderflowRatio) {
    for (double threshold : {4e-250, 4., 4e250}) {
        Randomizer random(527), reference(527);
        PhotonThinning policy(2, threshold);
        for (int i = 0; i < 100; ++i) {
            std::vector<Particle> particles{Make(Photon, threshold / 2)};
            const bool keep = reference.Rand() <= .25;
            policy.Run(particles, random);
            ASSERT_EQ(particles.size(), keep ? 1u : 0u);
            if (keep) EXPECT_NEAR(particles[0].Weight, 12, 1e-12);
        }
    }
    PhotonThinning policy(.1, 1e300);
    Randomizer random(6), reference(6);
    std::vector<Particle> particles{Make(Photon, 1e-300)};
    policy.Run(particles, random); // Ratio underflows, but p=1e-60 is valid.
    reference.Rand();
    EXPECT_TRUE(particles.empty());
    EXPECT_DOUBLE_EQ(random.Rand(), reference.Rand());
}
TEST(PhotonThinning, FixedIndependentHistoriesMatchAnalyticFirstAndSecondMoments) {
    // Frozen software check: 100000 one-photon histories at each probability,
    // seed 73102, six analytic standard errors for each raw moment. No reruns
    // or tolerance adjustment. These are independent Bernoulli histories.
    constexpr int histories = 100000;
    constexpr double weight = 3;
    for (double probability : {.1, .5, .9}) {
        SCOPED_TRACE(probability);
        PhotonThinning policy(1, 1);
        Randomizer random(73102);
        long double sum = 0, sumSquares = 0;
        for (int i = 0; i < histories; ++i) {
            std::vector<Particle> particles{Make(Photon, probability, weight)};
            policy.Run(particles, random);
            const long double value = particles.empty() ? 0 : particles[0].Weight;
            sum += value;
            sumSquares += value * value;
        }
        const double mean = sum / histories, second = sumSquares / histories;
        const double variance = weight * weight * (1 / probability - 1);
        const double secondVariance = std::pow(weight, 4) * (1 - probability) /
                                      std::pow(probability, 3);
        EXPECT_NEAR(mean, weight, 6 * std::sqrt(variance / histories));
        EXPECT_NEAR(second, weight * weight / probability,
                    6 * std::sqrt(secondVariance / histories));
    }
}

class OnePhoton : public RandomInteraction {
public:
    explicit OnePhoton(int& calls) : calls_(calls) {}
    double Rate(const Particle& particle) const override {
        return particle.Type == Electron ? 1 / units.Mpc : 0;
    }
    bool GetSecondaries(Particle& primary, std::vector<Particle>& secondaries,
                        Randomizer&) const override {
        ++calls_;
        Particle secondary = primary;
        secondary.SetType(Photon);
        secondary.Energy = 1e-300;
        secondaries.push_back(secondary);
        return true;
    }
    RandomInteraction* Clone() const override { return new OnePhoton(calls_); }
private:
    int& calls_;
};
TEST(PhotonThinning, EngineCountsUnrepresentableHistoryAsFailureWithoutRetry) {
    if (!cosmology.IsInitialized()) cosmology.Init(1);
    for (double variation : {0., .1}) {
        ParticleStack stack;
        Particle primary = Make(Electron, 1e12 * units.eV);
        primary.Time.setZ(.01);
        stack.AddPrimary(primary);
        Result result;
        PropagationEngine engine(stack, result, 17);
        PhotonThinning policy(10, 1);
        engine.SetThinning(&policy);
        int calls = 0;
        engine.AddInteraction(new OnePhoton(calls));
        engine.Run(variation);
        EXPECT_EQ(engine.FailureCount(), 1u);
        EXPECT_EQ(calls, 1);
    }
}
} // namespace
