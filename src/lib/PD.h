#ifndef MCRAY_PD_H
#define MCRAY_PD_H

#include "Interaction.h"
#include "Background.h"
#include <memory>

namespace Interactions {
using namespace mcray;

class PD : public RandomInteractionS {
public:
    explicit PD(BackgroundIntegral* aBackground);
    double Rate(const Particle& aParticle) const;
    RandomInteraction* Clone() const;
    bool SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const;
    void SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries,
                           double aS, Randomizer& aRandomizer) const;

private:
    struct Channel {
        int code;
        std::shared_ptr<const Function> sigma;
    };
    struct Isotope {
        std::shared_ptr<const Function> total;
        std::vector<Channel> channels;
    };
    PD(BackgroundIntegral* aBackground, const std::vector<Isotope>& aIsotopes);
    const Isotope& Data(ParticleType aType) const;
    int Branching(const Isotope& aData, double aS, Randomizer& aRandomizer) const;
    SmartPtr<BackgroundIntegral> fBackground;
    std::vector<Isotope> fIsotopes;
};
}
#endif
