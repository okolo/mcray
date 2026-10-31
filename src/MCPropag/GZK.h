//
// Created by ok on 01.08.15.
//

#ifndef MCPROPAG_GZK_H
#define MCPROPAG_GZK_H

#include "Interaction.h"
#include "Background.h"

namespace Interactions {
    using namespace MCPropag;

    class GZK : public RandomInteractionS {
    public:
        GZK(BackgroundIntegral* aBackground, int aRandSeed=0);

        virtual double Rate(const Particle &aParticle) const;

        virtual RandomInteraction *Clone() const;

        virtual bool SampleS(const Particle &aParticle, double &aS, Randomizer &aRandomizer) const;

        virtual void SampleSecondaries(Particle &aParticle, std::vector<Particle> &aSecondaries, double aS,
                                       Randomizer &aRandomizer) const;
        static void UnitTest();
        static void UnitTestEconserv(double aE, double aZ);
    private:
        Function* InitSigma(ParticleType aPrim);
        SmartPtr<BackgroundIntegral>	fBackground;
        SafePtr<Function>               fSigmaN;
        SafePtr<Function>               fSigmaP;
        double                          fMpSophia;
        double                          fMnSophia;
    };
}

#endif //MCPROPAG_GZK_H
