/*
 * PhotoDisintegration.h
 */

#ifndef MCRAY_PD_H
#define MCRAY_PD_H

#include "Interaction.h"
#include "Background.h"

namespace Interactions {
    using namespace mcray;

    class PD : public RandomInteractionS {
    public:
        PD(BackgroundIntegral* aBackground);

        virtual double Rate(const Particle &aParticle) const;

        virtual RandomInteraction *Clone() const;

        virtual bool SampleS(const Particle &aParticle, double &aS, Randomizer &aRandomizer) const;

        virtual void SampleSecondaries(Particle &aParticle, std::vector<Particle> &aSecondaries, double aS,
                                       Randomizer &aRandomizer) const;

    private:
        Function* InitSigma(int aPrimZ, int aPrimN, double aPrimM);
        int Branching(int aPrimZ, int aPrimN, double aPrimM, double aS, Randomizer &aRandomizer) const;
        
        std::vector<SafePtr<Function> > NuclSigma;

        ParticleType                    fParticleType;
        SmartPtr<BackgroundIntegral>	fBackground;
        SafePtr<Function>               fSigma;

        std::string                     TablesDir;
        std::string                     TableSumXS;
        std::string                     TableThinXS;
        std::string                     TableEps;

        int                             MaxNBranches = 50;
    };
}

#endif //MCRAY_GZK_H