/*
 * File:   PropagationEngine.h
 * Author: ok
 *
 * Created on August 26, 2011, 6:37 PM
 */

#ifndef PROPAGATIONENGINE_H
#define	PROPAGATIONENGINE_H

#include <vector>
#include "Interaction.h"
#include "Thinning.h"
#include "Filter.h"
#include "ParticleStack.h"
#include "Randomizer.h"
#include "Utils.h"
#include "Logger.h"

namespace MCPropag
{

class PropagationEngine {
public:
    PropagationEngine(ParticleStack &aStack, IResult &aResult, unsigned long int aRandSeed=0);
    PropagationEngine(PropagationEngine &aCloneFrom, unsigned long int aRandSeed=0);
    void Run(double aMaxRateVariation=0);
    void RunMultithread(double aMaxRateVariation=0, int aNumberOfThreads = 0);
    void RunMultithreadReleaseOnly(double aMaxRateVariation=0, int aNumberOfThreads = 0);
    inline void SetThinning(IThinning* aThinning) {fThinning = aThinning;}

    /// Set propagation filter (particles not passing the filter will be erased)
    inline void SetPropagationFilter(IFilter* aFilter) {fFilter = aFilter;}
    inline void AddInteraction(PCELInteraction aInt) {fContinuousEnergyLosses.push_back(aInt);}
    inline void AddInteraction(PDeflectionInteraction aInt) {fDeflectionInteractions.push_back(aInt);}
    inline void AddInteraction(PRandomInteraction aInt, double aRateMultiplier=1.) {//todo: implement calculation speedup for rare processes using rate multiplier in RunParticle (see implementation in RunParticleWithVariableRate)
        fRandomInteractions.push_back(aInt);
        fRandomInteractionRateMultipliers.push_back(aRateMultiplier);
    }
    //set minimal step for variable rate mode
    inline void SetMinimalStep(double aMinStep) { fMinStep = aMinStep; }
    inline void SetLogger(ILogger* aLogger){fLogger=aLogger;}
    double GetInteractionRate(const Particle &aParticle) const;
    IResult                             &fResult;
private:
    void RunParticle(Particle& aParticle);
    void move(Particle &aParticle, cosmo_time aDt, double aElossRate);

    //special mode for highly variable rates
    void RunParticleWithVariableRate(Particle& aParticle, double aMaxRelRateError);
    void HandleSecondaries(Particle& aPrimary, std::vector<Particle>& aSecondaries, double aWeightMultiplier);
    std::vector<PRandomInteraction>     fRandomInteractions;
    std::vector<PDeflectionInteraction> fDeflectionInteractions;
    std::vector<PCELInteraction>        fContinuousEnergyLosses;
    std::vector<double>                 fRandomInteractionRateMultipliers;//useful for very rare interactions
    IThinning                           *fThinning;
    SmartPtr<IFilter>					fFilter;
    ParticleStack                       &fStack;
    Randomizer                          fRandomizer;
    double                              fAccuracy;
    double                              fMinStep;
    std::vector<Particle>               fSecondaryParticles;
    std::vector<Particle>               fFinalParticles;
    static unsigned long int			fLastRandSeed;
    bool                                fInfoPrinted;
    ILogger*                            fLogger;
};

}
#endif	/* PROPAGATIONENGINE_H */
