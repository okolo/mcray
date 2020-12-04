/*
 * PropagationEngine.h
 *
 * Author:
 *       Oleg Kalashev
 *
 * Copyright (c) 2020 Institute for Nuclear Research, RAS
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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

namespace mcray
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
    IResult                             &fResult;
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

