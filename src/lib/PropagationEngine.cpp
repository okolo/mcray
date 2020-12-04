/*
 * PropagationEngine.cpp
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

#include <math.h>
#include <numeric>
#include "PropagationEngine.h"
#include <omp.h>

using namespace std;

namespace mcray
{

unsigned long int PropagationEngine::fLastRandSeed = 0;

PropagationEngine::PropagationEngine(ParticleStack &aStack, IResult &aResult, unsigned long int aRandSeed) :
		fThinning(0),
		fFilter(0),
        fStack(aStack),
        fResult(aResult),
        fRandomizer(aRandSeed ? aRandSeed:fLastRandSeed),
        fAccuracy(0.01),
        fMinStep(0.),
        fInfoPrinted(false),
        fLogger(0)
{
	if(!aRandSeed)
		fLastRandSeed++;

}

PropagationEngine::PropagationEngine(PropagationEngine &aEngine, unsigned long int aRandSeed) :
		fThinning(aEngine.fThinning),
		fFilter(aEngine.fFilter),
		fStack(aEngine.fStack),
		fResult(aEngine.fResult),
		fRandomizer(aRandSeed ? aRandSeed:fLastRandSeed),
		fAccuracy(aEngine.fAccuracy),
        fMinStep(aEngine.fMinStep),
        fLogger(aEngine.fLogger)
{
	if(!aRandSeed)
		fLastRandSeed++;
    for(vector<PCELInteraction>::iterator it = aEngine.fContinuousEnergyLosses.begin(); it != aEngine.fContinuousEnergyLosses.end(); it++)
    {
    	fContinuousEnergyLosses.push_back((*it)->Clone());
    }
    for(vector<PDeflectionInteraction>::iterator it = aEngine.fDeflectionInteractions.begin(); it != aEngine.fDeflectionInteractions.end(); it++)
    {
    	fDeflectionInteractions.push_back((*it)->Clone());
    }
    for(vector<PRandomInteraction>::iterator it = aEngine.fRandomInteractions.begin(); it != aEngine.fRandomInteractions.end(); it++)
    {
    	fRandomInteractions.push_back((*it)->Clone());
    }
    for(vector<double>::iterator it = aEngine.fRandomInteractionRateMultipliers.begin(); it != aEngine.fRandomInteractionRateMultipliers.end(); it++)
    {
        fRandomInteractionRateMultipliers.push_back(*it);
    }
}

void PropagationEngine::Run(double aMaxRateVariation) {
    try {
        Particle p(Electron, 0.);//will be reset by fStack.Pop
        while (fStack.Pop(p)) {
            try {
                while (true) {
                    if(fLogger)
                        fLogger->log(p); // record initial state
                    if (aMaxRateVariation > 0)
                        RunParticleWithVariableRate(p, aMaxRateVariation);
                    else
                        RunParticle(p);
                    if(fResult.AbortRequested())
                        return;
                    if (!fSecondaryParticles.size())
                        break;
                    p = fSecondaryParticles.back();
                    fSecondaryParticles.pop_back();
                }

            }catch(Exception* ex){
                fFinalParticles.clear();//don't save current paricle cascade
                std::string message = ex->Message() + "\tprimary " + p.ToString();
                LOG_ERROR(message);
                std::cerr << message << std::endl;
            }
            if(fFinalParticles.size()){
                fResult.Add(fFinalParticles);
                fFinalParticles.clear();
            }
        }
    }catch(Exception* ex)
    {
        LOG_ERROR(ex->Message());
        std::cerr << ex->Message() << std::endl;
    }
}

void PropagationEngine::RunMultithreadReleaseOnly(double aMaxRateVariation, int aNumberOfThreads)
{
#ifdef _DEBUG
	Run(aMaxRateVariation);
#else
	RunMultithread(aMaxRateVariation, aNumberOfThreads);
#endif
}

void PropagationEngine::RunMultithread(double aMaxRateVariation, int aNumberOfThreads)
{
    if(aNumberOfThreads>1)
        omp_set_num_threads(aNumberOfThreads);
    else if(aNumberOfThreads==1){
        Run(aMaxRateVariation);
        return;
    }
#pragma omp parallel
{
  SafePtr<PropagationEngine> pe;
#pragma omp critical (PropagationEngine)
  {
	   pe = new PropagationEngine(*this, fRandomizer.CreateIndependent());
  }

#pragma omp barrier

#pragma omp master
{
    if(!fInfoPrinted){
        int nthreads = omp_get_num_threads();
        cout << "Running in " << nthreads << " threads" << '\n';
        fInfoPrinted = true;
    }
}

  pe->Run(aMaxRateVariation);
}
}

void PropagationEngine::move(Particle &aParticle, cosmo_time aDt, double aElossRate)
{
    cosmo_time eLossHalf = exp(-0.5 * aElossRate * aDt);//relative energy loss during dt/2
    aParticle.Energy = aParticle.Energy * eLossHalf;//mean energy during period [t, t+dt]
    bool deflected = false;
    for(vector<PDeflectionInteraction>::iterator it = fDeflectionInteractions.begin(); it != fDeflectionInteractions.end(); it++)
    {
        if((*it)->Propagate(aDt, aParticle, fRandomizer))//should not modify aParticle.Energy and aParticle.Time
        {
            if(deflected)
                Exception::Throw("Particle treated by more than one deflection interaction");
            deflected = true;
        }
    }
    if(!deflected)
        aParticle.PropagateFreely(aDt);
    aParticle.Energy = aParticle.Energy * eLossHalf;
    if(fLogger)
        fLogger->log(aParticle);
}

double PropagationEngine::GetInteractionRate(const Particle &aParticle) const {
    double totRate = 0.;
    for(vector<PRandomInteraction>::const_iterator it = fRandomInteractions.begin(); it != fRandomInteractions.end(); it++)
    {
        double rate = (*it)->Rate(aParticle);
        totRate += rate;
    }
    return totRate;
}

void PropagationEngine::RunParticle(Particle& aParticle) 
{
    while((!fResult.AbortRequested()) && fResult.ContinuePropagation(aParticle))
    {
    	if(fFilter && !fFilter->Pass(aParticle))
    		return;
        double eLossRate = 0.;
        //double bRate = 0.;
        double randRate = 0.;
        cosmo_time maxDt = fResult.GetMinTravelTimeLeft(aParticle);
        ASSERT(maxDt>0);
        for(vector<PCELInteraction>::iterator it = fContinuousEnergyLosses.begin(); it != fContinuousEnergyLosses.end(); it++)
        {
            eLossRate += (*it)->Rate(aParticle);
        }
        //for(vector<PDeflectionInteraction>::iterator it = fDeflectionInteractions.begin(); it != fDeflectionInteractions.end(); it++)
        //{
        //    bRate += (*it)->Rate(aParticle);
        //}
        vector<double> rates;
        for(vector<PRandomInteraction>::iterator it = fRandomInteractions.begin(); it != fRandomInteractions.end(); it++)
        {
            double rate = (*it)->Rate(aParticle);
            randRate += rate;
            rates.push_back(rate);
        }
        cosmo_time dt = eLossRate>0 ? fAccuracy/eLossRate : maxDt;//fAccuracy<1
        if(dt>maxDt)
            dt = maxDt;

        double rand = -1;
        bool hasInteraction = false;
        if(randRate>0)
        {
			double freePath = fRandomizer.GetFreePath(1./randRate, rand);
			if(freePath<=dt)
			{
				dt = freePath;
				hasInteraction = true;
			}
        }

        move(aParticle, dt, eLossRate);


        for(vector<PCELInteraction>::iterator it = fContinuousEnergyLosses.begin(); it != fContinuousEnergyLosses.end(); it++)
        {
        	vector<Particle> secondaries;
            (*it)->GetSecondaries(aParticle,secondaries,dt,fRandomizer);
            HandleSecondaries(aParticle, secondaries, 1.);
        }

        if(hasInteraction)
        {
            //select interaction
            double r = fRandomizer.Rand() * randRate;
            double totRate = 0.;
            RandomInteraction* interaction = 0;
            for(unsigned int i=0; i<rates.size(); i++)
            {
                totRate += rates[i];
                if(totRate>=r)
                {
                    interaction = fRandomInteractions[i];
                    break;
                }
            }
            uint maxError_count = 3;
            bool completed = false;
            vector<Particle> secondaries;
            bool interactionOccurred = false;
            for(uint i=0; !completed; i++)
            {
                try {

                    interactionOccurred = interaction->GetSecondaries(aParticle, secondaries,
                                                    fRandomizer);//may return false if particle energy is close to process threshold

                    completed = true;
                }catch(Exception* ex)
                {
                    if(i==maxError_count) {
                        std::cerr << "recovery failed after " << maxError_count-1 << "attempts" << std::endl;
                        throw ex;
                    }
                    std::cerr << ex->Message() << "\n\tattempting to recover" << std::endl;
                    secondaries.clear();
                }
            }
            if(interactionOccurred) {
                HandleSecondaries(aParticle, secondaries, 1.);//todo: implement calculation speedup for rare processes using rate multiplier in RunParticle (see implementation in RunParticleWithVariableRate)
                if (interaction->KillsPrimary())
                    return;
            }
        }
    }
    fFinalParticles.push_back(aParticle);
}

void PropagationEngine::RunParticleWithVariableRate(Particle& aParticle, double aMaxRelRateError){
    const int maxNumberOfStepAdjustments = 50;//2e-50 ~ 1e-15 which is double type accuracy
    LOG_VERBOSE("RPWVR:0\t" + aParticle.ToString());
    while((!fResult.AbortRequested()) && fResult.ContinuePropagation(aParticle))
    {
        LOG_VERBOSE("RPWVR:1\t" + aParticle.ToString());
        if(fFilter && !fFilter->Pass(aParticle))
            return;

        double eLossRate1 = 0.;//total continuous energy loss rate at start point

        cosmo_time maxDt = fResult.GetMinTravelTimeLeft(aParticle);//maximal step
        ASSERT(maxDt>0);

        vector<double> celRates1(fContinuousEnergyLosses.size(),0.);// continuous energy loss rates for separate processes at start point
        int iRate = 0;
        for(vector<PCELInteraction>::iterator it = fContinuousEnergyLosses.begin(); it != fContinuousEnergyLosses.end(); it++)
        {
            double rate = (*it)->Rate(aParticle);
            celRates1[iRate++] = rate;
            eLossRate1 += rate;
        }

        cosmo_time dtCEL = eLossRate1 >0 ? fAccuracy/ eLossRate1 : maxDt;//fAccuracy<1
        if(dtCEL >maxDt)
            dtCEL = maxDt;
        vector<double> celRates2(fContinuousEnergyLosses.size(), 0.);
        //find optimal step for CEL interaction
        Particle part2 = aParticle;
        LOG_VERBOSE("RPWVR:2\t" + aParticle.ToString() + "\tdtCel = " + ToString(dtCEL/units.Mpc));
        for(int nSteps=0; true; nSteps++) {
            move(part2, dtCEL, eLossRate1);
            LOG_VERBOSE("RPWVR:2.1\t" + part2.ToString() + "\tdtCel = " + ToString(dtCEL/units.Mpc));
            double eLossRate2 = 0.;
            iRate = 0;
            for (vector<PCELInteraction>::iterator it = fContinuousEnergyLosses.begin();
                 it != fContinuousEnergyLosses.end(); it++) {
                celRates2[iRate] = (*it)->Rate(part2);
                eLossRate2 += celRates2[iRate];
                LOG_VERBOSE("RPWVR:2.2\t" + part2.ToString() + "\tRelDiffCel[" + ToString(iRate) + "] = " + ToString(MathUtils::RelDifference(celRates2[iRate],celRates1[iRate])));
                iRate++;
            }
            double relDif = MathUtils::RelDifference(eLossRate1,eLossRate2);
            LOG_VERBOSE("RPWVR:2.2\t" + part2.ToString() + "\tRelDiffCel = " + ToString(relDif));
            if(relDif<aMaxRelRateError)
                break;
            if(dtCEL<=fMinStep){
                LOG_WARNING("CEL rate variance " + ToString(relDif) + " (limitted by step size)")
                break;
            }
            if(nSteps==maxNumberOfStepAdjustments){
                Exception::Throw("Unable to achieve desired CEL rate variance. Minimal value " + ToString(relDif) + "\nSet minimal step and/or decrease desired accuracy");
            }
            else// adjust step and repeat
            {
                dtCEL *= (0.5 * aMaxRelRateError / relDif);
                part2 = aParticle;
            }
        };
        LOG_VERBOSE("RPWVR:3\t" + aParticle.ToString() + "\tdtCel = " + ToString(dtCEL/units.Mpc));
        double randRate1 = 0.;
        vector<double> randRates1(fRandomInteractions.size());
        iRate = 0;
        for(vector<PRandomInteraction>::iterator it = fRandomInteractions.begin(); it != fRandomInteractions.end(); it++)
        {
            double rate = (*it)->Rate(aParticle)*fRandomInteractionRateMultipliers[iRate];
            randRate1 += rate;
            randRates1[iRate++] = rate;
        }

        vector<double> randRates2(fRandomInteractions.size(),0);
        //double randRate2 = 0.;
        double dtRand = dtCEL;
        //vector<double> maxValidStep(fRandomInteractions.size(), dtCEL);
        iRate = 0;
        double relDif = -1.;
        for(int nSteps=0; true; nSteps++) {
            LOG_VERBOSE("RPWVR:4\t" + part2.ToString() + "\tdtRand/dtCEL = " + ToString(dtRand/dtCEL));
            for(; iRate < fRandomInteractions.size(); iRate++)
            {
                randRates2[iRate] = fRandomInteractions[iRate]->Rate(part2)*fRandomInteractionRateMultipliers[iRate];
                relDif = MathUtils::RelDifference(randRates1[iRate],randRates2[iRate]);
                if(relDif>aMaxRelRateError){
                    LOG_VERBOSE("RPWVR:4.1\tRelDiffCel[" + ToString(iRate) + "] = " + ToString(relDif));
                    double dtMult = 0.5*aMaxRelRateError/relDif;
                    dtRand *= dtMult;
                    for(int iPrevIntr=0; iPrevIntr <iRate; iPrevIntr++){// (using linear interpolation to avoid recalculation of well known rates)
                        randRates2[iPrevIntr] = randRates1[iPrevIntr] + (randRates2[iPrevIntr]-randRates1[iPrevIntr])*dtMult;
                    }
                    break;
                }
            }
            if(iRate==fRandomInteractions.size())//all tests passed
                break;

            if(dtRand<=fMinStep){
                LOG_WARNING("Rate variance " + ToString(relDif) + " (limitted by step size)")
                break;
            }
            if(nSteps == maxNumberOfStepAdjustments){
                Exception::Throw("Unable to achieve desired rate variance. Minimal value " + ToString(relDif) + "\nSet minimal step and/or decrease desired accuracy");
            }
            //else, adjust step and repeat
            part2 = aParticle;
            move(part2, dtRand, eLossRate1);
        };
        LOG_VERBOSE("RPWVR:5\t" + aParticle.ToString() + "\tdtRand = " + ToString(dtRand/units.Mpc));
        double randRate2 = std::accumulate(randRates2.begin(), randRates2.end(), 0.);

        double rand = -1;
        RandomInteraction* interaction = 0;
        double sec_weight = 1.;
        double sqrt_noninteracting_weight = 1.;
        double meanRandRate = 0.5*(randRate1+randRate2);
        double dt = dtRand;
        if(meanRandRate > 0)
        {
            double freePath = fRandomizer.GetFreePath(1./ meanRandRate, rand);
            if(freePath <= dtRand)
            {
                dt = freePath;
                //calculate mean rate and save it to randRates1 vector
                double sumRates = 0.;
                for(iRate = fRandomInteractions.size()-1; iRate>=0; iRate--){
                    randRates1[iRate] += (dt/dtRand*(randRates2[iRate]-randRates1[iRate]));
                    sumRates += randRates1[iRate];
                }
                //select interaction
                double r = fRandomizer.Rand() * sumRates;
                double totRate = 0.;
                for(size_t i=0; i< randRates1.size(); i++)
                {
                    totRate += randRates1[i];
                    if(totRate>=r)
                    {
                        interaction = fRandomInteractions[i];
                        sec_weight = 1./fRandomInteractionRateMultipliers[i];
                        break;
                    }
                }
            }
            else{
                double ln_noninteracting_weight = 0.;
                for(iRate = fRandomInteractions.size()-1; iRate>=0; iRate--){
                    if(fRandomInteractionRateMultipliers[iRate]!=1.){
                        ln_noninteracting_weight += (0.5*(randRates2[iRate]+randRates1[iRate])*(1.-1./fRandomInteractionRateMultipliers[iRate]));
                    }
                }
                sqrt_noninteracting_weight = exp(0.5*ln_noninteracting_weight);
            }

        }
        LOG_VERBOSE("RPWVR:6\t" + aParticle.ToString()+ "\tdt = " + ToString(dt/units.Mpc));
        aParticle.Weight *= sqrt_noninteracting_weight;
        move(aParticle, 0.5*dt, eLossRate1);
        for(vector<PCELInteraction>::iterator it = fContinuousEnergyLosses.begin(); it != fContinuousEnergyLosses.end(); it++)
        {
            vector<Particle> secondaries;
            (*it)->GetSecondaries(aParticle,secondaries, dt,fRandomizer);
            HandleSecondaries(aParticle, secondaries, 1.);//currently continuous energy loss weighting is not implemented
        }
        move(aParticle, 0.5*dt, eLossRate1);
        aParticle.Weight *= sqrt_noninteracting_weight;
        LOG_VERBOSE("RPWVR:6.1\t" + aParticle.ToString());

        if(interaction)
        {
            uint maxError_count = 3;
            bool completed = false;
            vector<Particle> secondaries;
            bool interactionOccurred = false;
            for(uint i=0; !completed; i++)
            {
                LOG_VERBOSE("RPWVR:7\t" + aParticle.ToString());
                try {

                    interactionOccurred = interaction->GetSecondaries(aParticle, secondaries,
                                                                      fRandomizer);//may return false if particle energy is close to process threshold

                    completed = true;
                }catch(Exception* ex)
                {
                    if(i==maxError_count) {
                        LOG_ERROR("recovery failed after " + ToString(maxError_count-1) + " attempts");
                        throw ex;
                    }
                    LOG_WARNING(ex->Message() + "\tattempting to recover");
                    secondaries.clear();
                }
            }
            if(interactionOccurred) {
                LOG_VERBOSE("RPWVR:8\t" + aParticle.ToString());
                HandleSecondaries(aParticle, secondaries, sec_weight);
                if (interaction->KillsPrimary()){
                    LOG_VERBOSE("RPWVR:-2\t" + aParticle.ToString());
                    return;
                }
            }
        }
        LOG_VERBOSE("RPWVR:9\t" + aParticle.ToString());
    }
    LOG_VERBOSE("RPWVR:10\t" + aParticle.ToString());
    fFinalParticles.push_back(aParticle);
    LOG_VERBOSE("RPWVR:-1\t" + aParticle.ToString());
}

void PropagationEngine::HandleSecondaries(Particle& aPrimary, vector<Particle>& aSecondaries, double aWeightMultiplier)
{
	long numberOfInteractions = aPrimary.Ninteractions + 1;
	if(fFilter)
		for(int i=aSecondaries.size()-1; i>=0; i--) {
            aSecondaries[i].Weight *= aWeightMultiplier;
            if (!fFilter->Pass(aSecondaries[i])) {
                LOG_VERBOSE("HS:filter\t" + (aSecondaries.begin() + i)->ToString());
                aSecondaries.erase(aSecondaries.begin() + i);
            }
		}
	if(fThinning)
		fThinning->Run(aSecondaries, fRandomizer);
	for(std::vector<Particle>::iterator  it = aSecondaries.begin(); it != aSecondaries.end(); it++)
	{
        it->SetParent(aPrimary);
        it->Ninteractions = numberOfInteractions;
        fSecondaryParticles.push_back(*it);
	}
}

}
