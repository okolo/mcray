/*
 * ParticleStack.h
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


#ifndef PARTICLESTACK_H
#define	PARTICLESTACK_H

#include "Particle.h"
#include "Output.h"
#include "Utils.h"
#include "Filter.h"
#include <vector>
#include <set>

namespace mcray
{

//Thread-safe stack
class ParticleStack {
public:
    ParticleStack(int aDebugOutputPeriod=-1);
	//TODO: remove this method and RawOutput::CopyToStack() when intermediate output is implemented via multiple IResult instances in PropagEngine
	//Currently used in OmegaCascade application
    void AddSecondaries(const std::vector<Particle> &aParticles);
    //void AddSecondary(const Particle &aParticle);
	void AddPrimary(const Particle &aParticle);
    bool Pop(Particle& aParticle);
    virtual ~ParticleStack();
    inline void SetOutputPrefix(std::string aPrefix) { fOutputPrefix = aPrefix; }
private:
    /*
    struct Less : public std::binary_function<Particle, Particle, bool>
        {
          bool
          operator()(const Particle& __x, const Particle& __y) const
          {
        	  return __x.Energy < __y.Energy;
          }
        };
    std::multiset<Particle,Less> fParticles;*/
    std::vector<Particle> fParticles;
	std::vector<Particle*> fPrimaryParticles;//stored here to guarantee validity of Particle::SourceParticle field
    time_t fLastDebugOutputTime;
    int fDebugOutputPeriod;
    std::string fOutputPrefix;
};

class IResult {
public:
	virtual bool AbortRequested() const = 0;
	virtual bool ContinuePropagation(const Particle& aParticle) const = 0;
	virtual void Add(std::vector<Particle> aParticles) = 0;
	virtual cosmo_time GetMinTravelTimeLeft(const Particle& aParticles) const = 0;
};

class Result : public IResult{
public:
	Result() : fFilter(0), fIncludeFilterInPropagation(false),fAutoFlushInterval(0),fParticleCounter(0), fStopRequested(false){};
	Result(IFilter*	aOutputFilter, bool aIsPropagationFilter = false) :
		fFilter(aOutputFilter), fIncludeFilterInPropagation(aIsPropagationFilter),fAutoFlushInterval(0),fParticleCounter(0), fStopRequested(false){};
	Result(double aEmin):
		fFilter(new EnergyBasedFilter(aEmin)), fIncludeFilterInPropagation(true),fAutoFlushInterval(0),fParticleCounter(0), fStopRequested(false){}
    //set to 0 to disable auto flush
    void SetAutoFlushInterval(uint aNumberOfPrimaryParticles);
	/*!
	@param[in]  aParticle particle
	@return true if farther propagation simulation is required for the particle, false otherwise
	*/
	bool ContinuePropagation(const Particle& aParticle) const
	{
		return GetMinTravelTimeLeft(aParticle)>0 && !(fIncludeFilterInPropagation &&  fFilter && !fFilter->Pass(aParticle));
	}
	inline const CosmoTime& T() const {return fT;}
	cosmo_time GetMinTravelTimeLeft(const Particle& aParticles) const;
	~Result();
	void SetEndTime(const CosmoTime& aT);
    void Add(std::vector<Particle> aParticles);
    void Add(const Particle& aParticle);
    void Flush(bool aFinal=true);
    void AddOutput(IOutput* aOutput) { fOutputs.push_back(aOutput); }
    inline void SetOutputFilter(IFilter* aFilter) {fFilter = aFilter;}
	void EnableFlushOnCtrlC();
	bool AbortRequested() const { return fStopRequested; }
private:
	static void SIGINT_handler(int);
	static std::vector<Result*> fCtrlCFlushList;
    CosmoTime			fT;//default end time is z=0 (see default constructor of CosmoTime)
    SmartPtr<IFilter>	fFilter;
    std::vector<SmartPtr<IOutput> >  fOutputs;
    bool fIncludeFilterInPropagation;
    uint fAutoFlushInterval;
    uint fParticleCounter;
	bool fStopRequested;
};

}
#endif	/* PARTICLESTACK_H */

