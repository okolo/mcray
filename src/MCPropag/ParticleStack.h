/*
 * File:   ParticleStack.h
 * Author: ok
 *
 * Created on August 26, 2011, 7:15 PM
 */

#ifndef PARTICLESTACK_H
#define	PARTICLESTACK_H

#include "Particle.h"
#include "Output.h"
#include "Utils.h"
#include "Filter.h"
#include "Observer.h"
#include <vector>
#include <set>

namespace MCPropag
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
	Observer fObserver;
};

class Result : public IResult{
public:
	Result() : fFilter(0), fIncludeFilterInPropagation(false),fAutoFlushInterval(0),fParticleCounter(0), fStopRequested(false), fObserver(-1){};
	Result(IFilter*	aOutputFilter, bool aIsPropagationFilter = false) :
		fFilter(aOutputFilter), fIncludeFilterInPropagation(aIsPropagationFilter),fAutoFlushInterval(0),fParticleCounter(0), fStopRequested(false), fObserver(-1){};
	Result(double aEmin):
		fFilter(new EnergyBasedFilter(aEmin)), fIncludeFilterInPropagation(true),fAutoFlushInterval(0),fParticleCounter(0), fStopRequested(false), fObserver(-1){}
    //set to 0 to disable auto flush
    void SetAutoFlushInterval(uint aNumberOfPrimaryParticles);
	/*!
	@param[in]  aParticle particle
	@return true if farther propagation simulation is required for the particle, false otherwise
	*/
	bool ContinuePropagation(const Particle& aParticle) const
	{
		return (GetMinTravelTimeLeft(aParticle)>0 && !fObserver.weInTarget(aParticle)) && !(fIncludeFilterInPropagation &&  fFilter && !fFilter->Pass(aParticle));
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
	Observer fObserver;
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
