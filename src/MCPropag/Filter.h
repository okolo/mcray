/*
 * Filter.h
 *
 *  Created on: Oct 6, 2013
 *      Author: ok
 */

#ifndef FILTER_H_
#define FILTER_H_
#include <math.h>
#include "Particle.h"
#include "Utils.h"

namespace MCPropag {

class IFilter : public Utils::ISmartReferencedObj {
public:
	//must be thread-safe
    virtual bool Pass(const Particle& aParticles) const = 0;
};

class EnergyBasedFilter : public  Utils::TSmartReferencedObj<IFilter>
{
public:
	EnergyBasedFilter(double aEmin, double aMaxDeflectionAngle=M_PI);
	bool Pass(const Particle& aParticle) const;
    virtual ~EnergyBasedFilter(){}
private:
    double fEmin;
    double fMaxDeflectionAngle;
};



class TypeSpecificEnergyBasedFilter : public  Utils::TSmartReferencedObj<IFilter>
{
public:
    TypeSpecificEnergyBasedFilter(std::vector<int> aFilterParticleType, std::vector<double> aFilterParticleEnergy, 
                                    double aMaxDeflectionAngle=M_PI);
    bool Pass(const Particle& aParticle) const;
    virtual ~TypeSpecificEnergyBasedFilter(){}
private:
    std::vector<int> fFilterParticleType;
    std::vector<double> fFilterParticleEnergy;
    double fMaxDeflectionAngle;
};



class ParticleTypeFilter : public  Utils::TSmartReferencedObj<IFilter>
{
public:
//aParticlesToRemove is array ending with ParticleTypeEOF
    ParticleTypeFilter(const ParticleType aParticlesToDiscard[] = 0);
    void DiscardParticle(ParticleType aType) {fFilter[aType] = false;}
    bool Pass(const Particle& aParticle) const { return fFilter[aParticle.Type]; }
    virtual ~ParticleTypeFilter(){}
private:
    bool fFilter[ParticleTypeEOF];
};

class JetOutputFilter : public  Utils::TSmartReferencedObj<IFilter>
{
public:
	JetOutputFilter(double aEmin, double aJetOpenningAngle=M_PI, double aObservationAngle=M_PI_2);
	bool Pass(const Particle& aParticle) const;
    virtual ~JetOutputFilter(){}
private:
    double fEmin;
    double fJetOpenningAngle;
    double fObservationAngle;
};

class JetRuntimeFilter : public  EnergyBasedFilter
{
public:
	JetRuntimeFilter(double aEmin, double aJetOpenningAngle=M_PI, double aObservationAngle=M_PI_2):
		EnergyBasedFilter(aEmin, aJetOpenningAngle + aObservationAngle){}
};

class JetFilter : public  Utils::TSmartReferencedObj<IFilter>
{
public:
	JetFilter(double aEmin, double aJetOpenningAngle=M_PI, double aObservationAngle=M_PI_2);
	bool Pass(const Particle& aParticle) const;
    virtual ~JetFilter(){}
private:
    double fEmin;
    double fJetOpenningAngleSin;
    double fObservationAngleCos;
    double fObservationAngleSin;
};

} /* namespace MCPropag */
#endif /* FILTER_H_ */
