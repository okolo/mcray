/* 
 * File:   Interaction.h
 * Author: ok
 *
 * Created on August 26, 2011, 4:58 PM
 */

#ifndef INTERACTION_H
#define	INTERACTION_H

#include <vector>
#include "Particle.h"
#include "Randomizer.h"
#include "Utils.h"

namespace MCPropag
{

///must be thread-safe
class Interaction : public Utils::SmartReferencedObj {
public:
    enum InteractionType
    {
        IntRandom,
        IntContEnergyLoss,
        IntDeflections
    };
    virtual ~Interaction(){};
    virtual InteractionType Type() const { return fType; };
    
    ///Mean interaction rate for random processes
    ///Energy loss rate 1/E dE/dt for continuous energy loss processes
    ///Characteristic time^-1 for deflection process
    virtual double Rate(const Particle& aParticle) const = 0;
protected:
    Interaction(InteractionType aType){fType=aType;}
private:
    InteractionType fType;
};
typedef Utils::SmartPtr<Interaction> PInteraction;

///must be thread-safe
class CELInteraction : public Interaction {
public:
	CELInteraction() : Interaction(IntContEnergyLoss) { }
	virtual void GetSecondaries(const Particle& aParticle, std::vector<Particle>& aSecondaries,
								cosmo_time aDeltaT, Randomizer& aRandomizer) const {};
	virtual CELInteraction* Clone() const = 0;
};
typedef Utils::SmartPtr<CELInteraction> PCELInteraction;

///must be thread-safe
class RandomInteraction : public Interaction {
public:
    RandomInteraction() : Interaction(IntRandom) { }
    //returns false if the process rate is zero
    virtual bool GetSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const = 0;
    virtual bool KillsPrimary() const
    {
    	return true;
    }
    virtual RandomInteraction* Clone() const = 0;
};
typedef Utils::SmartPtr<RandomInteraction> PRandomInteraction;

///must be thread-safe
class DeflectionInteraction : public Interaction {
public:
    DeflectionInteraction() : Interaction(IntDeflections) { }
    //increases aParticle.Time by deltaT and adjust particle position and momentum direction
    // Return false and leave aParticle fields unchanged if particle is subjected to the interaction
    virtual bool Propagate(cosmo_time deltaT, Particle &aParticle, Randomizer &aRandomizer) const = 0;
    virtual DeflectionInteraction* Clone() const = 0;
};
typedef Utils::SmartPtr<DeflectionInteraction> PDeflectionInteraction;

class RedShift : public CELInteraction
{
public:
	double Rate(const Particle& aParticle) const;
	CELInteraction* Clone() const { return new RedShift();}
};

class RandomInteractionS : public RandomInteraction{
public:
	bool GetSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const;
	virtual bool SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const = 0;
	virtual void SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const = 0;
};

}
#endif	/* INTERACTION_H */

