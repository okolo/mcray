/*
 * Interaction.h
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

#ifndef INTERACTION_H
#define	INTERACTION_H

#include <vector>
#include "Particle.h"
#include "Randomizer.h"
#include "Utils.h"

namespace mcray
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

//Auxiliary class for debugging etc.
class TypeChange : public RandomInteraction{
	ParticleType fPrimary;
	ParticleType fSecondary;
	double fRate;
public:
	TypeChange(ParticleType aPrimary, ParticleType aSecondary, double aRate);
	virtual double Rate(const Particle& aParticle) const;
	bool GetSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const;
	virtual RandomInteraction* Clone() const;
};

}
#endif	/* INTERACTION_H */

