/*
 * ElmagTest.h
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

#ifndef ELMAGTEST_H_
#define ELMAGTEST_H_

#include "Particle.h"
#include "Interaction.h"
#include "ICS.h"
#include "GammaPP.h"

namespace Test {
using namespace mcray;

/*!
 * Propagation parameters such as background model, minimal photon energy etc. are defined in Elmag code,
 * see user_sp102.f90 module user_variables)
 */
class ElmagProxy {
public:
	static double Rate(const Particle& aParticle);
	static double RateICScel(const Particle& aParticle);
	static bool GetSecondaries(const Particle& aParticle, std::vector<Particle>& aSecondaries);
	static bool SampleS(const Particle& aParticle, double& aS);
	static bool GetSecondaries(const Particle& aParticle, std::vector<Particle>& aSecondaries, double aS);
	static void Enable();
	static inline bool Enabled() { return fEnabled; }
private:
	static int ParticleTypeToInt(ParticleType aType);
	static bool fEnabled;
};

enum ElmagUsageMode
{
	ElmagDoNotUse = 0,
	ElmagCalculateRate = 1,
	ElmagSampleS = 2,
	ElmagSampleSecondaryEnergy = 4,
	ElmagCalculateCelRate = 8
};

class ElmagIcsPP : public RandomInteraction{
protected:
	ElmagIcsPP():
		fExternal(0),fElmagUsageMode(ElmagCalculateRate|ElmagSampleS|ElmagSampleSecondaryEnergy){}
public:
	//combination of ElmagUsageMode flags
	ElmagIcsPP(int aElmagUsageMode, RandomInteractionS* pExternalInteraction):
		fExternal(pExternalInteraction),fElmagUsageMode(aElmagUsageMode)
	{
		if(fElmagUsageMode!=ElmagDoNotUse)
			ElmagProxy::Enable();
	}
	virtual bool Filter(ParticleType aType) const = 0;
	double Rate(const Particle& aParticle) const
	{
		if(!Filter(aParticle.Type))
			return 0.;
		if(ElmagCalculateRate)
			return ElmagProxy::Rate(aParticle);
		else
			return fExternal->Rate(aParticle);
	}
	bool GetSecondaries(const Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const
	{
		bool result = false;
		if(!Filter(aParticle.Type))
			return result;
		double s = 0.;
		result = (fElmagUsageMode&ElmagSampleS)?ElmagProxy::SampleS(aParticle, s) : fExternal->SampleS(aParticle, s, aRandomizer);
		if(!result)
			return false;
		if(fElmagUsageMode&ElmagSampleSecondaryEnergy)
			return ElmagProxy::GetSecondaries(aParticle, aSecondaries,s);
		fExternal->SampleSecondaries(aParticle, aSecondaries, s, aRandomizer);
		return true;
	}
protected:
	SmartPtr<RandomInteractionS> fExternal;
	int							 fElmagUsageMode;
};

class ElmagIcs : public ElmagIcsPP{
public:
	ElmagIcs(){}
	ElmagIcs(int aElmagUsageMode, Interactions::ICS* pExternalInteraction):
		ElmagIcsPP(aElmagUsageMode, pExternalInteraction){};

	RandomInteraction* Clone() const { return new ElmagIcs(fElmagUsageMode, (Interactions::ICS*)fExternal->Clone()); }
	bool Filter(ParticleType aType) const
	{
		return aType == Electron || aType == Positron;
	}
};

class ElmagPP : public ElmagIcsPP{
public:
	ElmagPP(){}
	ElmagPP(int aElmagUsageMode, Interactions::GammaPP* pExternalInteraction):
		ElmagIcsPP(aElmagUsageMode, pExternalInteraction){};
	RandomInteraction* Clone() const { return new ElmagPP(fElmagUsageMode,(Interactions::GammaPP*)fExternal->Clone()); }
	virtual bool Filter(ParticleType aType) const
	{
		return aType == Photon;
	}
};

/*!
 * Continuous energy loss rate due to Inverse Compton Scattering (ICS) with soft secondary photons
 * (energy lower than minimal value defined in Elmag code see user_sp102.f90 module user_variables) is calculated in this class
 */
class ElmagIcsCEL : public CELInteraction
{
public:
	ElmagIcsCEL()
	{
		ElmagProxy::Enable();
	}
	double Rate(const Particle& aParticle) const
	{
		return ElmagProxy::RateICScel(aParticle);
	}
	CELInteraction* Clone() const { return new ElmagIcsCEL();}
};

class ElmagTests
{
public:
	static int Main(int argc, char** argv);
private:
	static std::string HumanReadableMode(int aMode);
	static void KachelriesTest(double aZmax, int aElmagUsageModePP, int aElmagUsageModeICS);
	static void RatesTest();
};

} //namespace Test

#endif /* ELMAGTEST_H_ */
