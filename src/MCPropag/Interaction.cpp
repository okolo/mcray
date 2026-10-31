/* 
 * File:   Interaction.cpp
 * Author: ok
 * 
 * Created on August 26, 2011, 4:58 PM
 */

#include "Interaction.h"
#include "Cosmology.h"

namespace MCPropag
{
double RedShift::Rate(const Particle& aParticle) const
{
	return cosmology.eLossRate(aParticle.Time.z());
}

bool RandomInteractionS::GetSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const
{
	double s = 0.;
	if(!SampleS(aParticle, s, aRandomizer))
		return false;
	debug.PrintLine(ToString(s));
	SampleSecondaries(aParticle, aSecondaries, s, aRandomizer);
	return true;
}

}
