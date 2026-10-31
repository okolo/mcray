/* 
 * File:   Thinning.cpp
 * Author: ok
 * 
 * Created on August 26, 2011, 7:01 PM
 */

#include "Thinning.h"
#include "Utils.h"

namespace MCPropag
{

EnergyBasedThinning::EnergyBasedThinning(double aAlpha, FilterMode aFilterMode) :
		fAlpha(aAlpha),
		fEnableFor(ParticleTypeEOF, aFilterMode != BlackList)
{
	ASSERT(aAlpha>=0 && aAlpha<=1);
}

void EnergyBasedThinning::Run(std::vector<Particle>& aParticles, Randomizer& aRandomizer)
{
	if(aParticles.size()<=1)
		return;
	double totEnergy = 0.;
	for(std::vector<Particle>::const_iterator it = aParticles.begin(); it != aParticles.end(); it++)
		if (fEnableFor[it->Type])
			totEnergy += it->Energy;
	if(totEnergy == 0.)
		return;//preserve all particles
	for(int i=aParticles.size()-1; i>=0; i--)
	{
		if (!fEnableFor[aParticles[i].Type])
			continue;
		double probability = pow(aParticles[i].Energy/totEnergy, fAlpha);//probability to survive
		if(aRandomizer.Rand()>probability) {
			LOG_VERBOSE("thinning:erase\t" + (aParticles.begin() + i)->ToString());
			aParticles.erase(aParticles.begin() + i);
		}
		else
			aParticles[i].Weight /= probability;
	}
}

PhotonThinning::PhotonThinning(double aAlpha, double aEnergy) :
		fAlpha(aAlpha),
		fEnergy(aEnergy)
{
	ASSERT(aAlpha>=0 && aAlpha<=10);
}
void PhotonThinning::Run(std::vector<Particle>& aParticles, Randomizer& aRandomizer)
{
	for(int i=aParticles.size()-1; i>=0; i--)
	{
		if (fAlpha == 0)
			continue;

		double probability = 1;
		if ((aParticles[i].Type == 2) and  (aParticles[i].Energy < fEnergy))
		{
			probability = pow(fEnergy, -fAlpha)*pow(aParticles[i].Energy, fAlpha);
			if(aRandomizer.Rand()>probability) {
				LOG_VERBOSE("thinning:erase\t" + (aParticles.begin() + i)->ToString());
				aParticles.erase(aParticles.begin() + i);
			}
			else
				aParticles[i].Weight /= probability;
		}

	}
}


}//end of namespace MCPropag
