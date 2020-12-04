/*
 * Thinning.cpp
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


#include "Thinning.h"

namespace mcray
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
		if(fEnableFor[it->Type])
			totEnergy += it->Energy;
	if(totEnergy == 0.)
		return;//preserve all particles
	for(int i=aParticles.size()-1; i>=0; i--)
	{
		if(!fEnableFor[aParticles[i].Type])
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

}//end of namespace mcray
