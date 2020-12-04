/*
 * Interaction.cpp
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

#include "Interaction.h"
#include "Cosmology.h"

namespace mcray
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

TypeChange::TypeChange(ParticleType aPrimary, ParticleType aSecondary, double aRate):
fPrimary(aPrimary),
fSecondary(aSecondary),
fRate(aRate){};
double TypeChange::Rate(const Particle& aParticle) const{
	return aParticle.Type == fPrimary ? fRate : 0.;
}
bool TypeChange::GetSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const{
	Particle sec = aParticle;
	sec.Type = fSecondary;
	sec.fCascadeProductionTime = aParticle.Time;
    aSecondaries.push_back(sec);
    return true;
}
RandomInteraction* TypeChange::Clone() const{
	return new TypeChange(fPrimary,fSecondary,fRate);
}

}
