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
#include <cmath>

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

PhotonThinning::PhotonThinning(double alpha, double threshold)
    : alpha_(alpha), threshold_(threshold)
{
    if (!std::isfinite(alpha) || alpha < 0 || alpha > 10)
        Utils::Exception::ThrowFatal("Photon thinning alpha must be finite and in [0,10]");
    if (!std::isfinite(threshold) || threshold <= 0)
        Utils::Exception::ThrowFatal("Photon thinning threshold must be finite and positive");
}

void PhotonThinning::Run(std::vector<Particle>& particles, Randomizer& randomizer)
{
    if (alpha_ == 0)
        return; // Disabled means no validation, modification or RNG draws.

    // Stable in-place compaction is linear even for a large secondary list.
    size_t retained = 0;
    for (size_t i = 0; i < particles.size(); ++i) {
        Particle& particle = particles[i];
        if (particle.Type == Photon) {
            if (!std::isfinite(particle.Energy) || particle.Energy <= 0 ||
                !std::isfinite(particle.Weight) || particle.Weight < 0)
                Utils::Exception::ThrowFatal("Photon thinning requires positive finite photon energy and finite nonnegative weight");
            if (particle.Energy < threshold_) {
                // Avoid the fork's product of powers, which can overflow even
                // when their ratio is moderate. Logs also handle E/Et underflow.
                const double ratio = particle.Energy / threshold_;
                const double logRatio = ratio > 0 ? std::log(ratio)
                    : std::log(particle.Energy) - std::log(threshold_);
                const double probability = std::exp(alpha_ * logRatio);
                if (!(probability > 0) || !std::isfinite(probability))
                    Utils::Exception::ThrowFatal("Photon thinning survival probability is not representable");
                const double weight = particle.Weight / probability;
                // Validate BEFORE sampling: dropping only overflowing survivors
                // would bias the estimator. The whole run must be failed.
                if (!std::isfinite(weight))
                    Utils::Exception::ThrowFatal("Photon thinning corrected weight is not representable");
                if (randomizer.Rand() > probability)
                    continue;
                particle.Weight = weight;
            }
        }
        if (retained != i)
            particles[retained] = particle;
        ++retained;
    }
    particles.erase(particles.begin() + retained, particles.end());
}

}//end of namespace mcray
