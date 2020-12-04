/*
 * Filter.cpp
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

#include "Filter.h"

namespace mcray
{

	JetOutputFilter::JetOutputFilter(double aEmin, double aJetOpenningAngle, double aObservationAngle) :
			fEmin(aEmin),
			fJetOpenningAngle(aJetOpenningAngle),
			fObservationAngle(aObservationAngle)
	{
	}

	bool JetOutputFilter::Pass(const Particle& aParticle) const
	{
		if(aParticle.Energy<=fEmin)
			return false;
		if(fObservationAngle<M_PI_2 && aParticle.ObservationAngle() > fObservationAngle)
			return false;
		if(fJetOpenningAngle<M_PI && aParticle.JetOpenningAngle() > fJetOpenningAngle)
			return false;
		return true;
	}

	EnergyBasedFilter::EnergyBasedFilter(double aEmin, double aMaxDeflectionAngle) :
			fEmin(aEmin),
			fMaxDeflectionAngle(aMaxDeflectionAngle)
	{
	}

	bool EnergyBasedFilter::Pass(const Particle& aParticle) const
	{
		if(aParticle.Energy<=fEmin)
			return false;
		if(fMaxDeflectionAngle<M_PI && aParticle.DeflectionAngle() > fMaxDeflectionAngle)
			return false;
		return true;
	}

	JetFilter::JetFilter(double aEmin, double aJetOpenningAngle, double aObservationAngle) :
			fEmin(aEmin),
			fJetOpenningAngleSin(aJetOpenningAngle<M_PI_2 ? sin(aJetOpenningAngle) : 1.),
			fObservationAngleCos(cos(aObservationAngle)),
			fObservationAngleSin(sin(aObservationAngle))
	{
	}

	bool JetFilter::Pass(const Particle& aParticle) const
	{
		if(aParticle.Energy<=fEmin)
			return false;
		if(aParticle.DeflectionAngle()>0.5*M_PI && fJetOpenningAngleSin<1.)
			return false;//exclude particles deflected by more than Pi/2
		double sinBeta = sin(aParticle.DeflectionAngle());
		double d = cosmology.z2d(aParticle.SourceParticle->Time.z());
		double l = cosmology.z2d(aParticle.LastDeflectionTime().z());

		if(sinBeta > d/l*fJetOpenningAngleSin)
			return false;//particle can not reach observer from any angle

		double r = sqrt(d*d+l*l-2.*d*l*fObservationAngleCos);//distance(at z=0) from the point of last deflection time to the source
															//assuming that particle arrives at the edge of PSF
		if(sinBeta > d/r*fObservationAngleSin)
			return false;//particle will reach observer from outside PSF

		return true;
	}

	ParticleTypeFilter::ParticleTypeFilter(const ParticleType aParticlesToDiscard[]) {
		for(int i=0; i<ParticleTypeEOF; i++){
			fFilter[i] = true;
		}
		if(aParticlesToDiscard)
			for(int i=0; aParticlesToDiscard[i]!=ParticleTypeEOF; i++) {
				fFilter[aParticlesToDiscard[i]] = false;
			}
	}
}
