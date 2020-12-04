/*
 * Filter.h
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

#ifndef FILTER_H_
#define FILTER_H_
#include <math.h>
#include "Particle.h"
#include "Utils.h"

namespace mcray {

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

} /* namespace mcray */
#endif /* FILTER_H_ */
