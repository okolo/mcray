/*
 * Deflection1D.cpp
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

#include "Deflection1D.h"

namespace Interactions {

Deflection1D::Deflection1D(Function* aB, Function* aLcor):
		fB(aB),fLcor(aLcor)
{
}

Deflection1D::~Deflection1D() {

}

bool Deflection1D::Propagate(cosmo_time aDeltaT, Particle& aParticle, Randomizer& aRandomizer, ILogger* aLogger) const
{
	double Lcor = fLcor->f(aParticle.Time.z());

	if(aParticle.ElectricCharge()!=0)
	{
		if(aParticle.CorrelatedBpath>=Lcor)
		{
			aParticle.Deflection2 += (aParticle.CorrelatedDeflection*aParticle.CorrelatedDeflection);
			aParticle.LastB = -1.;
		}
		if(aParticle.LastB < 0.)
		{
			aParticle.CorrelatedBpath = 0;
			aParticle.CorrelatedDeflection = 0;
		}
		cosmo_time corDeflDt = Lcor - aParticle.CorrelatedBpath;

		for(cosmo_time timeLeft = aDeltaT;  timeLeft>0;)
		{
			cosmo_time rate = Rate(aParticle);
			if(corDeflDt>timeLeft)
			{
				aParticle.CorrelatedBpath += timeLeft;
				aParticle.CorrelatedDeflection += rate*timeLeft;//Theta in radians
				break;
			}
			else
			{
				aParticle.CorrelatedBpath = 0.;
				double corDefl = aParticle.CorrelatedDeflection + rate*corDeflDt;//Theta in radians
				aParticle.Deflection2 += corDefl*corDefl;
				aParticle.CorrelatedDeflection = 0.;
				timeLeft -= corDeflDt;
				corDeflDt = Lcor;
				aParticle.LastB = -1.;
			}
		}
		aParticle.PropagateFreely(aDeltaT);//we don't modify Pdir in this class
		return true;
	}
	else
	{
		aParticle.CorrelatedBpath += aDeltaT;
		return false;
	}
}

//deflection rate in radians per unit length (internal units) assuming constant field
double Deflection1D::Rate(const Particle& aParticle) const
{
	int q = abs(aParticle.ElectricCharge());
	if(q)
	{
		if(aParticle.LastB<0.)
			aParticle.LastB = fabs(fB->f(aParticle.Time.z()));
		double deflectionRate = ((double)q)*0.52*M_PI/180./(aParticle.Energy/units.TeV)/(10.*units.kpc)*aParticle.LastB*1e15;//formula (14) arxiv/1106.5508v2 or (9) from astro-ph/9604098
		ASSERT_VALID_NO(deflectionRate);
		return deflectionRate;
	}
	return 0.;
}

Function* Deflection1D::RandomDirectionB::Clone() const
{
	return new RandomDirectionB(fAbsBvalue, fRandomizer.CreateIndependent());
}


double Deflection1D::RandomDirectionB::f(double aZ) const
{
	double cos = fRandomizer.Rand()*2.0 - 1;
	double sin = sqrt(1-cos*cos);
	return fAbsBvalue*sin;
}

double Deflection1D::ConstComovingCorrelationLength::f(double aZ) const
{
	return fCorL/(1.+aZ);
}

Function* Deflection1D::ConstComovingCorrelationLength::Clone() const
{
	return new RandomDirectionB(fCorL);
}

} /* namespace mcray */
