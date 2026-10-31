/*
 * Deflection1D.cpp
 *
 *  Created on: Mar 22, 2012
 *      Author: ok
 */

#include "Deflection1D.h"

namespace Interactions {

Deflection1D::Deflection1D(Function* aB, Function* aLcor):
		fB(aB),fLcor(aLcor)
{
}

Deflection1D::~Deflection1D() {

}

bool Deflection1D::Propagate(cosmo_time aDeltaT, Particle& aParticle, Randomizer& aRandomizer) const
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

} /* namespace MCPropag */
