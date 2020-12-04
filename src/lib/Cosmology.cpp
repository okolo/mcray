/*
 * Cosmology.cpp
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

#include "Cosmology.h"
#include "Units.h"
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_gamma.h"
#include <iostream>
#include "Particle.h"

namespace mcray {

Cosmology cosmology;

void Cosmology::Init(cosmo_time aMaxZ, cosmo_time aOmegaVac, cosmo_time aHubbleKm_sec_Mpc, int aAccuracy)
{
	ASSERT(t0==0);//init should be called only once
	ASSERT(aMaxZ>0);
	ASSERT(aOmegaVac>=0 && aOmegaVac<=1.);
	ASSERT(aAccuracy>100);
	Accuracy = aAccuracy;

	H = 1e5*aHubbleKm_sec_Mpc/Units::Mpc_in_cm*units.Tunit;
	fMaxZ = aMaxZ;

	Lv = aOmegaVac;
	Lm = 1.-Lv;
	sqrtLv = sqrt(Lv);
	sqrtLm = sqrt(Lm);

	if(Lv <= 0.)//Lyambda = 0
	{
		t0 = 2./3./H;
	}
	else
	{
		t0 = log((2. + 2.*sqrtLv - Lm)/Lm)/(3.*sqrtLv)/H;
		//t0 = 2./3./H/sqrtLv*log((1+sqrtLv)/sqrtLm);
	}
	initZ2T();
	initZ2D();
}

Cosmology::Cosmology()
{
	Lv = 0.;
	Lm = 0.;
	t0 = 0.;
	sqrtLm = 0.;
	sqrtLv = 0.;
	sqrtLv = 0.;
	Accuracy = 0;
}

void Cosmology::initZ2D()
{
	ASSERT(fMaxZ>0.);// must be initialized

	z.push_back(0);
	d.push_back(0);

	cosmo_time aMin = 1./(fMaxZ+1);
	cosmo_time da = (1.-aMin)/Accuracy;
	cosmo_time da2 = 0.5*da;
	cosmo_time a = 1.-da;
	cosmo_time prev = diffD(1.);
	cosmo_time sum = 0.;

	for(int i=1; i<=Accuracy; i++, a -= da)
	{
		sum += prev;
		sum += 4.*diffD(a+da2);
		prev = diffD(a);
		sum += prev;
		z.push_back(1./a-1.);
		d.push_back(da*sum/6.);
	}
	z2dFunc = new LinearFuncX<cosmo_time >(z,d,0,d[Accuracy]);
	d2zFunc = new LinearFuncX<cosmo_time >(d,z,0,z[Accuracy]);
}

void Cosmology::print() const
{
	std::cout << "#z\td\n";
	for(int i=0; i<Accuracy; i++)
	{
		std::cout << z[i] << "\t" << d[i]*units.Lunit/Units::Mpc_in_cm << "\n";
	}
	std::cout << "\n\n";
	std::cout << "#z\tt\n";
	for(int j=0; j<Accuracy; j++)
	{
		std::cout << zt[j] << "\t" << tt[j]*units.Lunit/Units::Mpc_in_cm << "\n";
	}
}

void Cosmology::UnitTest()
{
	cosmology.Init(10);
	cosmo_time z=1;
	Particle p(Photon, z);
	std::cout << "z=" << z << "\tt=" << p.Time.t()*units.Lunit/Units::Mpc_in_cm << " Mpc d="
			<< cosmology.z2d(z)*units.Lunit/Units::Mpc_in_cm << " Mpc\n\n";
	std::cout << "Age of universe " << cosmology.z2t(0)*units.Lunit/Units::Mpc_in_cm << " Mpc\n\n";
	cosmology.print();
}

void Cosmology::initZ2T()
{
	ASSERT(fMaxZ);// must be initialized

	zt.push_back(0);
	tt.push_back(0);

	cosmo_time dz = fMaxZ/Accuracy;
	cosmo_time dz2 = 0.5*dz;
	cosmo_time z = dz;
	cosmo_time prev = diffT(0.);
	cosmo_time sum = 0.;

	for(int i=1; i<=Accuracy; i++, z += dz)
	{
		sum += prev;
		sum += 4.*diffT(z-dz2);
		prev = diffT(z);
		sum += prev;
		zt.push_back(z);
		tt.push_back(dz*sum/6.);
	}
	z2tFunc = new LinearFuncX<cosmo_time>(zt,tt,0,(tt)[Accuracy]);
	t2zFunc = new LinearFuncX<cosmo_time>(tt,zt,0,(zt)[Accuracy]);
}

///The function is used to test red shift evolution in special case of the source with no evolution in comoving frame and alpha = 2
cosmo_time Cosmology::TestIntegral(cosmo_time r /*r===min(Emax/E, 1+Zmax)*/) const
{
	return TestIntegralF(r)-TestIntegralF(1.);
}

cosmo_time Cosmology::Hy2F1(cosmo_time a, cosmo_time b, cosmo_time c, cosmo_time z)
{
	if(fabs(z)<=1)
		return gsl_sf_hyperg_2F1(a,b,c,z);
	z = 1./z;
	cosmo_time mult1 = gsl_sf_gamma(b-a)*gsl_sf_gamma(c)/gsl_sf_gamma(b)/gsl_sf_gamma(c-a);
	mult1 *= pow(-z, a);
	cosmo_time mult2 = gsl_sf_gamma(a-b)*gsl_sf_gamma(c)/gsl_sf_gamma(a)/gsl_sf_gamma(c-b);
	mult2 *= pow(-z, b);
	cosmo_time f1 = gsl_sf_hyperg_2F1(a,a-c+1,a-b+1,z);
	cosmo_time f2 = gsl_sf_hyperg_2F1(b,b-c+1,b-a+1,z);
	cosmo_time result = mult1*f1 + mult2*f2;
	return result;
}

///The function is used to test red shift evolution in special case of the source with no evolution in comoving frame and alpha = 2
cosmo_time Cosmology::TestIntegralF(cosmo_time a/*a===1+z*/) const //F[1+z, alpha] = (1+z)^(-alpha)*(Lm*(1+z)^3 + Lv)^(-1/2)
{//TestIntegralF === Integrate[F[1+z, 1], z]
	cosmo_time a3 = a*a*a;
	cosmo_time arg = -(a3*Lm)/Lv;
	cosmo_time Hypergeometric2F1 = Hy2F1(0.6666666666666666,0.5,1.6666666666666667,
         arg);
	cosmo_time result = (sqrt(1 + (-1 + a3)*Lm)*
     (4*(-1 + Lm)*sqrt((-1 + Lm - a3*Lm)/
          (-1 + Lm)) +
       a3*Lm*Hypergeometric2F1))/(4.*a*Lv*Lv*
     sqrt((-1 + Lm - a3*Lm)/(-1 + Lm)));
	return result;
}

void Cosmology::SuppressionTest() const
{
	std::cerr << "red shift spectrum suppression test for m-alpha = -2 and Zmax = " << fMaxZ << ":\n";
//here we simulate how spectrum should be suppressed at highest energies due to red shifting
//The analytic estimate was made for spectrum Q(E,z) = E^-alpha * (1+z)^(3+m) * Theta(Emax-E), for m-alpha = -2
	cosmo_time norm = 1./TestIntegral(1.+fMaxZ);
	for(cosmo_time r=1; r<fMaxZ + 1; r+=0.05)
		std::cerr << r << "\t" << TestIntegral(r)*norm << "\n";
}

} /* namespace mcray */
