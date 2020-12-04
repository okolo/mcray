/*
 * ElmagTest.cpp
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

#include "ElmagTest.h"
#include "Units.h"
#include "Utils.h"
#include "ParticleStack.h"
#include "PropagationEngine.h"
#include "Thinning.h"
#include "TableBackgrounds.h"
#include "mpfrc++/mpreal.h"
#include "Deflection1D.h"

using namespace mcray;
using namespace Utils;
using namespace Backgrounds;
using namespace Interactions;

#ifdef NO_FORTRAN
void elmag_init_test_()
{
	throw "Elmag integration is disabled";
}

#define rate_bb_(e0,zz,icq) ASSERT(FALSE)
#define sample_photon_(e0,zz,sgam,ierr) ASSERT(FALSE)
#define zpair_(sgam)_ ASSERT(FALSE)
#define sample_electron_(e0,zz,sgam,ierr) ASSERT(FALSE)
#define zics_(e0,sgam) ASSERT(FALSE)
#define zloss_(e0,zz) ASSERT(FALSE)

#else

extern "C" void elmag_init_test_();

extern "C" double rate_bb_(double* e0,double* zz,int* icq);

extern "C" void sample_photon_(double* e0,double* zz,double* sgam,int* ierr);     // sample c.m. energy for PP interaction

extern "C" double zpair_(double* sgam);  // sample secondary e+ or e- energy fraction for PP interaction

extern "C" void sample_electron_(double* e0,double* zz,double* sgam,int* ierr);     // sample c.m. energy for ICS interaction

extern "C" double zics_(double* e0,double* sgam);  // sample secondary electron(positron) energy fraction for ICS interaction

extern "C" double zloss_(double* e0,double* zz); // relative energy loss (per cm) (due to under-threshold ICS photons)

//extern "C" void get_zics_consts_(double* e0,double* sgam,double* mm,double* emin,double* zmin,double* zmax);

#endif

namespace Test {

double ElmagProxy::Rate(const Particle& aParticle)
{
	if(!fEnabled)
		Exception::Throw("ElmagTest is not initialized");
	double e0 = aParticle.Energy*1e6*units.Eunit;
	double zz = aParticle.Time.z();
	int icq = ParticleTypeToInt(aParticle.Type);
	double rate_cm = rate_bb_(&e0,&zz,&icq);
	return rate_cm*units.Lunit;
}

double ElmagProxy::RateICScel(const Particle& aParticle)
{
	if(aParticle.Type != Electron && aParticle.Type != Positron)
		return 0.;
	if(!fEnabled)
		Exception::Throw("ElmagTest is not initialized");

	double e0 = aParticle.Energy*1e6*units.Eunit;
	double zz = aParticle.Time.z();
	double rate_cm = zloss_(&e0,&zz);
	return rate_cm*units.Lunit;
}

int ElmagProxy::ParticleTypeToInt(ParticleType aType)
{
	switch(aType)
	{
	case Electron: return -1;
	case Positron: return 1;
	case Photon: return 0;
	default:
		Exception::Throw("ElmagTest:ParticleTypeToInt unexpected argument");
		return -100;
	}
}

bool ElmagProxy::SampleS(const Particle& aParticle, double& aS)
{
	if(!fEnabled)
		Exception::Throw("ElmagTest is not initialized");
	double e0 = aParticle.Energy*1e6*units.Eunit;
	double zz = aParticle.Time.z();
	double sgam = 0.;
	int ierr = 0.;
	if(aParticle.Type == Photon)
	{
		sample_photon_(&e0,&zz,&sgam,&ierr);
		if(ierr)
		{
			std::cerr << "Elmag::sample_photon exited with error code " << ierr;
			return false;
		}
	}
	else if(aParticle.Type == Positron || aParticle.Type == Electron)
	{
		sample_electron_(&e0,&zz,&sgam,&ierr);
		if(ierr)
		{
			std::cerr << "Elmag::sample_electron exited with error code " << ierr;
			return false;
		}
	}
	aS = sgam*1e-12/units.Eunit/units.Eunit;
	Debug::PrintLine(ToString(aS));
	return true;
}

bool ElmagProxy::GetSecondaries(const Particle& aParticle, std::vector<Particle>& aSecondaries, double aS)
{
	if(!fEnabled)
		Exception::Throw("ElmagTest is not initialized");
	double e0 = aParticle.Energy*1e6*units.Eunit;
	double sgam = aS/(1e-12/units.Eunit/units.Eunit);
	if(aParticle.Type == Photon)
	{
		double r = zpair_(&sgam);
		Particle electron = aParticle;
		electron.Type = Electron;
		electron.Energy = aParticle.Energy*r;
		aSecondaries.push_back(electron);
		Particle positron = aParticle;
		positron.Type = Positron;
		positron.Energy = aParticle.Energy*(1.-r);
		aSecondaries.push_back(positron);
	}
	else if(aParticle.Type == Positron || aParticle.Type == Electron)
	{
		/*
		{
		const double Esec_min = 1.;//MeV
		double yMax = 1.-Esec_min/aParticle.Energy;
		double yMin = aParticle.Mass()*aParticle.Mass()/aS;
		ASSERT(yMin<=1);
		ASSERT(yMax>=yMin);
		}
		{
			double mm,emin,zmin,zmax;
			get_zics_consts_(&e0,&sgam,&mm,&emin,&zmin,&zmax);
			if(zmin>zmax)
			{
				mpfr::mpreal Esec_min = "1000000";//MeV
				mpfr::mpreal E = e0;
				mpfr::mpreal one = "1";
				mpfr::mpreal m = aParticle.Mass()*1000000;
				mpfr::mpreal s = sgam;

				mpfr::mpreal yMax = one-Esec_min/E;
				mpfr::mpreal yMin = m*m/s;
				ASSERT(yMin<=one);
				ASSERT(yMax>=yMin);
			}
		}
		*/
		double r = zics_(&e0,&sgam);
		Particle lepton = aParticle;
		lepton.Energy = aParticle.Energy*r;
		aSecondaries.push_back(lepton);
		Particle photon = aParticle;
		photon.Type = Photon;
		photon.Energy = aParticle.Energy*(1.-r);
		aSecondaries.push_back(photon);
	}
	return true;
}

bool ElmagProxy::GetSecondaries(const Particle& aParticle, std::vector<Particle>& aSecondaries)
{
	double s=0;
	if(!SampleS(aParticle, s))
		return false;
	return GetSecondaries(aParticle, aSecondaries, s);
}

void ElmagProxy::Enable()
{
	if(!fEnabled)
	{
		elmag_init_test_();
		fEnabled = true;
	}
}

bool ElmagProxy::fEnabled = false;

int ElmagTests::Main(int argc, char** argv)
{
	if(argc==0)
	{
		RatesTest();
		return 0;
	}
	const char* command = argv[0];
	if(!strcmp(command,"Kachelries"))
	{
		int modePP = ElmagCalculateRate|ElmagSampleS|ElmagSampleSecondaryEnergy;
		int modeICS = ElmagCalculateRate|ElmagSampleS|ElmagSampleSecondaryEnergy|ElmagCalculateCelRate;
		if(argc>1)
			sscanf(argv[1], "%d", &modePP);
		if(argc>2)
			sscanf(argv[2], "%d", &modeICS);
		KachelriesTest(0.02, modePP, modeICS);
		//KachelriesTest(0.15);
	}
	return 0;
}

std::string ElmagTests::HumanReadableMode(int aMode)
{
	std::string result = "elmag";
	if(aMode == ElmagDoNotUse)
		result += "None";
	if(aMode & ElmagCalculateRate)
		result += "Rate";
	if(aMode & ElmagSampleS)
		result += "S";
	if(aMode & ElmagSampleSecondaryEnergy)
		result += "E";
	if(aMode & ElmagCalculateCelRate)
		result += "Cel";
	return result;
}

void ElmagTests::KachelriesTest(double aZmax, int aElmagUsageModePP, int aElmagUsageModeICS)
{
	if(!cosmology.IsInitialized())
		cosmology.init(0.73, 71, 10);
	double logStepK = pow(10,0.05);
	double Emin = 1./units.Eunit;
	double alphaThinning = 0.5;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	double Bgauss = 1e-17;//Gauss
	double LcorMpc = 1.0;//Mpc
	ParticleStack particles;
	Result result(Emin);

	result.AddOutput(new SpectrumOutput(
			"Elmag_Kachelries_Zmax" + ToString(aZmax) + "_PP" + HumanReadableMode(aElmagUsageModePP) + "_ICS" + HumanReadableMode(aElmagUsageModeICS),
			Emin, pow(10,0.05)));

	CompoundBackground backgr;

	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	backgr.AddComponent(new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp));
	backgr.AddComponent(new Kneiske0309EBL());

	double stepZ = aZmax<0.2 ? aZmax/2 : 0.1;
	double epsRel = 1e-3;

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, aZmax, epsRel);

	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	pe.AddInteraction(new RedShift());

	if(aElmagUsageModeICS&ElmagCalculateCelRate)
		pe.AddInteraction(new ElmagIcsCEL());
	else
		pe.AddInteraction(new Interactions::IcsCEL(backgrI,Emin));

	pe.AddInteraction(new ElmagPP(aElmagUsageModePP, new Interactions::GammaPP(backgrI)));
	pe.AddInteraction(new ElmagIcs(aElmagUsageModeICS, new Interactions::ICS(backgrI,Emin)));
	pe.AddInteraction(new Deflection1D(Bgauss, LcorMpc));

	int nParticles = 1000;
	for(int i=0; i<nParticles; i++)
	{
		Particle gamma(Photon, aZmax);
		gamma.Energy = 1e8/units.Eunit;//MeV
		particles.Add(gamma);
	}

	pe.Run();
	//pe.RunMultithread();
}

void ElmagTests::RatesTest()
{
	if(!cosmology.IsInitialized())
		cosmology.init(0.73, 71, 10);

	double logStepK = pow(10,0.05);

	RedShift rs;
	ElmagIcs ics;
	ElmagIcsCEL icsCel;
	ElmagPP pp;
	double Emin=1;
	double Emax=1e10;

	std::ofstream ratesOut;
	ratesOut.open("elmag_rates",std::ios::out);

	ratesOut << "#E\tRedshiftRate\tIcsRate\tIcsCelRate\tPpRate\t[E]=1eV [Rate]=1/Mpc\n";
	Particle particle(Electron, 0.);

	for(double E=Emin; E<Emax; E*=logStepK)
	{
		particle.Energy = E/units.Eunit;//MeV
		particle.Type = Electron;
		double celRate = icsCel.Rate(particle);
		double icsRate = ics.Rate(particle);
		double redshiftRate = rs.Rate(particle);
		particle.Type = Photon;
		double ppRate = pp.Rate(particle);

		double mult = 1./units.Mpc;
		ratesOut << E*1e6 <<  "\t" << redshiftRate*mult <<  "\t" << icsRate*mult <<
				"\t" << celRate*mult <<  "\t" << ppRate*mult <<	std::endl;
	}
}

}//namespace Test {
