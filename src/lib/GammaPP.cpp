/*
 * GammaPP.cpp
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

#include "GammaPP.h"
#include <math.h>
#include "MathUtils.h"
#include "Output.h"
#include "ParticleStack.h"
#include "TableBackgrounds.h"
#include "PropagationEngine.h"

namespace Interactions {
using namespace mcray;

const double EMInteraction::AlphaEM = 7.2973525698e-3;//~1/137
const double EMInteraction::SigmaCoef = AlphaEM*AlphaEM*2*M_PI;

GammaPP::GammaPP(BackgroundIntegral* aBackground):
		fBackground(aBackground)
{
}

RandomInteraction* GammaPP::Clone() const
{
	return new GammaPP(fBackground->Clone());
}

GammaPP::~GammaPP() {
}

double GammaPP::SampleSecondaryFracE(double aS, double aRand) const
{
	double mass = Particle::Mass(Electron);
	double s = aS/(mass*mass);
	SecondaryEnergySamplingEquationPars pars = {s, aRand*DifSigmaIntegral(0.5, s) };
    double x_lo = 0.5*(1. - sqrt(1. - 4./s));
    double x_hi = 0.5;
    const double relError = 1e-3;
    const int max_iter = 100;
    double r = 0;
    try{
    	//TODO: use MathUtils::SampleLogDistribution sampling method
    	r = MathUtils::SolveEquation(&SecondaryEnergySamplingEquation, x_lo, x_hi, &pars, relError, max_iter, gsl_root_fsolver_bisection);
    }catch(Exception* ex)
    {
#ifdef _DEBUG
    	fDebugOutput = true;
    	double dx = 0.01*(x_hi-x_lo);
    	for(double x=x_lo; x<=x_hi; x+=dx)
    	{
    		SecondaryEnergySamplingEquation(x, &pars);
    	}
    	ASSERT(0);
#endif
    	throw ex;
    }
    return r;
}

bool GammaPP::SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const
{
	ASSERT(aParticle.Type == Photon);
	aS = 0.;
	return (fBackground->GetRateAndSampleS(fSigma, aParticle, aRandomizer, aS)!=0);
}

void GammaPP::SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const
{
	ASSERT(aParticle.Type == Photon);
	Particle sec = aParticle;
	Particle sec2 = aParticle;
	double rand = aRandomizer.Rand();
	sec.Type = rand>0.5 ? Electron : Positron;
	sec2.Type = rand>0.5 ? Positron : Electron;
	double r = SampleSecondaryFracE(aS, aRandomizer.Rand());
	sec.Energy = aParticle.Energy*r;
	sec2.Energy = aParticle.Energy - sec.Energy;
	aSecondaries.push_back(sec);
	aSecondaries.push_back(sec2);
}

GammaPP::Sigma::Sigma()
{
	double m = Particle::Mass(Electron);
	fThresholdS=4.*m*m;
}

double GammaPP::Sigma::f(double s) const
{
	double beta2 = 1. - fThresholdS/s;
	if(beta2<=0)
		return 0;
	double beta = sqrt(beta2);

	double result = SigmaCoef*((3.-beta2*beta2)*log((1.+beta)/(1.-beta))-2.*beta*(2.-beta2))/s;
	ASSERT_VALID_NO(result);
	return result;
}

double GammaPP::Rate(const Particle& aParticle) const
{
	if(aParticle.Type != Photon)
		return 0.;
	return fBackground->GetRateS(fSigma, aParticle);
}

bool GammaPP::fDebugOutput = false;

double GammaPP::SecondaryEnergySamplingEquation(double r, void* aParams)
{
	SecondaryEnergySamplingEquationPars* p=(SecondaryEnergySamplingEquationPars*)aParams;
	double result = DifSigmaIntegral(r, p->s) - p->Integral;
#ifdef _DEBUG
	if(GammaPP::fDebugOutput)
	{
		std::cerr << "f( " << r << " ) = " << result << std::endl;
	}
#endif
	return result;
}

double GammaPP::DifSigmaIntegral(double r, double s)
{
	double beta = sqrt(1.-4./s);
	double rMin = 0.5*(1.0 - beta);
	if(r<=rMin)
		return 0.;
	return DifSigmaIntegralUnnorm(r, s) - DifSigmaIntegralUnnorm(rMin, s);
}

double GammaPP::DifSigmaIntegralUnnorm(double r, double s)
{
/*
Formula for unnormalized diff cross section was taken from http://lanl.arxiv.org/abs/1106.5508v1 eq. (10)
DifSigmaIntegral is integral of diff cross section from r_min = 1/2*(1 - (1 - 4/s)^(1/2)) to r
Mathematica program used to produce the output:
G[r_, s_] = 1/r*(r^2/(1 - r) + 1 - r + 4/s/(1 - r) - 4/s/s/r/(1 - r)/(1 - r))/(1 + (2 - 8/s)*4/s)
F[r_, s_] = Integrate[G[t, s], {t, 1/2*(1 - (1 - 4/s)^(1/2)), r}, Assumptions -> r > 1/2*(1 - (1 - 4/s)^(1/2)) && r < 0.5 && s > 4]
//the exact expression below has small (~1e-16 in mathematica) nonzero value for r=rMin (mathematica bug?)
//therefore we subtract F[rMin, s] from F[r, s]
*/
	double s2 = s*s;
	double r2 = r*r;
	double r3 = r2*r;
	double beta = sqrt(1.-4./s);
	double log_r_minus = log(1 - r);
	double log_r = log(r);
	double result = 0.;
	/*
	///test (moved to DifSigmaIntegralTest)
	double rMin = 0.5*(1.0 - beta);
	if(r<=rMin)
		return 0.;
	double a=r-rMin;
	if(a<1e-3)
	{
		double b_1=beta-1.;
		double b_P1 = beta+1.;
		result = (16*a*(-2 + s))/
			    (b_1*b_1*b_P1*b_P1*
			      (-32 + 8*s + s2)) +
			   (32*a*a*(-4*beta + beta*s))/
			    (b_1*b_1*b_1*b_P1*b_P1*b_P1*
			      (-32 + 8*s + s2)) +
			   (256*a*a*a*(28 - 11*s + s2))/
			    (3.*b_1*b_1*b_1*b_1*b_P1*b_P1*b_P1*b_P1*s*
			      (-32 + 8*s + s2));
	}
	else
	{*/
		result = -((4 - 8*r + r*s2 - 3*r2*s2 +
	       2*r3*s2 - 4*r*beta +
	       4*r2*beta - r*s2*beta +
	       r2*s2*beta + 8*r*log_r_minus -
	       8*r2*log_r_minus - 4*r*s*log_r_minus +
	       4*r2*s*log_r_minus - r*s2*log_r_minus +
	       r2*s2*log_r_minus - 8*r*log_r + 8*r2*log_r +
	       4*r*s*log_r - 4*r2*s*log_r + r*s2*log_r -
	       r2*s2*log_r +
	       (-1 + r)*r*(-8 + 4*s + s2)*log(1 - beta) -
	       (-1 + r)*r*(-8 + 4*s + s2)*log(1 + beta))/
	     ((-1 + r)*r*(-32 + 8*s + s2)));
	//}
	ASSERT_VALID_NO(result);
	return result;
}

double GammaPP::DifSigmaIntegralTest(double r, double s, bool aSeries, double& a)
{
/*
Formula for unnormalized diff cross section was taken from http://lanl.arxiv.org/abs/1106.5508v1 eq. (10)
DifSigmaIntegral is integral of diff cross section from r_min = 1/2*(1 - (1 - 4/s)^(1/2)) to r
Mathematica program used to produce the output:
G[r_, s_] = 1/r*(r^2/(1 - r) + 1 - r + 4/s/(1 - r) - 4/s/s/r/(1 - r)/(1 - r))/(1 + (2 - 8/s)*4/s)
F[r_, s_] = Integrate[G[t, s], {t, 1/2*(1 - (1 - 4/s)^(1/2)), r}, Assumptions -> r > 1/2*(1 - (1 - 4/s)^(1/2)) && r < 0.5 && s > 4]
//the exact expression below has small (~1e-16 in mathematica) nonzero value for r=rMin (mathematica bug?)
*/
	double s2 = s*s;
	double r2 = r*r;
	double r3 = r2*r;
	double beta = sqrt(1.-4./s);
	double log_r_minus = log(1 - r);
	double log_r = log(r);
	double rMin = 0.5*(1.0 - beta);
	if(r<=rMin)
		return 0.;
	a=r-rMin;
	double result = 0.;
	if(aSeries)///Series[F[1/2*(1 - (1 - 4/s)^(1/2)) + a, s], {a, 0, 3}] excluding part depending on s only
	{
		double b_1=beta-1.;
		double b_P1 = beta+1.;
		result = (16*a*(-2 + s))/
			    (b_1*b_1*b_P1*b_P1*
			      (-32 + 8*s + s2)) +
			   (32*a*a*(-4*beta + beta*s))/
			    (b_1*b_1*b_1*b_P1*b_P1*b_P1*
			      (-32 + 8*s + s2)) +
			   (256*a*a*a*(28 - 11*s + s2))/
			    (3.*b_1*b_1*b_1*b_1*b_P1*b_P1*b_P1*b_P1*s*
			      (-32 + 8*s + s2));
	}
	else
	{//F[r_, s_]
		result = -((4 - 8*r + r*s2 - 3*r2*s2 +
	       2*r3*s2 - 4*r*beta +
	       4*r2*beta - r*s2*beta +
	       r2*s2*beta + 8*r*log_r_minus -
	       8*r2*log_r_minus - 4*r*s*log_r_minus +
	       4*r2*s*log_r_minus - r*s2*log_r_minus +
	       r2*s2*log_r_minus - 8*r*log_r + 8*r2*log_r +
	       4*r*s*log_r - 4*r2*s*log_r + r*s2*log_r -
	       r2*s2*log_r +
	       (-1 + r)*r*(-8 + 4*s + s2)*log(1 - beta) -
	       (-1 + r)*r*(-8 + 4*s + s2)*log(1 + beta))/
	     ((-1 + r)*r*(-32 + 8*s + s2)));
	}
	ASSERT_VALID_NO(result);
	return result;
}

void GammaPP::UnitTest()
{
	double Zmax = 2e-6;
	double Emin = 10/units.Eunit;
	//double alphaThinning = 0;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	result.AddOutput(new SpectrumOutput("out", Emin, pow(10, 0.05)));

	cosmology.Init();

	CompoundBackground backgr;
	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;

	backgr.AddComponent(new PlankBackground(cmbTemp, 1e-4*cmbTemp, 1e4*cmbTemp));
	//backgr.AddComponent(new PlankBackground(cmbTemp, 6.2e-10/units.Eunit, 6.4e-10/units.Eunit));//only central value
	//backgr.AddComponent(new Backgrounds::Kneiske1001EBL());
	//backgr.AddComponent(new Backgrounds::GaussianBackground(6.3e-10/units.Eunit, 1e-11/units.Eunit, 1., 413*units.Vunit));
	//double n = BackgroundUtils::CalcIntegralDensity(backgr)/units.Vunit;

	double stepZ = 0.05;
	double logStepK = pow(10,0.05);
	double epsRel = 1e-3;

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel);
	GammaPP* pp = new GammaPP(backgrI);
	RedShift* rs = new RedShift();
	double E = Emin;
	std::cerr << "#E\tRedshiftRate\tPPRate\t\t[E]=1eV [Rate]=1/Mpc\n";

	for(double E=1; E<1e14; E*=logStepK)
	{
		Particle photon(Photon, 0.);
		photon.Energy = E/units.Eunit;//MeV
		double rate = pp->Rate(photon);
		double redshiftRate = rs->Rate(photon);
		double mult = 1./units.Lunit*Units::Mpc_in_cm;
		std::cerr << E*1e6 <<  "\t" << redshiftRate*mult <<  "\t" << rate*mult <<	std::endl;
	}
	double s = 8.;
	double x_lo = 0.5*(1. - sqrt(1. - 4./s));
	double dx = 0.01*(0.5-x_lo);
	std::cerr << "\n\n\n\n\n";
	for(double x=x_lo; x<=0.5; x+=dx)
	{
		double a;
		double seriesRes = DifSigmaIntegralTest(x,s,true,a);
		double origRes = DifSigmaIntegralTest(x,s,false,a);
		std::cerr << a << "\t" << origRes << "\t" << seriesRes << std::endl;
	}


	PropagationEngine pe(particles, result, 2011);
	//EnergyBasedThinning thinning(alphaThinning);
	//pe.SetThinning(&thinning);
	pe.AddInteraction(rs);
	pe.AddInteraction(pp);

	int nParticlesPerBin = 10000;
	double Emax = 1e9/units.Eunit;
	/// build initial state
	//for(double E=Emin; E<Emax; E*=logStepK)
	//{
	E=Emax;
		Particle photon(Photon, Zmax);
		photon.Energy = E;//MeV
		//electron.Weight = exp(-Emax/E);
		//if(electron.Weight>0)
		for(int i=0; i<nParticlesPerBin; i++)
			particles.AddPrimary(photon);
	//}
	std::cerr << std::setprecision(15);
	pe.Run();
}

void GammaPP::UnitTestSampling(double aE, double aZ)
{
	std::string out = "sample_ph_E";
	out += ToString(aE*1e6);//eV
	out += "_Z";
	out += ToString(aZ);
	std::ofstream secPPout;
	secPPout.open("sample_sec_pp",std::ios::out);
	std::ofstream secPPoutVarS;
	secPPoutVarS.open("sample_sec_pp_varS",std::ios::out);

	debug.SetOutputFile(out);

	if(!cosmology.IsInitialized())
		cosmology.Init();

	CompoundBackground backgr;
	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;

	backgr.AddComponent(new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp));
	//backgr.AddComponent(new PlankBackground(cmbTemp, 6.2e-10/units.Eunit, 6.4e-10/units.Eunit));//only central value

	double stepZ = 0.05;
	double logStepK = pow(10,0.05);
	double epsRel = 1e-3;

	Backgrounds::Kneiske0309EBL* kn0309 = new Backgrounds::Kneiske0309EBL();

	backgr.AddComponent(kn0309);

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, aZ+1, epsRel);
	Interactions::GammaPP pp(backgrI);
	Particle photon(Photon, aZ);
	photon.Energy = aE;//MeV
	photon.Time.setZ(aZ);
	photon.Type = Photon;
	Randomizer rnd;
	double S=2.;
	for(int i=0; i<10000; i++)
	{
		std::vector<Particle> secondaries;
		if(pp.GetSecondaries(photon, secondaries, rnd))//this outputs sampled s to debug file
		{
			ASSERT(secondaries.size()==2);
			secPPoutVarS << secondaries[0].Energy/aE << "\n" << secondaries[1].Energy/aE << "\n";
		}
		double r = pp.SampleSecondaryFracE(S, rnd.Rand());
		secPPout << r << '\n';
		secPPout << 1.-r << std::endl;
	}
	secPPout.close();
}

void GammaPP::UnitTest(bool aPrintRates, bool aPropagate, bool aSplit, bool aBestFitEBL)
{
	double Zmax = 2e-4;
	double Emax = 3e8/units.Eunit;
	double Emin = 1/units.Eunit;
	//double alphaThinning = 0;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	if(aPropagate)
	{
		std::string out = "out_pp_";
		out += aSplit?"split" : "no_split";
		result.AddOutput(new SpectrumOutput(out, Emin, pow(10, 0.05)));
	}
	if(!cosmology.IsInitialized())
		cosmology.Init();

	CompoundBackground backgr;
	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;

	backgr.AddComponent(new PlankBackground(cmbTemp, 1e-5*cmbTemp, 1e3*cmbTemp));
	//backgr.AddComponent(new PlankBackground(cmbTemp, 6.2e-10/units.Eunit, 6.4e-10/units.Eunit));//only central value

	double stepZ = 0.05;
	double logStepK = pow(10,0.05);
	double epsRel = 1e-3;

	TableBackground* ebl = aBestFitEBL ? ((TableBackground*)new Backgrounds::Kneiske0309EBL()) :
			((TableBackground*)new Backgrounds::Kneiske1001EBL());
	PBackgroundIntegral backgrIebl;
	PRandomInteraction ppEbl;
	if(!aSplit)
	{
		backgr.AddComponent(ebl);
	}
	else
	{
		backgrIebl = new ContinuousBackgroundIntegral(*ebl, stepZ, logStepK, Zmax, epsRel);
		ppEbl = new GammaPP(backgrIebl);
	}

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel);
	PRandomInteraction pp = new GammaPP(backgrI);
	PCELInteraction rs = new RedShift();
	double E = Emin;
	if(aPrintRates)
	{
		std::ofstream ratesOut;
		std::string ratesFile = "pp_rates_";
		ratesFile += aBestFitEBL ? "BestFitEBL_" : "MinimalEBL_";
		ratesFile += aSplit ? "split" : "no_split";
		ratesOut.open(ratesFile.c_str(),std::ios::out);

		ratesOut << "#E\tRedshiftRate\tPPRate\t\t[E]=1eV [Rate]=1/Mpc\n";

		for(double E=Emin; E<Emax; E*=logStepK)
		{
			Particle photon(Photon, 0.);
			photon.Energy = E/units.Eunit;//MeV
			double ppRate = pp->Rate(photon);
			double redshiftRate = rs->Rate(photon);
			if(aSplit)
			{
				ppRate += ppEbl->Rate(photon);
			}
			double mult = 1./units.Lunit*Units::Mpc_in_cm;
			ratesOut << E*1e6 <<  "\t" << redshiftRate*mult <<  "\t" << ppRate*mult << std::endl;
		}
	}
	if(!aPropagate)
		return;

	PropagationEngine pe(particles, result, 2011);
	//EnergyBasedThinning thinning(alphaThinning);
	//pe.SetThinning(&thinning);
	//pe.AddInteraction(rs);
	pe.AddInteraction(pp);
	//ics->fDeleteElectron = true;//check in single interaction mode
	if(aSplit)
	{
		pe.AddInteraction(ppEbl);
	}

	int nParticlesPerBin = 10000;
	/// build initial state
	//for(double E=Emin; E<Emax; E*=logStepK)
	//{
	E=Emax;
		Particle photon(Photon, Zmax);
		photon.Energy = E;//MeV
		//electron.Weight = exp(-Emax/E);
		//if(electron.Weight>0)
		for(int i=0; i<nParticlesPerBin; i++)
			particles.AddPrimary(photon);
	//}

	pe.RunMultithread();
	if(aSplit)
		delete ebl;
}

} /* namespace Interactions */

