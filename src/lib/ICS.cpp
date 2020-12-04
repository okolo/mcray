/*
 * ICS.cpp
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
#include "ICS.h"
#include <math.h>
#include "MathUtils.h"
#include "Output.h"
#include "ParticleStack.h"
#include "TableBackgrounds.h"
#include "PropagationEngine.h"
#include "PrecisionTests.h"

namespace Interactions {
using namespace mcray;

ICS::ICS(BackgroundIntegral* aBackground, double aGammaEmin):
		fBackground(aBackground),
		fGammaEmin(aGammaEmin),
		fDeleteElectron(false)
{
}

ICS::~ICS() {
}

double ICS::SampleSecondaryGammaFracE(double aS, double rMin, double rand) const
{
	double mass = Particle::Mass(Electron);
	double s = aS/(mass*mass);
	SecondaryEnergySamplingEquationPars pars;
	pars.s = s;
	pars.Integral = rand*PartialSigma(rMin, s);
    double x_lo = rMin;
    double x_hi = 1.-1./s;
    const double relError = 1e-3;
    const int max_iter = 100;
    double r = MathUtils::SolveEquation(&SecondaryEnergySamplingEquation, x_lo, x_hi, &pars, relError, max_iter, gsl_root_fsolver_bisection);
    return r;
}

bool ICS::SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const
{
	ASSERT(aParticle.Type == Electron || aParticle.Type == Positron);
	double rMin = fGammaEmin/aParticle.Energy;
	Sigma sigma(rMin);
	aS = 0.;
	return (fBackground->GetRateAndSampleS(sigma, aParticle, aRandomizer, aS)!=0);
}

void ICS::SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const
{
	double rMin = fGammaEmin/aParticle.Energy;
	double r = SampleSecondaryGammaFracE(aS, rMin, aRandomizer.Rand());
	Particle secLepton = aParticle;
	Particle secPhoton = aParticle;
	secPhoton.Type = Photon;
    secPhoton.Energy = aParticle.Energy*r;
    secLepton.Energy = aParticle.Energy - secPhoton.Energy;
	aSecondaries.push_back(secPhoton);
	if(!fDeleteElectron)
		aSecondaries.push_back(secLepton);
}

ICS::Sigma::Sigma(double rMin) : f_rMin(rMin)
{
	double mass = Particle::Mass(Electron);
	fM2 = mass*mass;
	fXmin = fM2/(1.-rMin);
}

double ICS::Sigma::f(double s) const
{
	return PartialSigma(f_rMin, s/fM2);
}

double ICS::Rate(const Particle& aParticle) const
{
	if(aParticle.Type != Electron && aParticle.Type != Positron)
		return 0.;
	double rMin = fGammaEmin/aParticle.Energy;
	if(rMin>=1)
		return 0.;
	Sigma sigma(rMin);
	return fBackground->GetRateS(sigma, aParticle);
}

double ICS::SecondaryEnergySamplingEquation(double r, void* aParams)
{
	SecondaryEnergySamplingEquationPars* p=(SecondaryEnergySamplingEquationPars*)aParams;
	return PartialSigma(r, p->s) - p->Integral;
}

double ICS::PartialSigma(double rMin, double s)
{
/*
Formula for partial cross section was taken from http://lanl.arxiv.org/abs/1106.5508v1 eq. (9)
*/
	double m=Particle::Mass(Electron);

	//return SigmaCoef/m/m*Test::PrecisionTests::PartialSigmaAccurate(rMin,s);


	double yMin=1./s;//s here is in units of m^2
	double yMax = 1-rMin;
	if(yMin>=yMax)
	{
		ASSERT((yMax-yMin)/(yMax+yMin)<1e-10);
		return 0.;
	}

	double y_1 = 1.-yMin;
	double y_1_2 = y_1*y_1;
	double a=yMax-yMin;
	double result = 0;
	if(a>1e-3)
	{
		result = yMin*(yMax-yMin)/y_1*(log(yMax/yMin)/(yMax-yMin)*(1.-4.*yMin*(1.+yMin)/y_1_2)+4.*(yMin/yMax+yMin)/y_1_2+0.5*(yMax+yMin));
		ASSERT_VALID_NO(result);
	}
	else
	{
		double x=yMin;
		double x2=x*x;
		double x3=x2*x;
		double x4=x3*x;
		double x_1=x-1.;//negative
		double x_1_2=x_1*x_1;
		double x_1_3=x_1_2*x_1;
		double a2=a*a;
		double a3=a2*a;
		double a4=a3*a;
		double a5=a4*a;

		//result = a/x_1*(-1.-x2+a/x/x_1*(-1.-3.*x+x2-x3+a/x/x_1*(-1./3.-2.*x+x2+a/x*(0.25+2.5*x-0.75*x2+0.2*a/x*(-1.-14.*x+3*x2)))));

		result = (a4*(1 + 10*x - 3*x2))/(4.*x_1_3*x3) +
		   (a*(1 + x2))/(1 - x) +
		   (a5*(-1 - 14*x + 3*x2))/(5.*x_1_3*x4) +
		   (a3*(-1 - 6*x + 3*x2))/(3.*x_1_3*x2) +
		   (a2*(-1 - 3*x + x2 - x3))/(2.*x_1_2*x);

		ASSERT_VALID_NO(result);
	}

	return SigmaCoef/m/m*result;
}

double PartialSigmaSeries(double rMin, double s)
{
/*
Formula for partial cross section was taken from http://lanl.arxiv.org/abs/1106.5508v1 eq. (9)
*/
	double yMin=1./s;//s here is in units of m^2
	double yMax = 1-rMin;
	if(yMin>=yMax)
	{
		ASSERT((yMax-yMin)/(yMax+yMin)<1e-10);
		return 0.;
	}

	double a=yMax-yMin;
	double result = 0;
	{
		double x=yMin;
		double x2=x*x;
		double x3=x2*x;
		double x4=x3*x;
		double x_1=x-1.;//negative
		double x_1_2=x_1*x_1;
		double x_1_3=x_1_2*x_1;
		double a2=a*a;
		double a3=a2*a;
		double a4=a3*a;
		double a5=a4*a;

		//result = a/x_1*(-1.-x2+a/x/x_1*(-1.-3.*x+x2-x3+a/x/x_1*(-1./3.-2.*x+x2+a/x*(0.25+2.5*x-0.75*x2+0.2*a/x*(-1.-14.*x+3*x2)))));

		result = (a4*(1 + 10*x - 3*x2))/(4.*x_1_3*x3) +
		   (a*(1 + x2))/(1 - x) +
		   (a5*(-1 - 14*x + 3*x2))/(5.*x_1_3*x4) +
		   (a3*(-1 - 6*x + 3*x2))/(3.*x_1_3*x2) +
		   (a2*(-1 - 3*x + x2 - x3))/(2.*x_1_2*x);

		ASSERT_VALID_NO(result);
	}
	return result;
}

double PartialSigmaExact(double rMin, double s)
{
/*
Formula for partial cross section was taken from http://lanl.arxiv.org/abs/1106.5508v1 eq. (9)
*/
	double yMin=1./s;//s here is in units of m^2
	double yMax = 1-rMin;
	if(yMin>=yMax)
	{
		ASSERT((yMax-yMin)/(yMax+yMin)<1e-10);
		return 0.;
	}

	double y_1 = 1.-yMin;
	double y_1_2 = y_1*y_1;
	//double a=yMax-yMin;
	double result = 0;

	result = yMin*(yMax-yMin)/y_1*(log(yMax/yMin)/(yMax-yMin)*(1.-4.*yMin*(1.+yMin)/y_1_2)+4.*(yMin/yMax+yMin)/y_1_2+0.5*(yMax+yMin));
	ASSERT_VALID_NO(result);

	return result;
}

void ICS::UnitTestDistribution(double aEmin, double aE, double aZ)
{
	std::string out = "distrib_s_ics_Emin";
	out += ToString(aEmin*1e6);//eV
	out += "_E";
	out += ToString(aE*1e6);//eV
	out += "_Z";
	out += ToString(aZ);

	std::ofstream distribOut;
	distribOut.open(out.c_str(),std::ios::out);
	FunctionCallLoggerX<double> logger(distribOut);

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
	PBackgroundIntegral backgrIebl;
	PCELInteraction icsCelEbl,icsCelFullEbl;
	PRandomInteraction icsEbl,icsFullEbl;

	backgr.AddComponent(kn0309);

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, aZ+1, epsRel);
	Interactions::ICS ics(backgrI,aEmin);
	Particle electron(Electron, aZ);
	electron.Energy = aE;//MeV
	Randomizer rnd;

	std::vector<Particle> secondaries;
	MathUtils::SetLogger(&logger);
	ics.GetSecondaries(electron, secondaries, rnd);//this outputs distribution

	distribOut.close();
}

void ICS::UnitTestSampling(double aEmin, double aE, double aZ)
{
	std::string out = "sample_e_Emin";
	out += ToString(aEmin*1e6);//eV
	out += "_E";
	out += ToString(aE*1e6);//eV
	out += "_Z";
	out += ToString(aZ);

    debug.SetOutputFile(out);

	std::ofstream secICSout;
	secICSout.open("sample_sec_gamma_ics",std::ios::out);

	std::ofstream secICSoutVarS;
	secICSoutVarS.open("sample_sec_gamma_ics_varS",std::ios::out);

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
	PBackgroundIntegral backgrIebl;
	PCELInteraction icsCelEbl,icsCelFullEbl;
	PRandomInteraction icsEbl,icsFullEbl;

	backgr.AddComponent(kn0309);

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, aZ+1, epsRel);
	Interactions::ICS ics(backgrI,aEmin);
	Particle electron(Electron,aZ);

	electron.Energy = aE;//MeV

	Randomizer rnd;
	double S=2.;
	for(int i=0; i<10000; i++)
	{
		std::vector<Particle> secondaries;
		if(ics.GetSecondaries(electron, secondaries, rnd))//this outputs sampled s to debug file
		{
			ASSERT(secondaries.size()==2);
			ASSERT(secondaries[0].Type == Photon);
			secICSoutVarS << secondaries[0].Energy/aE << "\n";
		}
		double r = ics.SampleSecondaryGammaFracE(S, aEmin/aE, rnd.Rand());
		secICSout << r << "\n";
	}
	secICSout.close();
}

void ICS::UnitTest(bool aPrintRates, bool aPropagate, bool aSplit, double Emin, double Emax, double Zmax)
{
	//double alphaThinning = 0;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	if(aPropagate)
	{
		std::string out = "ics_";
		out += aSplit?"split" : "no_split";
		out += "_Emin";
		out += ToString(Emin);
		out += "_Emax";
		out += ToString(Emax);
		out += "_Zmax";
		out += ToString(Zmax);
		result.AddOutput(new SpectrumOutput("out_" + out, Emin, pow(10, 0.05)));
        debug.SetOutputFile("debug_" + out);
	}
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
	PBackgroundIntegral backgrIebl;
	PCELInteraction icsCelEbl,icsCelFullEbl;
	PRandomInteraction icsEbl,icsFullEbl;
	if(!aSplit)
	{
		backgr.AddComponent(kn0309);
	}
	else
	{
		backgrIebl = new ContinuousBackgroundIntegral(*kn0309, stepZ, logStepK, Zmax, epsRel);
		icsCelEbl = new Interactions::IcsCEL(backgrIebl,Emin);
		icsEbl = new Interactions::ICS(backgrIebl,Emin);
		icsFullEbl = new Interactions::ICS(backgrIebl,0);
		icsCelFullEbl = new Interactions::IcsCEL(backgrIebl,1e20);
	}

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel);
	PCELInteraction icsCel = new Interactions::IcsCEL(backgrI,Emin);
	PRandomInteraction ics = new Interactions::ICS(backgrI,Emin);
	PRandomInteraction icsFull = new Interactions::ICS(backgrI,0);
	PCELInteraction icsCelFull = new Interactions::IcsCEL(backgrI,1e20);
	PCELInteraction rs = new RedShift();
	double E = Emin;
	if(aPrintRates)
	{
		std::ofstream ratesOut;
		std::string ratesFile = "ics_rates_";
		ratesFile += aSplit ? "split" : "no_split";
		ratesOut.open(ratesFile.c_str(),std::ios::out);

		ratesOut << "#E\tRedshiftRate\tIcsRate\tIcsCelRate\tIcsFullRate\tIcsCelFull\t\t[E]=1eV [Rate]=1/Mpc\n";

		for(double E=Emin; E<Emax; E*=logStepK)
		{
			Particle electron(Electron, 0.);
			electron.Energy = E/units.Eunit;//MeV
			double celRate = icsCel->Rate(electron);
			double icsRate = ics->Rate(electron);
			double redshiftRate = rs->Rate(electron);
			double icsFullRate = icsFull->Rate(electron);
			double celFull = icsCelFull->Rate(electron);
			if(aSplit)
			{
				celRate += icsCelEbl->Rate(electron);
				icsRate += icsEbl->Rate(electron);
				icsFullRate += icsFullEbl->Rate(electron);
				celFull += icsCelFullEbl->Rate(electron);
			}
			double mult = 1./units.Lunit*Units::Mpc_in_cm;
			ratesOut << E*1e6 <<  "\t" << redshiftRate*mult <<  "\t" << icsRate*mult <<
					"\t" << celRate*mult <<  "\t" << icsFullRate*mult <<
					"\t" << celFull*mult <<	std::endl;
		}
	}
	if(!aPropagate)
		return;

	PropagationEngine pe(particles, result);
	//EnergyBasedThinning thinning(alphaThinning);
	//pe.SetThinning(&thinning);
	//pe.AddInteraction(rs);
	pe.AddInteraction(ics);
	//ics->fDeleteElectron = true;//check in single interaction mode
	//pe.AddInteraction(icsCel);
	if(aSplit)
	{
		//pe.AddInteraction(icsCelEbl);
		pe.AddInteraction(icsEbl);
	}

	int nParticlesPerBin = 100000;

	/// build initial state
	//for(double E=Emin; E<Emax; E*=logStepK)
	//{
	E=Emax;
		Particle electron(Electron, Zmax);
		electron.Energy = E;//MeV
		//electron.Weight = exp(-Emax/E);
		//if(electron.Weight>0)
		for(int i=0; i<nParticlesPerBin; i++)
			particles.AddPrimary(electron);
	//}

	pe.RunMultithread();
	if(aSplit)
		delete kn0309;
}

void ICS::UnitTestMono(bool aMono, int aAc, bool aRatesOnly)
{
	if(!cosmology.IsInitialized())
		cosmology.Init();
	//double rMin = 0.00031827314821827185;//3.18e-4;
	//double s=1.0003183787358734;//1./(0.999-rMin);
	//double exactSigma=PartialSigmaExact(rMin,s);
	//double seriesSigma=PartialSigmaSeries(rMin,s);
	//double sigmaAc = Test::PrecisionTests::PartialSigmaAccurate(rMin,s);

	double Zmax = 1e-6;
	double Emin = 10/units.Eunit;
	double alphaThinning = 0;//1 conserves number of particles on the average; alpha = 0 disables thinning
	double logStepK = pow(10,0.05/aAc);
	double epsRel = 1e-3;

	double conc = 413*units.Vunit;//413 cm^{-3}
	double centralE = 6.3e-10/units.Eunit;
	ConstFunction k(centralE);
	ConstFunction c(conc);

	GaussianBackground backgr(centralE, 0.1*centralE, 1., conc);

	//double concTest = BackgroundUtils::CalcIntegralDensity(backgr, epsRel)/units.Vunit;


	double stepZ = Zmax<0.2 ? Zmax/2 : 0.1;

	PBackgroundIntegral backgrM = new MonochromaticBackgroundIntegral(c, k, logStepK, epsRel);
	PBackgroundIntegral backgrC = new ContinuousBackgroundIntegral(backgr, stepZ, pow(logStepK,0.1), Zmax, epsRel);
	BackgroundIntegral* backgrI = aMono ? backgrM : backgrC;

	Interactions::IcsCEL* icsCel = new Interactions::IcsCEL(backgrI,Emin);
	Interactions::ICS* ics = new Interactions::ICS(backgrI,Emin);
	Interactions::ICS icsFull(backgrI,0);
	Interactions::IcsCEL icsCelFull(backgrI,1e20);
	RedShift* rs = new RedShift();
	double E = Emin;
	std::ofstream ratesOut;
	std::string ratesFile = "ics_rates_";
	ratesFile += aMono ? "mono" : "gauss";
	ratesFile += "_ac" + ToString(aAc);
	ratesOut.open(ratesFile.c_str(),std::ios::out);
	ratesOut << "#E\tRedshiftRate\tIcsRate\tIcsCelRate\tIcsFullRate\tIcsCelFull\t\t[E]=1eV [Rate]=1/Mpc\n";

	for(double E=Emin; E<1e8*Emin; E*=logStepK)
	{
		Particle electron(Electron, 0.);
		electron.Energy = E/units.Eunit;//MeV
		double celRate = icsCel->Rate(electron);
		double icsRate = ics->Rate(electron);
		double redshiftRate = rs->Rate(electron);
		double icsFullRate = icsFull.Rate(electron);
		double celFull = icsCelFull.Rate(electron);
		double mult = 1./units.Lunit*Units::Mpc_in_cm;
		ratesOut << E*1e6 <<  "\t" << redshiftRate*mult <<  "\t" << icsRate*mult <<
				"\t" << celRate*mult <<  "\t" << icsFullRate*mult <<
				"\t" << celFull*mult <<	std::endl;
	}
	ratesOut.close();
	if(aRatesOnly)
		return;

	ParticleStack particles;
	Result result(Emin);
	result.AddOutput(new SpectrumOutput(aMono ? "out_ics_mono" : "out_ics_gauss", Emin, pow(10, 0.05)));
	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	pe.AddInteraction(rs);
	pe.AddInteraction(ics);
	//ics->fDeleteElectron = true;//check in single interaction mode
	pe.AddInteraction(icsCel);
	//pe.AddInteraction(new Interactions::GammaPP(backgrI));

	int nParticlesPerBin = 1000;
	double Emax = 1e6/units.Eunit;
	/// build initial state
	//for(double E=Emin; E<Emax; E*=logStepK)
	//{
	E=Emax;
		Particle electron(Electron, Zmax);
		electron.Energy = E;//MeV
		//electron.Weight = exp(-Emax/E);
		//if(electron.Weight>0)
		for(int i=0; i<nParticlesPerBin; i++)
			particles.AddPrimary(electron);
	//}

	pe.Run();
}

void IcsCEL::UnitTest()
{
	double Emin=1;
	double step = pow(10.,0.1);
	for(double E = Emin; E<10*Emin; E*=step)
	{
		double rMaxGamma = Emin/E;
		SigmaE sigma(rMaxGamma);
		std::cout << "# Emin = " << Emin << " ; E = " << E << "\n";
		sigma.Print(std::cout, 1000, true, 0, 1e5*E);
		std::cout << "\n\n";
	}
}

IcsCEL::IcsCEL(BackgroundIntegral* aBackground, double aMaxGammaE):
		fBackground(aBackground),
		fMaxGammaE(aMaxGammaE)
{
}

IcsCEL::~IcsCEL()
{

}

CELInteraction* IcsCEL::Clone() const
{
	return new IcsCEL(fBackground->Clone(), fMaxGammaE);
}

RandomInteraction* ICS::Clone() const
{
	return new ICS(fBackground->Clone(), fGammaEmin);
}

IcsCEL::SigmaE::SigmaE(double rGammaMax)
: f_rGammaMax(rGammaMax)
{
	f_m2 = Particle::Mass(Electron);
	f_m2*=f_m2;
}

double IcsCEL::SigmaE::f(double sDim) const
{
	double result = 0;
	double s=sDim/f_m2;//s in units of m^2
	double s_1 = s-1.;
	if(s_1<=0)
		return 0.;
	double a = 1.-1./s;
	if(a>f_rGammaMax)
	{
		a = f_rGammaMax;
		double a2=a*a;
		double s_1_2 = s_1*s_1;
// Expression below is obtained from output of Mathematica:
// F[s_, r_] := (r + 1/r + 4/(s - 1)*(1 - 1/r) + 4/(s - 1)^2*(1 - 1/r)^2)
// Integrate[F[s, r]*(1 - r), {r, 1 - a, 1}, Assumptions -> a < 1 && a > 0]
// F[s,r] is expression inside brackets in formula (23) of astro-ph/9604098
		if(a<1e-3)
		{
			double s2=s*s;
			result = a2*(1.+a/s_1*(-4./3.+a/s_1*(0.25*(9.-6.*s+s2)+0.2*a*(13.-6.*s+s2)+a2/6.*(17.-6.*s+s2))));
		}
		else
		{
			result = -((a*(2.*a2*a*s_1_2 - 6.*(-7. + s)*(1 + s) +
				3.*a*(-5. + s*(-10. + 3.*s)) - a2*(5. + s*(2. + 5.*s))) +
				6.*(-1 + a)*(-7. + s)*(1 + s)*log(1. - a))/(6.*(-1 + a)*s_1_2));
		}
		ASSERT_VALID_NO(result);
	}
	else
	{
		// Expression below is obtained from output of Mathematica:
		// F[s_, r_] := (r + 1/r + 4/(s - 1)*(1 - 1/r) + 4/(s - 1)^2*(1 - 1/r)^2)
		// Integrate[F[s, r]*(1 - r), {r, 1/s, 1}, Assumptions -> s > 1]
		// F[s,r] is expression inside brackets in formula (23) of astro-ph/9604098
		if(s_1<1e-3)
		{
			result = s_1*s_1*(2./3. + s_1*(-7./5.+s_1*(49./20.-404./105.*s_1)));
		}
		else
		{
			result = -1.*(((s_1*(2. + s*(-5. + s*(-3. + s*(-71. + 5.*s)))))/(6.*s*s*s) -
			       (-7 + s)*(1 + s)*log(s))/s_1/s_1);
		}
		ASSERT_VALID_NO(result);
	}
	result *= (SigmaCoef/(sDim-f_m2));
	ASSERT_VALID_NO(result);
	return result;
}

double IcsCEL::SigmaE::Xmin() const
{
	return f_m2;
}

double IcsCEL::Rate(const Particle& aParticle) const
{
	if(aParticle.Type != Electron && aParticle.Type != Positron)
		return 0.;

	double rMaxGamma = fMaxGammaE/aParticle.Energy;
	SigmaE sigma(rMaxGamma);
	return fBackground->GetRateS(sigma, aParticle);
}

} /* namespace Interactions */

