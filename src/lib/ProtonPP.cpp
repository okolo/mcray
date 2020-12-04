/*
 * ProtonPP.cpp
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

#include "ProtonPP.h"
#include <fstream>

//used for unit test only
#include "PropagationEngine.h"

namespace Interactions {
using namespace mcray;

ProtonPPSigma::ProtonPPSigma()
{
	double Me = Particle::Mass(Electron);
	double Mp = Particle::Mass(Proton);
	fMpMe_1 = 1./(Me*Mp);
	fThresholdS=2.*Me + Mp;
	fThresholdS *= fThresholdS;
	fMp2 = Mp*Mp;
	fCoef = AlphaEM*AlphaEM*AlphaEM/Me/Me; //alpha*Z^2*r_0^2 == alpha^3*Z^2/m_e^2
}

double ProtonPPSigma::sigmaLab(double k) const
{
	double result = 0;
	if(k<4)//using formula (A1) from M. J. Chodorowski, A. A. Zdziarski, and M. Sikora, Astrophys. J. 400, 181 (1992)
	{
		double eta = (k - 2.)/(k+2);
		double eta2 = eta*eta;
		double kk = (k - 2.)/k;
		//2pi/3*fCoef * kk^3 * (1+1/2 eta + 23/40 eta^2 + 37/120 eta^3 + 61/192 eta^4)
		result = 2.0943951023932*fCoef*kk*kk*kk*(1.+0.5*eta+0.575*eta2+
				0.308333333333333*eta2*eta+0.317708333333333*eta2*eta2);
	}
	else
	{//using formula (A2) from M. J. Chodorowski, A. A. Zdziarski, and M. Sikora, Astrophys. J. 400, 181 (1992)
		double ln2k = log(2.*k);
		double ln2k2 = ln2k*ln2k;
		double coef2 = 6.*ln2k + 0.666666666666667*ln2k2*ln2k-ln2k2-3.28986813369645*ln2k + 0.549047873167415;
		double coef4 = -(0.1875*ln2k+0.125);
		double coef6 = -(0.0125868055555556*ln2k-0.00557002314814815);
		double r2 = 4./k/k;
		double r4 = r2*r2;
		result = fCoef*(3.11111111111111*ln2k-8.07407407407407+r2*coef2+r4*coef4+r2*r4*coef6);
	}
	ASSERT_VALID_NO(result);
	return result;
}

double ProtonPPSigma::f(double s) const
{
	if(s<=fThresholdS)
		return 0.;

	double k = 0.5*(s - fMp2)*fMpMe_1; // photon energy in proton rest frame in units of electron mass
	return sigmaLab(k);
}

ProtonPP::ProtonPP(BackgroundIntegral* aBackground):
fBackground(aBackground)
{
}

ProtonPP::~ProtonPP() {
}

double ProtonPP::Rate(const Particle& aParticle) const
{
	if(aParticle.Type != Proton)
		return 0.;
	return fBackground->GetRateS(fSigma, aParticle);
}

bool ProtonPP::SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const
{
	ASSERT(aParticle.Type == Proton);
	aS = 0.;
	return (fBackground->GetRateAndSampleS(fSigma, aParticle, aRandomizer, aS)!=0);
}

double ProtonPP::maxPPPemaxError = 0.;

void ProtonPP::SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const
{
	ASSERT(aParticle.Type == Proton);
	Particle secProton = aParticle;
	Particle secElectron = aParticle;
	Particle secPositron = aParticle;
	secElectron.Type = Electron;
	secPositron.Type = Positron;
	double E = aParticle.Energy;
	double Mp = Particle::Mass(Proton);
	double k = 0.5*(aS/Mp-Mp);
	double Me = Particle::Mass(Electron);
	double betaP = aParticle.beta();
	double gammaP = E / Mp;
	double mFrac = Mp/Me;
	double km = k/Me;
	double c1 = mFrac*km-2.;

	///formula (3.28) J. W. Motz, H. A. Olsen, and H. W. Koch, Rev. Mod. Phys. 41, 581 (1969).
	double pMaxProtonLab = (km*(mFrac*(km+mFrac)-2)+(km+mFrac)*sqrt(c1*c1-4.*mFrac*mFrac))/mFrac/(2.*km+mFrac)*Me;
	ASSERT_VALID_NO(pMaxProtonLab);
	double EmaxProtonLab = sqrt(pMaxProtonLab*pMaxProtonLab+Mp*Mp);
	double EminProton = gammaP*(EmaxProtonLab - betaP*pMaxProtonLab);

	///formula (3.43a) J. W. Motz, H. A. Olsen, and H. W. Koch, Rev. Mod. Phys. 41, 581 (1969).
	double c2 = km+mFrac;
	double c3 = mFrac*c2*(km-1.);
	double c4 = km*sqrt(km*km + mFrac*mFrac*km*(km-2.));
	double c5 = (c2*c2-km*km);

	double EeLabMaxAntiparallelToK = (c3-c4)/c5*Me;
	double EeLabMaxParallelToK = (c3+c4)/c5*Me;

	double PeLabMaxAntiparallelToK = (EeLabMaxAntiparallelToK>1.001*Me) ? sqrt(EeLabMaxAntiparallelToK*EeLabMaxAntiparallelToK - Me*Me) : 0.;
	double PeLabMaxParallelToK = (EeLabMaxParallelToK>1.001*Me) ? sqrt(EeLabMaxParallelToK*EeLabMaxParallelToK - Me*Me) : 0.;
	double Emin = gammaP*(EeLabMaxParallelToK - betaP*PeLabMaxParallelToK);
	double Emax = gammaP*(EeLabMaxAntiparallelToK + betaP*PeLabMaxAntiparallelToK);
	if(Emin>Emax)
	{//this may occur if EeLabMaxAntiparallelToK << Me, i.e. close to threshold
		ASSERT(EeLabMaxAntiparallelToK<1.001*Me);
		double tmp = Emin;
		Emin = Emax;
		Emax = tmp;
	}

	ASSERT(Emin>Me);
	ASSERT(Emax<E);

	if(E < (EminProton+Emax+Emin))
	{
		//ASSERT(EminProton+Emax+Emin<E);
		double correctedEmax = E - EminProton - Emin;
		double relError = (Emax-correctedEmax)/Emax;
		if(relError>maxPPPemaxError)
		{
			std::cerr << "PPP Emax error " << relError << std::endl;
			maxPPPemaxError = relError;
		}
		Emax = correctedEmax;
	}

	double E1 = SampleE(Emin, Emax, aRandomizer);
	double E2max = E - E1 - EminProton;
	if(E2max>Emax)
		E2max = Emax;
	double E2 = SampleE(Emin, E2max, aRandomizer);
	if(aRandomizer.Rand()>0.5)
	{//Randomly choose which particle has energy E1 and which E2
		double tmp = E1;
		E1 = E2; E2 = tmp;
	}
	secElectron.Energy = E1;
	secPositron.Energy = E2;
	secProton.Energy = E - E1 - E2;
	aSecondaries.push_back(secElectron);
	aSecondaries.push_back(secPositron);
	aSecondaries.push_back(secProton);
}

void ProtonPP::UnitTest()
{
	double Zmax = 1;
	double Emin = 1e6*units.eV;
	double Emax = 1e19*units.eV;
	int Nparticles = 1000;
	std::string outputDir = "testProtonPP";
	unsigned int kAc = 1;
	double epsRel = 1e-3;
	double stepZ = Zmax<0.05 ? Zmax/2 : 0.025;
	double logStepK = pow(10,0.05/kAc);
	if(!cosmology.IsInitialized())
		cosmology.Init(Zmax + 10);
	double cmbTemp = 2.73*units.K;
	double alphaThinning = 0.5;
	CompoundBackground backgr;
	backgr.AddComponent(new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp, 0., Zmax + 1.));
	//IR/O component
	backgr.AddComponent(new GaussianBackground(0.1*units.eV, 0.01*units.eV, 5, 1./units.cm3, 0, Zmax + 1.));
	PBackgroundIntegral backgrI(new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel));
	int seed = 2015;
	Result result(new EnergyBasedFilter(Emin, M_PI), true);

	SmartPtr<RawOutput> pOutput = new RawOutput(outputDir, false);
	result.AddOutput(pOutput);
	ParticleStack particles;
	PropagationEngine pe(particles, result, seed);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);

	Particle proton(Proton, Zmax);
	int nSteps = log(Emax/Emin)/log(logStepK)+1.;
	double mult = pow(Emax/Emin, 1./nSteps);
	{
		std::ofstream rateOut;
		rateOut.open((outputDir + "/PPP").c_str(),std::ios::out);
		PRandomInteraction i = new ProtonPP(backgrI);
		pe.AddInteraction(i);
		int step=0;
		for(proton.Energy=Emax/100; step<=nSteps; proton.Energy*=mult, step++)
			rateOut << proton.Energy/units.eV << "\t" << i->Rate(proton)*units.Mpc << "\n";
		rateOut.close();
	}
	Randomizer rand;
	CosmoTime tEnd;
	proton.Energy = Emax;
	pOutput->SetOutputDir(outputDir + "/z0");
	result.SetEndTime(tEnd);

	for(int i=0; i<Nparticles; i++)
	{
		particles.AddPrimary(proton);
	}
	pe.RunMultithreadReleaseOnly();
}

double ProtonPP::SampleE(double aEmin, double aEmax, Randomizer& aRandomizer) const
{//assuming dif sigma is proportional to E^{-7/4}
	double r = aRandomizer.Rand();
	double E = pow(pow(aEmin,-0.75)*(1.-r) + pow(aEmax,-0.75)*r, -4./3.);
	ASSERT(E>=aEmin && E<=aEmax);
	return E;
}

RandomInteraction* ProtonPP::Clone() const
{
	return new ProtonPP(fBackground->Clone());
}

ProtonPPcel::SigmaE::SigmaE()
{
	fInelCoef = 4.*Particle::Mass(Electron)/Particle::Mass(Proton);
}

double ProtonPPcel::SigmaE::f(double s) const
{
	double k = 0.5*(s - fMp2)*fMpMe_1; // photon energy in p-rest frame in units of electron mass
	if(k<=2.)
		return 0.;

	//using formulas (3.7)-(3.10) from M. J. Chodorowski, A. A. Zdziarski, and M. Sikora, Astrophys. J. 400, 181 (1992)
	double inelasticity = 0.;
	if(k<1000.)
	{
		double ln = log(k-1);
		double lnI = 1.;
		for(int i=0; i<4; i++, lnI *= ln)
			inelasticity += lnI*a[i];
	}
	else
	{
		double ln = log(k);
		double lnI = 1.;
		for(int i=0; i<4; i++, lnI *= ln)
			inelasticity += lnI*b[i];
		inelasticity /= (3.11111111111111*(ln + M_LN2)-8.07407407407407);
	}
	inelasticity *= (fInelCoef/k);
	ASSERT_VALID_NO(inelasticity);
	return ProtonPPSigma::sigmaLab(k)*inelasticity;
}

//constants given by (3.8) from M. J. Chodorowski, A. A. Zdziarski, and M. Sikora, Astrophys. J. 400, 181 (1992)
const double ProtonPPcel::SigmaE::a[] = {1, 0.3958, 0.1, 0.00781};

//constants given by (3.10) from M. J. Chodorowski, A. A. Zdziarski, and M. Sikora, Astrophys. J. 400, 181 (1992)
const double ProtonPPcel::SigmaE::b[] = {-8.778, 5.512, -1.614, 0.666666666666667};

ProtonPPcel::ProtonPPcel(BackgroundIntegral* aBackground):
fBackground(aBackground)
{
}

ProtonPPcel::~ProtonPPcel() {
}

void ProtonPPcel::UnitTest()
{
	std::ofstream rateOut;
	rateOut.open("ppp_test",std::ios::out);
	SigmaE sigmaE;
	ProtonPPSigma sigma;
	double Mp = Particle::Mass(Proton);
	double Me = Particle::Mass(Electron);
	double coef = 3.88593900935001e-07/Mp/Me;
	for(double k = 2.; k<=1e4; k *= 1.58489319246111)
	{//compare sigK to dashed curve from figure 2 of M. J. Chodorowski, A. A. Zdziarski, and M. Sikora, Astrophys. J. 400, 181 (1992)
		double s = Mp*(Mp + 2.*k*Me);
		double sig = sigma(s)/coef;
		double sigK = sigmaE(s)/coef;
		rateOut << k << "\t" << sig << "\t" << sigK << "\n";
	}
	rateOut.close();
}

double ProtonPPcel::Rate(const Particle& aParticle) const
{
	if(aParticle.Type != Proton)
		return 0.;

	return fBackground->GetRateS(fSigma, aParticle);
}

void ProtonPPcel::GetSecondaries(const Particle& aParticle, std::vector<Particle>& aSecondaries,
				cosmo_time aDeltaT, Randomizer& aRandomizer) const
{

}

CELInteraction* ProtonPPcel::Clone() const
{
	return new ProtonPPcel(fBackground->Clone());
}

} /* namespace Interactions */
