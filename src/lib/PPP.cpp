/*
 * PPP.h
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


#include "PPP.h"
#include "MathUtils.h"
#include "Output.h"
#include "ParticleStack.h"
#include "TableBackgrounds.h"
#include "PropagationEngine.h"
#include "ProtonPP.h"
#include "Test.h"

namespace Interactions {
using namespace mcray;
using namespace Utils;

PPP::PPP(BackgroundIntegral* aBackground, double aScaleFactor):
fBackground(aBackground),
f_ScaleFactor(aScaleFactor),
fMp(Particle::Mass(Proton)),
fMe(Particle::Mass(Electron))
{
	ASSERT(aScaleFactor>0);
}

RandomInteraction* PPP::Clone() const
{
	return new PPP(fBackground->Clone(), f_ScaleFactor);
}

PPP::~PPP() {
}

PPP::SigmaTot::SigmaTot():
fMp(Particle::Mass(Proton)),
fMe(Particle::Mass(Electron)),
fThresholdS(fMp*(4.*fMe+fMp)),
fMaxS(1e300),//1e2*fMp*fMp
fSigmaConst(AlphaEM*AlphaEM*AlphaEM/fMe/fMe)//alpha*r_0
{

}

double PPP::SigmaTot::f(double s) const
{
	double k = 0.5*(s/fMp-fMp)/fMe;//photon energy in proton rest frame in units of electron mass
	if(k<=2)
		return 0;
	double result = 0.;
	if(k<4)
	{//formula (A1) from M. J. Chodorowski, A. A. Zdziarski, and M. Sikora, Astrophys. J. 400, 181 (1992)
		double eta = (k-2.)/(k+2.);
		result = (k-2.)/k;
		result *= (result*result*2./3.*M_PI*(1.+eta*(0.5+eta*(23./40.+eta*(37./120.+eta*61./192.)))));
	}
	else
	{//formula (A2) from M. J. Chodorowski, A. A. Zdziarski, and M. Sikora, Astrophys. J. 400, 181 (1992)
		double ln2k = log(2.*k);
		double _kk = 4./k/k;
		double _kk2 = _kk*_kk;
		result = 28./9.*ln2k-218./27.+_kk*(ln2k*(6.-M_PI*M_PI/3.+ln2k*(-1+2./3.*ln2k))+0.549054066848226);
		result -= (_kk2*(0.1875*ln2k+0.125+_kk*(0.0125868055555556*ln2k-0.00557002314814815)));
	}
	result *= fSigmaConst;
	ASSERT_VALID_NO(result);
	return result;
}

PPP::DifSigma_Eprime::DifSigma_Eprime():
fE(-1),
fK(-1),
fCoef(-1)
{
	fKern_Theta.function = sigma;
	fKern_Theta.params = this;
}

double PPP::s_epsRel = 1e-4;

double PPP::DifSigma_Eprime::f(double E) const
{
	((PPP::DifSigma_Eprime*)this)->fE = E;
	double result = ((Utils::MathUtils)fMath).Integration_qag(fKern_Theta,-1.,1.,0, s_epsRel, (size_t)(1./s_epsRel));
	ASSERT_VALID_NO(result);
	return result*fCoef;
}

void PPP::DifSigma_Eprime::SetK(double aK)
{
	fK = aK;
	fCoef = 0.5*fSigmaConst/(aK*aK*aK);
}


double PPP::DifSigma_Eprime::sigma(double aCosTheta, void* data)
{
	const PPP::DifSigma_Eprime* sig = (const PPP::DifSigma_Eprime*)data;
	return PPP::DiffSigma(sig->fK, sig->fE, aCosTheta);
}

//d(sigma)/d(e)/d(eps)
double PPP::DiffSigma(double k, double em, double cosT)
{//formula (10) of Blumenthal, Phys Rev D (1970) vol 1, number 6, page 1596
// all quantities are given in units of electron mass
// em - electron energy in proton rest frame
// k - photon energy in proton rest frame
// cosT - angle between electron and photon in proton rest frame
// constant multiplier alpha*r_0^2/2k^3 is omitted
	if(em<=1.)
		return 0.;
	double em2 = em*em;
	double ep = k-em;//positron energy in proton rest frame (we neglect proton recoiling momentum)
	if(ep<=1.)
		return 0.;
	double ep2 = ep*ep;
	double pm2 = em2-1.;
	double pm = sqrt(pm2);//electron momentum in proton rest frame
	double pp = sqrt(ep2-1.);//positron momentum in proton rest frame
	double e=em-cosT*pm;//delta_ - electron energy in the LAB (related to observer on Earth) frame divided by proton gamma factor
	double e2 = e*e;
	double e4 = e2*e2;
	double k2 = k*k;
	double sin2T = 1.-cosT*cosT;
	double T2 = k2+pm2-2.*k*pm*cosT;
	double T = sqrt(T2);
	double Y = 2./pm2*log((ep*em+pp*pm+1)/k);
	double ep_pp = (ep-pp<1e-8)? 0.5/ep/ep : (ep-pp);
	double yp =	log((ep+pp)/ep_pp)/pp;
	double deltaTp = log((T+pp)/(T-pp));
	double result = -4.*sin2T*(2.*em2+1.)/(pm2*e4) + (5.*em2-2.*ep*em+3.)/(pm2*e2) + (pm2-k2)/(T2*e2) + 2.*ep/(pm2*e);
	result += (Y/(pm*pp)*(2.*em*sin2T*(3.*k+pm2*ep)/e4 + (2.*em2*(em2+ep2)-7.*em2-3.*ep*em-ep2+1.)/e2 + k*(em2-em*ep-1.)/e));
	result -= (deltaTp/(pp*T)*(2./e2-3.*k/e-k*(pm2-k2)/(T2*e)) + 2.*yp/e);
	result *= (pp*pm);
	if(result<0)
		return 0.;
	ASSERT_VALID_NO(result);
	return result;
}

//d(sigma)/d(e)/d(eps)
double PPP::DifSigma_r::Kern::f(double eps) const
{//formula (10) of Blumenthal, Phys Rev D (1970) vol 1, number 6, page 1596
 // all quantities are given in units of electron mass
 // eps = E_ - electron energy in proton rest frame
 // k - photon energy in proton rest frame
 // e=delta_ - electron energy in the LAB (related to observer on Earth) frame divided by proton gamma factor
 // cos(Theta)=(E_-e)/p_; constant multiplier alpha*r_0^2/2k^3 is omitted
 // and extra multiplier dcosTheta/de = 1/p_ is introduced to take into account change of integration variables cosTheta -> e
	double em = eps;//electron energy in proton rest frame
	ASSERT(em>1);
	double em2 = em*em;
	double ep = k-em;//positron energy in proton rest frame (we neglect proton recoiling momentum)
	ASSERT(ep>1);
	double ep2 = ep*ep;
	double pm2 = em2-1.;
	double pm = sqrt(pm2);//electron momentum in proton rest frame
	double pp = sqrt(ep2-1.);//positron momentum in proton rest frame
	double cosT = (em-e)/pm; //cos(Theta) the angle between photon and electron in p-rest frame
	double sin2T = 1.-cosT*cosT;
	if(cosT>=1.)
	{
		sin2T = (2.*e*pm-1.)/pm2;
		if(sin2T<0)
			return 0;
		ASSERT(sin2T<1e-3);
		cosT = 1.-0.5*sin2T;
	}
	double T2 = k2+pm2-2.*k*pm*cosT;
	double T = sqrt(T2);
	double Y = 2./pm2*log((ep*em+pp*pm+1)/k);
	double ep_pp = (ep-pp<1e-8)? 0.5/ep/ep : (ep-pp);
	double yp =	log((ep+pp)/ep_pp)/pp;
	double deltaTp = log((T+pp)/(T-pp));
	double result = -4.*sin2T*(2.*em2+1.)/(pm2*e4) + (5.*em2-2.*ep*em+3.)/(pm2*e2) + (pm2-k2)/(T2*e2) + 2.*ep/(pm2*e);
	result += (Y/(pm*pp)*(2.*em*sin2T*(3.*k+pm2*ep)/e4 + (2.*em2*(em2+ep2)-7.*em2-3.*ep*em-ep2+1.)/e2 + k*(em2-em*ep-1.)/e));
	result -= (deltaTp/(pp*T)*(2./e2-3.*k/e-k*(pm2-k2)/(T2*e)) + 2.*yp/e);
	result *= pp;
	if(result<0)
		return 0.;//TODO find the reason
	//ASSERT_VALID_NO(result);
	return result;
}

void PPP::DifSigma_r::Kern::SetParams(double _k, double _e)
{
	k = _k;
	k2 = k*k;
	e = _e;//electron energy in lab frame divided by proton gamma factor
	e2 = e*e;
	e4 = e2*e2;
	fEpsMin = 0.5*(e + 1./e);
	//ASSERT(Xmin()<=Xmax());//if Xmin>Xmax PPP::DifSigma::f(double r) will return 0
	//(this may occure for large k (k>m_p, since approximation of zero proton kinetic energy is not valid anymore)
}

void PPP::DifSigma_r::SetS(double aS)
{
	fK = 0.5*(aS/fMp-fMp)/fMe;//photon energy in proton rest frame in units of electron mass
	double m = fMp/fMe;//proton mass in units of electron mass
	double tmp = fK + m - 1.;
	tmp *= tmp;
	//maximal momentum of electron directed opposite to photon momentum
	//in proton rest frame in units of electron mass
	double maxP = (sqrt(fK*(fK*(-1 + m) - 2*m)*(-1 + m)*tmp) +
				   fK*(-1 + fK + m - fK*m))/((-1 + m)*(-1 + 2*fK + m));
	f_rMax = (sqrt(maxP*maxP+1.)+maxP)/m;//converting to lab frame and dividing by proton energy
	ASSERT(f_rMax<1);
	//maximal momentum of electron directed along photon momentum
	//in proton rest frame in units of electron mass
	double minP = ((-1 + fK)*fK*(-1 + m) + sqrt(fK*(fK*(-1 + m) - 2*m)*(-1 + m)*
												tmp))/((-1 + m)*(-1 + 2*fK + m));

	//converting to lab frame and dividing by proton energy
	f_rMin = (minP<1e4 ? (sqrt(minP*minP+1.)-minP) : 0.5/minP/minP)/m;
}

double PPP::DifSigma_r::f(double r) const
{
	double e = fMp/fMe*r;//  E_e/gamma_p in units of electron mass
	Kern fKern;
	fKern.SetParams(fK, e);
	double epsMin = fKern.Xmin();
	double epsMax = fKern.Xmax();
	double relDiff = (epsMax-epsMin)/epsMax;
	if(relDiff<1e-9)
		return 0.;
	const size_t limit = (size_t)(1./s_epsRel);
	double relError = ( relDiff > 1e-5)? s_epsRel : (relDiff<3e-7?0.1:0.01);
	double result = ((MathUtils)fMath).Integration_qag((gsl_function)fKern, epsMin, epsMax, 0, relError, limit);
	ASSERT_VALID_NO(result);
	result *= (0.5*fSigmaConst*fMp/fMe/fK/fK/fK);
	return result;
}

double PPP::Rate(const Particle& aParticle) const
{
	if(aParticle.Type != Proton)
		return 0.;
	return fBackground->GetRateS(fSigma, aParticle) / f_ScaleFactor;
}

bool PPP::SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const
{
	ASSERT(aParticle.Type == Proton);
	aS = 0.;
	return (fBackground->GetRateAndSampleS(fSigma, aParticle, aRandomizer, aS)!=0);
}

void PPP::SampleSecondariesMean(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const
{//// sampling secondary electrons/positrons is expensive so thinning is performed before sampling
 /// via usage of scale factor parameter and generating either electron or positron (not both)
	ASSERT(aParticle.Type == Proton);

	Particle sec = aParticle;
	sec.Type = aRandomizer.Rand()>0.5 ? Electron : Positron;
	DifSigma_r difSigma;
	double r = difSigma.SampleR(aS, s_epsRel, math, aRandomizer.Rand());
	sec.Energy *= r;
	sec.Weight *= (f_ScaleFactor*2);//2 takes into account that pair is produced
	sec.fCascadeProductionTime = aParticle.Time;
	aSecondaries.push_back(sec);

	//don't add proton to the list of secondaries since this.KillsPrimary == false
}

	double PPP::DifSigma_r::SampleR(double aS, double aEpsRel, const Utils::MathUtils& aMath, double aRand)
	{
		double r = 0.;
		SetS(aS);
		if(fK<20000)
		{
			double sigmaTot;
			aMath.SampleLogDistribution(*this, aRand, r, sigmaTot, Xmin(), Xmax(), aEpsRel);
			ASSERT(r>=f_rMin && r<=f_rMax);
		}
		else
		{ //to avoid rounding error using analytical estimate dSigma/dr ~ r^{-7/4} for fK>>1
			// < r > = m_e/m_p //this is clear from plotting r^2 dSigma/dr in log scale (try e.g. PPP::UnitTest)
			double rMinAdj = f_rMin;
			for(int nIt=1; true; nIt++)//approximately solving eq.
			{
				double rMinNew = pow(3.*fMp/fMe*(pow(f_rMax,0.25)-pow(rMinAdj,0.25))+pow(f_rMax,-0.75),-4./3.);
				ASSERT(rMinNew>0 && rMinNew<f_rMax);
				double relErr = 2.*fabs(rMinNew-rMinAdj)/(rMinNew+rMinAdj);
				rMinAdj = rMinNew;
				if(relErr<aEpsRel)
					break;
				ASSERT(nIt<50);
			}
			//double rMaxAdj = fMe/fMp/3./pow(f_rMin,0.75) + pow(f_rMin,0.25);   //adjust rMin to have < r > = m_e/m_p (assuming rMin<<rMax)
			//ASSERT(rMaxAdj<1.);
			//ASSERT(rMaxAdj/f_rMin<1e-3);
			r = pow(pow(rMinAdj,-0.75)*(1.-aRand)+aRand*pow(f_rMax,-0.75) , -4./3.);
			ASSERT(r>=rMinAdj && r<=f_rMax);
		}
		return r;
	}


void PPP::SampleSecondariesExact(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const
{
	ASSERT(aParticle.Type == Proton);
	double m = aParticle.Mass();
	double k = 0.5*(aS/m-m)/fMe;//photon energy in proton rest frame in units of electron mass
	ASSERT(k>2);
	((PPP::DifSigma_Eprime&)fDifSigma_Eprime).SetK(k);

	double Eprime, cosTheta, sigmaTot;
	((MathUtils&)math).SampleLogDistribution(fDifSigma_Eprime, aRandomizer.Rand(), Eprime, sigmaTot, fDifSigma_Eprime.Xmin(), fDifSigma_Eprime.Xmax(), 10.*s_epsRel);
	CosThetaSamplingEq cosProb(Eprime, k);
	((MathUtils&)math).SampleDistribution(cosProb, aRandomizer.Rand(), cosTheta, sigmaTot, cosProb.Xmin(), cosProb.Xmax(), s_epsRel);
	Particle electron = aParticle;
	electron.Type = Electron;
	electron.Energy *= ((Eprime-cosTheta*sqrt(Eprime-1.))*fMe/m);
	Particle positron = aParticle;
	positron.Type = Positron;
	Eprime = k-Eprime;//we neglect proton kinetic energy
	positron.Energy *= ((Eprime+cosTheta*sqrt(Eprime-1.))*fMe/m);//we assume that positron is directed against electron??? TODO: this is probably wrong
	aParticle.Energy -= (electron.Energy + positron.Energy);
	if(aParticle.Energy < m)
		Exception::Throw("PPP::SampleSecondariesExact error: Ep<m");
	aParticle.fCascadeProductionTime = aParticle.Time;
	electron.Weight *= f_ScaleFactor;
	positron.Weight *= f_ScaleFactor;
	aSecondaries.push_back(electron);
	aSecondaries.push_back(positron);
}

void PPP::SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const
{
	//SampleSecondariesExact(aParticle, aSecondaries, aS, aRandomizer); //wrong calculation (inelasticity growing with S instead of dropping)
	//also in SampleSecondariesExact we assume that positron is directed against electron which is probably wrong

	//this should to be correct on the average although not precise for individual interactions
	//since electron and positron energies are sampled independently of each other for individual interactions
	SampleSecondariesMean(aParticle, aSecondaries, aS, aRandomizer);
}

double PPP::CosThetaSamplingEq::f(double aCosT) const
{
	return PPP::DiffSigma(fK,fE,aCosT);
}

void PPP::UnitTest()
{
	Utils::MathUtils math;
	SigmaTot sigma;
	DifSigma_r difSigma;
	double Mp(Particle::Mass(Proton));
	double Me(Particle::Mass(Electron));
	std::ofstream sigmaPPout;
	sigmaPPout.open("ppp_diff_sigma",std::ios::out);
	std::cout << "k/m_e\tsigTot/mbarn\t2.*(sigTot-sigIntTot)/(sigTot+sigIntTot)" << std::endl;
	for(double k=2.2; k<3e4; k*=10.)
	{
		double s = Mp + k*Me;
		s*=s;
		//if(s>sigma.Xmax())
		//	continue;
		double sigTot = sigma(s);
		difSigma.SetS(s);
		double sigIntTot = math.Integration_qag((gsl_function)difSigma, difSigma.Xmin(), difSigma.Xmax(), 0, 1e-5, 10000);
		std::cout << k << "\t" << sigTot/units.barn*1e3 << "\t" << 2.*(sigTot-sigIntTot)/(sigTot+sigIntTot) << std::endl;
		sigmaPPout << "# k=" << k << std::endl;
		difSigma.Print(sigmaPPout,100,true);
		sigmaPPout << "\n" << std::endl;
	}
	UnitTestSampling(1e12, 0);
	UnitTestEconserv(1e20, 0.1);
}

void PPP::UnitTestSampling(double aE, double aZ)
{
	std::string out = "debug_sample_ppp_E";
	out += ToString(aE*1e6);//eV
	out += "_Z";
	out += ToString(aZ);
	std::ofstream secPPout;
	secPPout.open("sample_sec_ppp",std::ios::out);
	//std::ofstream secPPoutVarS;
	//secPPoutVarS.open("sample_sec_ppp_varS",std::ios::out);

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
	Interactions::PPP pp(backgrI);
	Particle proton(Proton, aZ);
	proton.Energy = aE;//MeV
	proton.Time.setZ(aZ);
	proton.Type = Proton;
	Randomizer rnd;
	/*
	double S=2.;
	for(int i=0; i<10000; i++)
	{
		std::vector<Particle> secondaries;
		if(pp.GetSecondaries(proton, secondaries, rnd))//this outputs sampled s to debug file
		{		std::vector<Particle> secondariesS;
		pp.SampleSecondaries(proton, secondariesS, S rnd);
			ASSERT(secondaries.size()==3);
			//print electron and proton energies
			secPPoutVarS << secondaries[0].Energy/aE << "\n" << secondaries[1].Energy/aE << "\n";
		}
	}*/

	double Mp(Particle::Mass(Proton));
	double Me(Particle::Mass(Electron));
	SigmaTot sigma;
	for(double k=2.2; k<3e5; k*=10.)
	{
		double s = Mp + k*Me;
		s*=s;
		if(s>sigma.Xmax())
			continue;
		for(int i=0; i<100; i++)
		{
			std::vector<Particle> secondaries;
			pp.SampleSecondaries(proton, secondaries, s, rnd);
			ASSERT(secondaries.size()==1);//only electron or positron is sampled
			secPPout << k << "\t" << secondaries[0].Energy/aE << "\n";
			//std::cout << k << "\t" << secondaries[0].Energy/aE << "\n";
		}
		secPPout << "\n" << std::endl;
		//std::cout << std::endl;
	}

	secPPout.close();
}

void PPP::UnitTestEconserv(double aE, double aZ){
	uint scaleFactorPPP = 10;
	int nParticles = 100;
	Test::EConservationTest test(aZ,1e6,2015);
	test.alphaThinning = 0.9;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	test.Engine().AddInteraction(new PPP(test.Backgr(),scaleFactorPPP));
	test.Engine().AddInteraction(new ProtonPPcel(test.Backgr()));
	Particle p(Proton, aZ);
	p.Energy = aE*units.eV;
	test.SetPrimary(p, nParticles);
	test.Run();
}

} /* namespace Interactions */
