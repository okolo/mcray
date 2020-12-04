/*
 * Test.cpp
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

#include "Test.h"
#include <cstdlib>
#include "ParticleStack.h"
#include "Utils.h"
#include "Cosmology.h"
#include "PropagationEngine.h"
#include <iostream>
#include "Output.h"
#include "TableBackgrounds.h"
#include "ICS.h"
#include "GammaPP.h"
#include "ProtonPP.h"
#include "PPP.h"
#include "GZK.h"
#include "Deflection3D.h"
#include "NeutronDecay.h"
#include "SteckerEBL.h"
#include "PhotoDisintegration.h"
#include "Stecker16Background.h"

#ifdef ELMAG_TEST
#include "ElmagTest.h"
#endif

using namespace std;
using namespace mcray;
using namespace Utils;
using namespace Backgrounds;
using namespace Interactions;

namespace Test {

void BuildInitialState(ParticleStack& particles, double aZmax)
{
	Particle gamma(Photon, aZmax);
	gamma.Energy = 1e7/units.Eunit;//MeV

	particles.AddPrimary(gamma);
}

void TestList::EngineTest(double aAlphaThinning)
{
	class SimpleSplitting : public RandomInteraction
	{
	public:
		SimpleSplitting(double aRate):fRate(aRate){}
		RandomInteraction* Clone() const { return new SimpleSplitting(fRate); }
		virtual ~SimpleSplitting(){};
		double Rate(const Particle& aParticle) const {return fRate;}
		bool GetSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const
		{
			Particle product = aParticle;
			product.Energy = 0.5 * product.Energy;
			aSecondaries.push_back(product);
			aSecondaries.push_back(product);
			return true;
		}
	private:
		double fRate;
	};
	double Zmax = 0.01;
	double Emin = 1.;
	double Emax = 1000.;
	int Ninteractions = 100;

	ParticleStack particles;
	Result result(Emin);
	string outFile = "engine_alpha_" + ToString(aAlphaThinning);
	result.AddOutput(new SpectrumOutput(outFile, Emin, pow(10, 0.05)));

	cosmology.Init();

	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(aAlphaThinning);
	pe.SetThinning(&thinning);
	pe.AddInteraction(new SimpleSplitting(Ninteractions/cosmology.z2t(Zmax)));
	Randomizer r;

	/// build initial state
	double t1 = cosmology.z2t(Zmax);
	double t2 = cosmology.z2t(0.);
	int Nparticles = 10000000;
	double mult = pow(10, 0.05);
	int nBins = (int)(log(Emax/Emin)/log(mult) + 0.5);
	int nParticlesPerBin = Nparticles/nBins;
	double E=Emin*mult;
	for(int i=0; i<nBins; i++, E*=mult)
	{
		for(int j=0; j<nParticlesPerBin; j++)
		{
			double t = t1+r.Rand()*(t2-t1);
			Particle p(Photon, cosmology.t2z(t));
			p.Energy = E;
			p.Weight = 1/E;///dN/dlogE
			particles.AddPrimary(p);
		}
	}
	pe.Run();
}

void TestList::FrameworkTest()
{
	double Zmax = 2;
	double Emin = 1./units.Eunit;
	//double alphaThinning = 0;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	result.AddOutput(new SpectrumOutput("out", Emin, pow(10, 0.05)));

	cosmology.Init();

	CompoundBackground backgr;
	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;

	backgr.AddComponent(new PlankBackground(cmbTemp, 1e-5*cmbTemp, 1e3*cmbTemp));
	//backgr.AddComponent(new Backgrounds::Kneiske1001EBL());
	//double stepZ = 0.1;
	//double logStepK = pow(10,0.05);
	//double epsRel = 1e-3;

	//BackgroundIntegral backgrI(backgr, stepZ, logStepK, Zmax, epsRel);

	PropagationEngine pe(particles, result, 2011);
	//EnergyBasedThinning thinning(alphaThinning);
	//pe.SetThinning(&thinning);
	pe.AddInteraction(new RedShift());
	//pe.AddInteraction(new Interactions::ICS(backgrI,Emin));
	//pe.AddInteraction(new Interactions::IcsCEL(backgrI,Emin));
	//pe.AddInteraction(new Interactions::GammaPP(backgrI));

	/// build initial state
	BuildInitialState(particles, Zmax);

	pe.Run();
}

int TestList::Main(int argc, char** argv)
{
	if(argc==0)
	{
		FrameworkTest();
		return 0;
	}
	const char* command = argv[0];
	std::map<std::string,TestProc>::iterator tit = fTests.find(command);
	if(tit!=fTests.end()){
		TestProc proc = tit->second;
		return proc();
	}



	if(strcmp(command,"elmag")==0)
	{
#ifdef ELMAG_TEST
		ElmagTests::Main(argc-1, argv+1);
#else
		Exception::Throw("Elmag tests are not linked in this version");
#endif
	}
	else if(strcmp(command,"FERMI")==0)
		FermiTests::Main(argc-1, argv+1);
	else if(strcmp(command,"cosmology")==0)
		Cosmology::UnitTest();
	else if(strcmp(command,"background")==0)
		BackgroundTest();
	else if(strcmp(command,"backgroundI")==0)
		ContinuousBackgroundIntegral::UnitTest();
	else if(strcmp(command,"GZK")==0)
		//GZK::UnitTest();
		GZK::UnitTestEconserv(1e20, 0.1);
	else if(strcmp(command,"ICS")==0)
	{
		//ICS::UnitTest(false, true, true,1./units.Eunit, 1e6/units.Eunit, 1e-6);
		//ICS::UnitTest(false, true, false,1./units.Eunit, 1e6/units.Eunit, 1e-6);
		//ICS::UnitTest(false, true, true,100./units.Eunit, 1e6/units.Eunit, 1e-6);
		//ICS::UnitTest(false, true, false,100./units.Eunit, 1e6/units.Eunit, 1e-6);

		ICS::UnitTest(true, false, false,1./units.Eunit, 1e14/units.Eunit, 1e-6);//print rates only
		ICS::UnitTest(true, false, true,1./units.Eunit, 1e14/units.Eunit, 1e-6);//print rates only

		//ICS::UnitTest(true, false, false);
		//ICS::UnitTest(true, false, true);

		//ICS::UnitTestMono(true, 1, true);
		//ICS::UnitTestMono(false, 1, true);
		//ICS::UnitTestMono(false, 10, true);
		//ICS::UnitTestMono(false, 100, true);
		//ICS::UnitTestMono(false, 1000, true);
		//ICS::UnitTest();
	}
	else if(strcmp(command,"ICSsample")==0)
	{
		ICS::UnitTestSampling(1./units.Eunit, 1e8/units.Eunit, 1e-4);
	}
	else if(strcmp(command,"ICSdistr")==0)
	{
		ICS::UnitTestDistribution(1./units.Eunit, 1e8/units.Eunit, 1e-4);
	}
	else if(strcmp(command,"ICScel")==0)
		IcsCEL::UnitTest();
	else if(strcmp(command,"PP")==0)
	{
		//GammaPP::UnitTest();
		//GammaPP::UnitTest(true, true, true);
		//GammaPP::UnitTest(true, true, false);
		GammaPP::UnitTest(true, false, false, false);//rates only
		GammaPP::UnitTest(true, false, false, true);//rates only
	}
	else if(strcmp(command,"PPsample")==0)
	{
		GammaPP::UnitTestSampling(1e8/units.Eunit, 1e-4);
	}
	else if(strcmp(command,"PPP")==0)
	{
		PPP::UnitTest();
		//ProtonPP::UnitTest();
	}
	else if(strcmp(command,"split")==0)
	{
		SplittedBackgroundTest(true);
		SplittedBackgroundTest(false);
	}
	else if(strcmp(command,"thinning")==0)
	{
		ThinningTest(1);
		ThinningTest(0.5);
		ThinningTest(0);
	}
	else if(strcmp(command,"Kachelries")==0)
	{// KachelriesTest(double aEmax_eV, double aZmax, double aAlphaThinning, bool aSplitted, int aN, bool aNoPP, bool aNoICS, ParticleType aPrimary);
		if(argc!=12 && argc!=13)//E${E}_z${z}_N${N}_prim${PrimaryPart}_thin${alpha}
		{
			Exception::Throw("usage: --test Kachelries Emin/eV E/eV z N primary alpha_thinning NoPP NoICSdiscr NoICScel NoRedshift split [output/dir]\nprimary: 0 - photon; -1 - e-; 1 - e+");
		}
		double E,z,alpha,Emin;
		int N,NoPP,NoICSdiscr,NoICScel,NoRedshift,particle,split;
		int pos = 1;
		if(sscanf(argv[pos++],"%lg", &Emin)!=1)
			Exception::Throw("Failed to read Emin parameter");
		if(sscanf(argv[pos++],"%lg", &E)!=1)
			Exception::Throw("Failed to read E parameter");
		if(sscanf(argv[pos++],"%lg", &z)!=1)
			Exception::Throw("Failed to read z parameter");
		if(sscanf(argv[pos++],"%d", &N)!=1)
			Exception::Throw("Failed to read N parameter");
		if(sscanf(argv[pos++],"%d", &particle)!=1)
			Exception::Throw("Failed to read primary parameter");
		ParticleType pt = ParticleTypeEOF;
		switch(particle)
		{
		case 0:	pt = Photon; break;
		case 1:	pt = Positron; break;
		case -1: pt = Electron; break;
		default:
			Exception::Throw("Invalid 'primary' parameter");
			break;
		}
		if(sscanf(argv[pos++],"%lg", &alpha)!=1)
			Exception::Throw("Failed to read alpha_thinning parameter");
		if(sscanf(argv[pos++],"%d", &NoPP)!=1)
			Exception::Throw("Failed to read NoPP parameter");
		if(sscanf(argv[pos++],"%d", &NoICSdiscr)!=1)
			Exception::Throw("Failed to read NoICSdiscr parameter");
		if(sscanf(argv[pos++],"%d", &NoICScel)!=1)
			Exception::Throw("Failed to read NoICScel parameter");
		if(sscanf(argv[pos++],"%d", &NoRedshift)!=1)
			Exception::Throw("Failed to read NoRedshift parameter");
		if(sscanf(argv[pos++],"%d", &split)!=1)
			Exception::Throw("Failed to read split parameter");
		const char* outputDir = argc==13 ? argv[pos++] : 0;

		KachelriesTest(Emin, E, z, alpha, split!=0, N, NoPP!=0, NoICSdiscr!=0, NoICScel!=0, NoRedshift!=0, pt, outputDir);
	}
	else if(strcmp(command,"mono")==0)
	{
		//MonochromaticAndGausianBackgroundTest(true,true,false,false);
		MonochromaticAndGausianBackgroundTest(true,false,true,false);
		//MonochromaticAndGausianBackgroundTest(true,false,false,true);
		//MonochromaticAndGausianBackgroundTest(false,true,false,false);
		MonochromaticAndGausianBackgroundTest(false,false,true,false);
		//MonochromaticAndGausianBackgroundTest(false,false,false,true);
	}
	else if(strcmp(command,"OpenMP")==0)
	{
		MultithreadTest(false);
		MultithreadTest(true);
	}
	else if(strcmp(command,"bi")==0)
	{
		BichromaticTest(true);
	}
	else if(strcmp(command,"sampling")==0)
	{
		SamplingTest();
	}
	else if(strcmp(command,"engine")==0)
	{
		EngineTest(0.5);
	}
	return 0;
}

void TestList::SamplingTest()
{
	ParamlessFunction aDistrib(sin);
	double aRand = 0.5;
	double aOutputX;
	double aOutputIntegral;
	double xMin = 0;
	double xMax = M_PI;
	double aRelError = 1e-4;
	//SampleDistribution(const Function& aDistrib, double aRand, double& aOutputX, double& aOutputIntegral, double aInitialStep, double xMin, double xMax, double aRelError);
	MathUtils::SampleDistribution(aDistrib, aRand, aOutputX, aOutputIntegral, xMin, xMax, aRelError);

	aOutputX/=M_PI;
	aOutputX*=M_PI;


	class SamplingTestFunc : public Function
	{
		GaussianBackground fTestFunc;
	public:
		SamplingTestFunc():
			fTestFunc(1e10, 1e8, 99, 1)
		{

		}
		double f(double aX) const
		{
			return fTestFunc.n(aX, 0.)/aX;
		}
		virtual double Xmin() const {return fTestFunc.MinE();}
		virtual double Xmax() const {return fTestFunc.MaxE();}
	};
	std::ofstream samplingOut;
	samplingOut.open("samplingOld",std::ios::out);
	std::ofstream samplingOutNew;
	samplingOutNew.open("sampling",std::ios::out);
	Randomizer r(0);
	SamplingTestFunc f;
	std::ofstream funcOut;
	funcOut.open("comulative.f",std::ios::out);
	double step = pow(f.Xmax()/f.Xmin(),0.001);
	double sum = 0;
	funcOut << f.Xmin() << "\t0\n";
	MathUtils math;
	for(double x = f.Xmin(); x<f.Xmax(); x*=step)
	{
		sum+=math.Integration_qag(f,x,x*step,1e-300,1e-5,1000);
		funcOut << x << "\t" << sum << "\n";
	}
	funcOut.close();

	int nSteps = 100;
	xMin=10;
	xMax=1e20;
	aRelError = 1e-3;
	double x;
	double I;
	int nTrials = 100000;
	for(int i=0; i<nTrials; i++)
	{
		double rand = r.Rand();
		math.SampleLogscaleDistribution<double>(f, rand, x, I, nSteps, xMin, xMax, aRelError);
		samplingOut << x << "\n";
		MathUtils::SampleLogDistribution(f, rand, x, I, xMin, xMax, aRelError);
		samplingOutNew << x << "\n";

	}
	samplingOut.close();
	samplingOutNew.close();
	std::cout << "prepare randomly generated comulative distribution file with command:\n"
			<<"sort -g sampling | awk ' {SUM=SUM+1; printf(\"%g\\t%d\\n\",$1,SUM)}' > comulative.r\n"
			<<"and compare it to comulative.f using gnuplot command\n\t"
			<<"pl 'comulative.r' u 1:(" << 1.0/nTrials <<  "*$2) w l, 'comulative.f'" << endl;

}

void TestList::BackgroundTest()
{
	double zMax = 7;
	//cosmology.init(0.73, 71, 10);

	CompoundBackground backgr;
	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	IBackground* plank = new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp, 0, zMax);
	IBackground* kneiske = new Kneiske0309EBL();
	double conc = 413*units.Vunit;//413 cm^{-3}
	double centralE = 6.3e-10/units.Eunit;
	GaussianBackground backgrGauss(centralE, 0.1*centralE, 1., conc);
	backgr.AddComponent(plank);
	backgr.AddComponent(kneiske);
	ofstream plankout;
	plankout.open("backgr_plank",ios::out);
	BackgroundTest(*plank, plankout);
	ofstream kneiskeout;
	kneiskeout.open("backgr_kneiskeBF",ios::out);
	BackgroundTest(*kneiske, kneiskeout);
	ofstream gaussout;
	gaussout.open("backgr_gauss",ios::out);
	BackgroundTest(backgrGauss, gaussout);
	ofstream sumout;
	sumout.open("backgr_sum",ios::out);
	BackgroundTest(backgr, sumout);
}

void TestList::BackgroundTest(const IBackground& aBackground, std::ostream& aOut)
{
	int kIntervals = 100;
	int zIntervals = aBackground.MaxZ()-aBackground.MinZ();
	if(zIntervals<1)
		zIntervals = 2;
	double zStep = (aBackground.MaxZ()-aBackground.MinZ())/zIntervals;//usually zStep=1;
	double kStep = pow(aBackground.MaxE()/aBackground.MinE(), 1./kIntervals);
	double z = aBackground.MinZ();
	aOut << "# " << aBackground.Name() << " background\n#col1: E/eV\t\tcol2: k dn(k)/dk /cm^{-3}\n";
	for(int iZ=0; iZ<=zIntervals; iZ++, z+=zStep)
	{
		aOut << "#z=" << z << "\n";
		double k = aBackground.MinE();
		for(int iK=0; iK<=kIntervals; iK++, k*=kStep)
			aOut << k*units.Eunit*1e6 << "\t" << aBackground.n(k,z)/units.Vunit << "\n";
		aOut << "\n\n";
	}
}


void TestList::SplittedBackgroundTest(bool aSplitted) {
	if(!cosmology.IsInitialized())
		cosmology.Init();
	double Zmax = 0.15;
	double Emin = 1./units.Eunit;
	double alphaThinning = 1;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	result.AddOutput(new SpectrumOutput(aSplitted?"out_split" : "out_no_split", Emin, pow(10, 0.05)));

	CompoundBackground backgr;

	//double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	IBackground* b1 = new GaussianBackground(1e-6/units.Eunit, 1e-7/units.Eunit, 1., 1*units.Vunit);
	IBackground* b2 = new GaussianBackground(6.3e-10/units.Eunit, 1e-11/units.Eunit, 1., 413*units.Vunit);
	//IBackground* b1 = new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp);
	//IBackground* b2 = new Kneiske1001EBL();
	//double n = BackgroundUtils::CalcIntegralDensity(*b1)/units.Vunit;

	backgr.AddComponent(b1);
	backgr.AddComponent(b2);


	double stepZ = Zmax<0.2 ? Zmax/2 : 0.1;
	double logStepK = pow(10,0.05);
	double epsRel = 1e-3;

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel);
	PBackgroundIntegral backgrI1 = new ContinuousBackgroundIntegral(*b1, stepZ, logStepK, Zmax, epsRel);
	PBackgroundIntegral backgrI2 = new ContinuousBackgroundIntegral(*b1, stepZ, logStepK, Zmax, epsRel);


	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	pe.AddInteraction(new RedShift());
	if(aSplitted)
	{
		pe.AddInteraction(new Interactions::ICS(backgrI1,Emin));
		pe.AddInteraction(new Interactions::IcsCEL(backgrI1,Emin));
		pe.AddInteraction(new Interactions::GammaPP(backgrI1));
		pe.AddInteraction(new Interactions::ICS(backgrI2,Emin));
		pe.AddInteraction(new Interactions::IcsCEL(backgrI2,Emin));
		pe.AddInteraction(new Interactions::GammaPP(backgrI2));
	}
	else
	{
		pe.AddInteraction(new Interactions::ICS(backgrI,Emin));
		pe.AddInteraction(new Interactions::IcsCEL(backgrI,Emin));
		pe.AddInteraction(new Interactions::GammaPP(backgrI));
	}

	int nParticles = 10000;
	for(int i=0; i<nParticles; i++)
	{
		Particle gamma(Photon, Zmax);
		gamma.Energy = 1e8/units.Eunit;//MeV
		particles.AddPrimary(gamma);
	}

	pe.Run();
}

void TestList::ThinningTest(double aAlpha)
{
	if(!cosmology.IsInitialized())
		cosmology.Init();
	double logStepK = pow(10,0.05);
	double Zmax = 0.01;
	double Emin = 1./units.Eunit;
	double alphaThinning = aAlpha;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	result.AddOutput(new SpectrumOutput("out_alpha" + ToString(aAlpha), Emin, logStepK));

	CompoundBackground backgr;

	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	//IBackground* b1 = new GaussianBackground(1e-6/units.Eunit, 1e-7/units.Eunit, 1., 1*units.Vunit);
	//IBackground* b2 = new GaussianBackground(6.3e-10/units.Eunit, 1e-11/units.Eunit, 1., 413*units.Vunit);
	IBackground* b1 = new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp);
	//IBackground* b2 = new Kneiske1001EBL();
	//double n = BackgroundUtils::CalcIntegralDensity(*b1)/units.Vunit;

	backgr.AddComponent(b1);
	//backgr.AddComponent(b2);


	double stepZ = Zmax<0.2 ? Zmax/2 : 0.1;
	double epsRel = 1e-3;

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel);


	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	pe.AddInteraction(new RedShift());

	{
		pe.AddInteraction(new Interactions::ICS(backgrI,Emin));
		pe.AddInteraction(new Interactions::IcsCEL(backgrI,Emin));
		pe.AddInteraction(new Interactions::GammaPP(backgrI));
	}

	int nParticles = 1000;
	for(int i=0; i<nParticles; i++)
	{
		Particle gamma(Photon, Zmax);
		gamma.Energy = 1e8/units.Eunit;//MeV
		particles.AddPrimary(gamma);
	}

	pe.Run();
}

void TestList::MonochromaticAndGausianBackgroundTest(bool aMonochromatic, bool aICS, bool aICScel, bool aPP)
{
	if(!cosmology.IsInitialized())
		cosmology.Init();
	double logStepK = pow(10,0.05);
	double Zmax = 0.01;
	double Emin = 1./units.Eunit;
	double alphaThinning = 0.5;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	std::string out = aMonochromatic ? "out_mono" : "out_gauss";
	if(aICS)
		out = out + "_ICS";
	if(aICScel)
			out = out + "_ICScel";
	if(aPP)
			out = out + "_PP";

	result.AddOutput(new SpectrumOutput(out , Emin, logStepK));

	double conc = 413*units.Vunit;//413 cm^{-3}
	double centralE = 6.3e-10/units.Eunit;
	ConstFunction k(centralE);
	ConstFunction c(conc);

	GaussianBackground backgr(centralE, 0.1*centralE, 1., conc);

	double stepZ = Zmax<0.2 ? Zmax/2 : 0.1;
	double epsRel = 1e-3;

	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	pe.AddInteraction(new RedShift());
	PBackgroundIntegral backgrM = new MonochromaticBackgroundIntegral(c, k, logStepK, epsRel);
	PBackgroundIntegral backgrC = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel);
	BackgroundIntegral* backgrI = aMonochromatic ? backgrM : backgrC;

	if(aICS)
		pe.AddInteraction(new Interactions::ICS(backgrI,Emin));
	if(aICScel)
		pe.AddInteraction(new Interactions::IcsCEL(backgrI,Emin));
	if(aPP)
		pe.AddInteraction(new Interactions::GammaPP(backgrI));

	int nParticles = 3000;
	for(int i=0; i<nParticles; i++)
	{
		Particle startParticle(aPP?Photon:Electron, Zmax);
		startParticle.Energy = 1e9/units.Eunit;//MeV
		particles.AddPrimary(startParticle);
	}

	//pe.RunMultithread();
	pe.Run();
}

//
void TestList::BichromaticTest(bool aGauss)
{
	if(!cosmology.IsInitialized())
		cosmology.Init();
	double logStepK = pow(10,0.05);
	double Zmax = 1.;
	double Emin = 1./units.Eunit;
	double alphaThinning = 0.5;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	std::string out = "bichromatic_";
	out += aGauss ? "gauss" : "mono";

	result.AddOutput(new SpectrumOutput(out , Emin, logStepK));

	double centralE1 = 6.3e-10/units.Eunit;//6.3e-4eV
	double n1 = 413*units.Vunit;//413cm^-3
	ConstFunction k1(centralE1);
	ConstFunction c1(n1);
	GaussianBackground g1(centralE1, 0.1*centralE1, 1., n1);

	double centralE2 = 1e-6/units.Eunit;//1 eV
	double n2 = 1.*units.Vunit;//1 cm^-3
	ConstFunction k2(centralE2);
	ConstFunction c2(n2);
	GaussianBackground g2(centralE2, 0.1*centralE2, 1., n2);

	double epsRel = 1e-3;

	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	//pe.AddInteraction(new RedShift());
	PBackgroundIntegral I1,I2;
	if(aGauss)
	{
		I1 = new ContinuousBackgroundIntegral(g1, 0.2, logStepK, Zmax, epsRel);
		I2 = new ContinuousBackgroundIntegral(g2, 0.2, logStepK, Zmax, epsRel);
	}
	else
	{
		I1 = new MonochromaticBackgroundIntegral(c1, k1, logStepK, epsRel);
		I2 = new MonochromaticBackgroundIntegral(c2, k2, logStepK, epsRel);
	}
	pe.AddInteraction(new Interactions::ICS(I1,Emin));
	pe.AddInteraction(new Interactions::IcsCEL(I1,Emin));
	pe.AddInteraction(new Interactions::GammaPP(I1));
	pe.AddInteraction(new Interactions::ICS(I2,Emin));
	pe.AddInteraction(new Interactions::IcsCEL(I2,Emin));
	pe.AddInteraction(new Interactions::GammaPP(I2));

	int nParticles = 1000;
	for(int i=0; i<nParticles; i++)
	{
		Particle startParticle(Photon, Zmax);
		startParticle.Energy = 1e9/units.Eunit;//MeV
		particles.AddPrimary(startParticle);
	}

	pe.RunMultithread();
	//pe.Run();
}

void TestList::MultithreadTest(bool aMultithread)
{
	if(!cosmology.IsInitialized())
		cosmology.Init();
	double logStepK = pow(10,0.05);
	double Zmax = 0.01;
	double Emin = 1./units.Eunit;
	double alphaThinning = 0;//0.5;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	result.AddOutput(new SpectrumOutput(aMultithread ? "out_multithread" : "out_singlethread" , Emin, logStepK));

	double conc = 413*units.Vunit;//413 cm^{-3}
	double centralE = 6.3e-10/units.Eunit;
	ConstFunction k(centralE);
	ConstFunction c(conc);

	GaussianBackground backgr(centralE, 0.1*centralE, 1., conc);

	double stepZ = Zmax<0.2 ? Zmax/2 : 0.1;
	double epsRel = 1e-3;

	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	pe.AddInteraction(new RedShift());
	PBackgroundIntegral backgrM = new MonochromaticBackgroundIntegral(c, k, logStepK, epsRel);
	PBackgroundIntegral backgrC = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel);

	pe.AddInteraction(new Interactions::ICS(backgrM,Emin));
	pe.AddInteraction(new Interactions::IcsCEL(backgrM,Emin));
	pe.AddInteraction(new Interactions::GammaPP(backgrM));
	pe.AddInteraction(new Interactions::ICS(backgrC,Emin));
	pe.AddInteraction(new Interactions::IcsCEL(backgrC,Emin));
	pe.AddInteraction(new Interactions::GammaPP(backgrC));

	int nParticles = 30;
	for(int i=0; i<nParticles; i++)
	{
		Particle gamma(Photon, Zmax);
		gamma.Energy = 1e9/units.Eunit;//MeV
		particles.AddPrimary(gamma);
	}

	if(aMultithread)
		pe.RunMultithread();
	else
		pe.Run();
}

//this test is used to compare calculations of this code with Elmag http://uk.arxiv.org/pdf/1106.5508v2
void TestList::KachelriesTest(double aEmin_eV, double aEmax_eV, double aZmax,
		double aAlpha, bool aSplitted, int nParticles,
		bool noPP, bool noICS, bool noICScel, bool aNoRedshift,
		ParticleType primary, const char* aOutputDir)
{
	int kAc = 1;
	double Emax = aEmax_eV;//eV
	Emax *= (1e-6/units.Eunit);
	if(!cosmology.IsInitialized())
		cosmology.Init();// Elmag 2.0 uses H=70 km/sec/Mpc (I checked this by measuring dt/dz at z=0
	                                 // in subroutine cascade(icq,e00,weight0,z_in)
									 // H = z_in/((1.d0-t)*t_0)/3.2407792903e-20
									 // Omega_vac = 0.7 as it was stated in http://uk.arxiv.org/pdf/1106.5508v2 page 6 (section 2.6. Cosmology)
	double logStepK = pow(10,0.05/kAc);
	double Emin = aEmin_eV*1e-6/units.Eunit;
	double alphaThinning = aAlpha;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	std::string outputDir;
	if(aOutputDir)
		outputDir = aOutputDir;
	else
	{
		outputDir = Particle::Name(primary);
		outputDir += ("_E" + ToString(aEmax_eV) + "_z" + ToString(aZmax) + "_N" + ToString(nParticles));
		if(noPP)
			outputDir += "_NoPP";
		if(noICS)
			outputDir += "_NoICSdiscr";
		if(noICScel)
			outputDir += "_NoICScel";
		if(aNoRedshift)
			outputDir += "_NoRedshift";
		if(aSplitted)
			outputDir += "_splited";
	}

	result.AddOutput(new RawOutput(outputDir, true));
	result.AddOutput(new SpectrumOutput(outputDir + "/spec", Emin, pow(10,0.05)));

	CompoundBackground backgr;

	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	//IBackground* b1 = new GaussianBackground(1e-6/units.Eunit, 1e-7/units.Eunit, 1., 1*units.Vunit);
	//IBackground* b2 = new GaussianBackground(6.3e-10/units.Eunit, 1e-11/units.Eunit, 1., 413*units.Vunit);
	IBackground* b1 = new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp);
	IBackground* b2 = new ElmagKneiskeBestFit(false/*using wrong background intentionally to compare to Elmag*/);
	//double n = BackgroundUtils::CalcIntegralDensity(*b1)/units.Vunit;

	backgr.AddComponent(b1);
	backgr.AddComponent(b2);


	double stepZ = aZmax<0.05 ? aZmax/2 : 0.025;
	double epsRel = 1e-3;

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, aZmax, epsRel);
	PBackgroundIntegral backgrI1 = new ContinuousBackgroundIntegral(*b1, stepZ, logStepK, aZmax, epsRel);
	PBackgroundIntegral backgrI2 = new ContinuousBackgroundIntegral(*b2, stepZ, logStepK, aZmax, epsRel);

	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	if(!aNoRedshift)
		pe.AddInteraction(new RedShift());
	Particle electron(Electron, aZmax);
	Particle photon(Photon, aZmax);
	int nSteps = log(Emax/Emin)/log(logStepK)+1.;
	double mult = pow(Emax/Emin, 1./nSteps);
	int step;
	if(aSplitted)
	{
		if(!noICS)
		{
			std::ofstream rateOut;
			rateOut.open((outputDir + "/ICS").c_str(),std::ios::out);
			PRandomInteraction i1,i2;
			pe.AddInteraction((i1 = new Interactions::ICS(backgrI1,Emin)));
			pe.AddInteraction((i2 = new Interactions::ICS(backgrI2,Emin)));
			for(step=0, electron.Energy=Emin; step<=nSteps; electron.Energy*=mult, step++)
			{
				double r1 = i1->Rate(electron)*units.Mpc;
				double r2 = i2->Rate(electron)*units.Mpc;
				rateOut << electron.Energy * units.Eunit * 1e6 << "\t" << r1+r2 << "\t" << r1 << "\t" << r2 << "\n";
			}
			rateOut.close();
		}
		if(!noICScel)
		{
			std::ofstream rateOut;
			rateOut.open((outputDir + "/ICScel").c_str(),std::ios::out);
			PCELInteraction i1,i2;
			pe.AddInteraction((i1 = new Interactions::IcsCEL(backgrI1,Emin)));
			pe.AddInteraction((i2 = new Interactions::IcsCEL(backgrI2,Emin)));
			for(step=0, electron.Energy=Emin; step<=nSteps; electron.Energy*=mult, step++)
			{
				double r1 = i1->Rate(electron)*units.Mpc;
				double r2 = i2->Rate(electron)*units.Mpc;
				rateOut << electron.Energy * units.Eunit * 1e6 << "\t" << r1+r2 << "\t" << r1 << "\t" << r2 << "\n";
			}
			rateOut.close();
		}
		if(!noPP)
		{
			std::ofstream rateOut;
			rateOut.open((outputDir + "/PP").c_str(),std::ios::out);
			PRandomInteraction i1,i2;
			pe.AddInteraction((i1 = new Interactions::GammaPP(backgrI1)));
			pe.AddInteraction((i2 = new Interactions::GammaPP(backgrI2)));
			for(step=0, photon.Energy=Emin; step<=nSteps; photon.Energy*=mult, step++)
			{
				double r1 = i1->Rate(photon)*units.Mpc;
				double r2 = i2->Rate(photon)*units.Mpc;
				rateOut << photon.Energy * units.Eunit * 1e6 << "\t" << r1+r2 << "\t" << r1 << "\t" << r2 << "\n";
			}
			rateOut.close();
		}
	}
	else
	{
		if(!noICS)
		{
			std::ofstream rateOut;
			rateOut.open((outputDir + "/ICS").c_str(),std::ios::out);
			PRandomInteraction i = new Interactions::ICS(backgrI,Emin);
			pe.AddInteraction(i);
			for(step=0, electron.Energy=Emin; step<=nSteps; electron.Energy*=mult, step++)
				rateOut << electron.Energy * units.Eunit * 1e6 << "\t" << i->Rate(electron)*units.Mpc << "\n";
			rateOut.close();
		}
		if(!noICScel)
		{
			std::ofstream rateOut;
			rateOut.open((outputDir + "/ICScel").c_str(),std::ios::out);
			PCELInteraction i = new Interactions::IcsCEL(backgrI,Emin);
			pe.AddInteraction(i);
			for(step=0, electron.Energy=Emin; step<=nSteps; electron.Energy*=mult, step++)
				rateOut << electron.Energy * units.Eunit * 1e6 << "\t" << i->Rate(electron)*units.Mpc << "\n";
			rateOut.close();
		}
		if(!noPP)
		{
			std::ofstream rateOut;
			rateOut.open((outputDir + "/PP").c_str(),std::ios::out);
			PRandomInteraction i = new Interactions::GammaPP(backgrI);
			pe.AddInteraction(i);
			for(step=0, photon.Energy=Emin; step<=nSteps; photon.Energy*=mult, step++)
				rateOut << photon.Energy * units.Eunit * 1e6 << "\t" << i->Rate(photon)*units.Mpc << "\n";
			rateOut.close();
		}
	}

	for(int i=0; i<nParticles; i++)
	{
		//if(i%100==0)
		//	std::cerr << i << " of " << nParticles << " ready" << std::endl;
		Particle gamma(primary, aZmax);
		gamma.Energy = Emax;//MeV
		particles.AddPrimary(gamma);
		//pe.RunMultithread();
	}
#ifdef _DEBUG
	pe.Run();
#else
	pe.RunMultithread();
	//pe.Run();
#endif
}

void FermiTests::Main(int argc, char** argv)
{
	if(argc==0)
	{
		PrintUsage();
	}
	const char* command = argv[0];
	if(argc == 6 && !strcmp(command,"Absorption"))
	{
		int aInfiniteB;
		double aEmax;
		double aSpecAlpha;
		int aEBLmodel;
		double aDistanceMpc;
		int par=1;
		if(sscanf(argv[par++], "%d", &aInfiniteB)!=1)
			Exception::Throw("Failed to read InfiniteB");
		if(sscanf(argv[par++], "%lg", &aSpecAlpha)!=1)
			Exception::Throw("Failed to read SpecAlpha");
		if(sscanf(argv[par++], "%d", &aEBLmodel)!=1)
			Exception::Throw("Failed to read EBLmodel");
		if(sscanf(argv[par++], "%lg", &aEmax)!=1)
			Exception::Throw("Failed to read Emax");
		if(sscanf(argv[par++], "%lg", &aDistanceMpc)!=1)
			Exception::Throw("Failed to read DistanceMpc");
		AbsorptionTest(aInfiniteB!=0, aSpecAlpha, aEBLmodel, aEmax*1e-6/units.Eunit, aDistanceMpc);
	}
	else if(argc == 5 && !strcmp(command,"DiffuseToPoint"))
	{
		double aSpecAlpha;
		int aEBLmodel;
		double aEmax;
		int par=1;
		int aEvolutionModel = -1;
		//if(sscanf(argv[par++], "%d", &aInfiniteB)!=1)
		//	Exception::Throw("Failed to read InfiniteB");
		if(sscanf(argv[par++], "%lg", &aSpecAlpha)!=1)
			Exception::Throw("Failed to read SpecAlpha");
		if(sscanf(argv[par++], "%d", &aEBLmodel)!=1)
			Exception::Throw("Failed to read EBLmodel");
		if(sscanf(argv[par++], "%lg", &aEmax)!=1)
			Exception::Throw("Failed to read Emax");
		if(sscanf(argv[par++], "%d", &aEvolutionModel)!=1)
			Exception::Throw("Failed to read EvolutionModel");
		DiffuseToPointTest(aSpecAlpha, aEBLmodel, aEmax*1e-6/units.Eunit, aEvolutionModel);
	}
	else
		PrintUsage();
}

void FermiTests::PrintUsage()
{
	cerr << "FERMI test commands syntax:" << std::endl <<
	"Absorption InfiniteB SpecAlpha EBLmodel Emax/eV DistanceMpc" << std::endl << "\tor" << std::endl <<
	"DiffuseToPoint SpecAlpha EBLmodel Emax/eV EvolutionModel" << std::endl;
}

void FermiTests::InitCosmology()
{
	if(!cosmology.IsInitialized())
		cosmology.Init();
}

void FermiTests::InitSpec(int nParticles, ParticleStack& aStack, double aSpecAlpha, double aEmin, double aEmax, double aZ, double aWeight, Result* aInitialOutput)
{
	double logStep = pow(aEmax/aEmin, 1./(double)nParticles);
	Particle gamma(Photon, aZ);
	double logScaleAlpha = aSpecAlpha - 1.;
	double E = aEmin*logStep;
	for(int i=0; i<nParticles; i++, E*=logStep)
	{
		gamma.Energy = E;
		gamma.Weight = aWeight*pow(E/aEmin, -logScaleAlpha);
		aStack.AddPrimary(gamma);
		if(aInitialOutput)
			aInitialOutput->Add(gamma);
	}
}

Particle FermiTests::SamplePowerLawSpec(double aSpecAlpha, double aEmin, double aEmax, double aZ, double aWeight, double aRand)
{
	double logScaleAlpha = aSpecAlpha - 1.;
	double E = aEmin * pow(aEmax/aEmin, aRand);//sampling uniform in logscale
	Particle gamma(Photon, aZ);
	gamma.Energy = E;
	gamma.Weight = aWeight*pow(E/aEmin, -logScaleAlpha);//power law spectrum
	return gamma;
}

IBackground* FermiTests::CreateEBL(int aEBLmodel)
{
	IBackground* b2 = 0;
	if(aEBLmodel==8)
		b2 = new Kneiske1001EBL();
	else if(aEBLmodel==9)
		b2 = new Kneiske0309EBL();
	else if(aEBLmodel!=0)
		Exception::Throw("Unsupported EBLmodel parameter: expected values 0-None, 8-minimal, 9-best fit");
	return b2;
}

const double FermiTests::MinE = 10000;//MeV

void FermiTests::AbsorptionTest(bool aInfiniteB, double aSpecAlpha, int aEBLmodel, double aEmax, double aDistanceMpc)
{
	int nParticles = 30000;
	InitCosmology();
	double aZmax = cosmology.d2z(aDistanceMpc*units.Mpc);

	double EminOutput = MinE/units.Eunit;
	double EminCalc = 0.99*EminOutput;
	double alphaThinning = 0.5;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;

	std::string modelStr = ToString(aDistanceMpc) + "Mpc_B" + (aInfiniteB ? "1" : "0") +
			 + "_Emax" + ToString(aEmax/units.eV) + "_Alpha" + ToString(aSpecAlpha);
	Result resultFin(EminCalc);
	resultFin.AddOutput(new SpectrumOutput(modelStr + "_final_spec" , EminOutput, pow(10,0.05)));
	SmartPtr<IntegralSpectrumOutput> finalIntSpec = new IntegralSpectrumOutput(modelStr + "_final_int_spec", 100./units.Eunit);
	resultFin.AddOutput(finalIntSpec);

	Result resultIni(EminCalc);
	resultIni.AddOutput(new SpectrumOutput(modelStr + "_ini_spec" , EminOutput, pow(10,0.05)));
	SmartPtr<IntegralSpectrumOutput> iniIntSpec = new IntegralSpectrumOutput(modelStr + "_ini_int_spec", 100./units.Eunit);
	resultIni.AddOutput(iniIntSpec);

	PropagationEngine pe(particles, resultFin);
	EnergyBasedThinning thinning(alphaThinning);
	InitPropagEngine(pe, aInfiniteB, aEBLmodel, &thinning, aZmax, EminCalc);

	InitSpec(nParticles, particles, aSpecAlpha, EminOutput, aEmax, aZmax, 1., &resultIni);
	//Randomizer r;
	//for(; nParticles>0; nParticles--)
	//{
	//	Particle p = SamplePowerLawSpec(aSpecAlpha, EminOutput, aEmax, aZmax, 1., r.Rand());
	//	particles.Add(p);
	//	resultIni.Add(p);
	//}


	resultIni.Flush();

	//main calculation cycle
	pe.RunMultithread();
	//pe.Run();

	//output integral and diff spectra
	resultFin.Flush();

	//output absorption factor
	SmartPtr<TableFunction> fIni = iniIntSpec->IntegralSpectrum();
	SmartPtr<TableFunction> fFin = finalIntSpec->IntegralSpectrum();
	double step = pow(10,0.05);
	std::ofstream absOut;
	absOut.open((modelStr + "_abs").c_str(),std::ios::out);

	//fIni->Print(std::cout);

	double aEmaxeV = aEmax*1e6*units.Eunit;
	for(double E=fIni->Xmin(); E<=fIni->Xmax(); E*=step)
	{
		double initialF = fIni->f(E);
		ASSERT(initialF>=0);
		if(initialF==0.)
			continue;
		double absorbtion = fFin->f(E)/fIni->f(E);
		absOut << aEmaxeV <<  "\t" << aEBLmodel <<
				"\t" << aSpecAlpha << "\t" << (aInfiniteB?1:0) <<
				"\t" << aDistanceMpc << "\t" << E << "\t" << absorbtion << std::endl;
	}
	absOut.close();
}

void FermiTests::InitPropagEngine(PropagationEngine& pe, bool aInfiniteB, int aEBLmodel, IThinning* aThinning, double aZmax, double EminCalc)
{
	bool splitted = false;
	double logStepK = pow(10,0.05);
	CompoundBackground backgr;

	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	IBackground* b1 = new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp);
	IBackground* b2 = CreateEBL(aEBLmodel);

	backgr.AddComponent(b1);//transfers ownership
	if(b2)
		backgr.AddComponent(b2);

	double stepZ = aZmax<0.2 ? aZmax/2 : 0.1;
	double epsRel = 1e-3;

	PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, aZmax, epsRel);
	PBackgroundIntegral backgrI1 = new ContinuousBackgroundIntegral(*b1, stepZ, logStepK, aZmax, epsRel);
	PBackgroundIntegral backgrI2 = b2?new ContinuousBackgroundIntegral(*b2, stepZ, logStepK, aZmax, epsRel):0;


	pe.SetThinning(aThinning);
	pe.AddInteraction(new RedShift());

	if(splitted)
	{
		if(!aInfiniteB)
		{
			pe.AddInteraction(new Interactions::ICS(backgrI1,EminCalc));
			pe.AddInteraction(new Interactions::IcsCEL(backgrI1,EminCalc));
			if(b2)
			{
				pe.AddInteraction(new Interactions::ICS(backgrI2,EminCalc));
				pe.AddInteraction(new Interactions::IcsCEL(backgrI2,EminCalc));
			}
		}
		pe.AddInteraction(new Interactions::GammaPP(backgrI1));
		if(b2)
			pe.AddInteraction(new Interactions::GammaPP(backgrI2));
	}
	else
	{
		if(!aInfiniteB)
		{
			pe.AddInteraction(new Interactions::ICS(backgrI,EminCalc));
			pe.AddInteraction(new Interactions::IcsCEL(backgrI,EminCalc));
		}
		pe.AddInteraction(new Interactions::GammaPP(backgrI));
	}
}

void FermiTests::DiffuseToPointTest(double aSpecAlpha, int aEBLmodel, double aEmax, int aEvolutionModel)
{
	time_t tStart = time(0);
	int nTotalParticles = 30000;
	double ZmaxFull = 3;
	double ZmaxPoint = 0.06;

	InitCosmology();
	double EminOutput = MinE/units.Eunit;
	double EminCalc = 0.99*EminOutput;
	double alphaThinning = 0.5;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particlesPointB0;
	ParticleStack particlesPointB1;
	ParticleStack particlesFull;

	std::string modelStr = //ToString(aInfiniteB ? "B1_" : "B0_") +
			ToString("Emax") + ToString(aEmax/units.eV) + "_Alpha" + ToString(aSpecAlpha) + "_Evol" + ToString(aEvolutionModel);
	Result resultFull(EminCalc);
	resultFull.AddOutput(new SpectrumOutput(modelStr + "_full_spec" , EminOutput, pow(10,0.05)));
	SmartPtr<IntegralSpectrumOutput> fullIntSpec = new IntegralSpectrumOutput(modelStr + "_full_int_spec", 100./units.Eunit);
	resultFull.AddOutput(fullIntSpec);

	Result resultPointB0(EminCalc);
	resultPointB0.AddOutput(new SpectrumOutput(ToString("B0_") + modelStr + "_point_spec" , EminOutput, pow(10,0.05)));
	SmartPtr<IntegralSpectrumOutput> pointB0IntSpec = new IntegralSpectrumOutput(ToString("B0_") + modelStr + "_point_int_spec", 100./units.Eunit);
	resultPointB0.AddOutput(pointB0IntSpec);

	Result resultPointB1(EminCalc);
	resultPointB1.AddOutput(new SpectrumOutput(ToString("B1_") + modelStr + "_point_spec" , EminOutput, pow(10,0.05)));
	SmartPtr<IntegralSpectrumOutput> pointB1IntSpec = new IntegralSpectrumOutput(ToString("B1_") + modelStr + "_point_int_spec", 100./units.Eunit);
	resultPointB1.AddOutput(pointB1IntSpec);

	PropagationEngine pePointB0(particlesPointB0, resultPointB0);
	PropagationEngine pePointB1(particlesPointB1, resultPointB1);
	PropagationEngine peFull(particlesFull, resultFull);
	EnergyBasedThinning thinning(alphaThinning);
	InitPropagEngine(peFull, false, aEBLmodel, &thinning, ZmaxFull, EminCalc);
	InitPropagEngine(pePointB0, false, aEBLmodel, &thinning, ZmaxFull, EminCalc);
	InitPropagEngine(pePointB1, true, aEBLmodel, &thinning, ZmaxFull, EminCalc);

	int nEnergyBins=100;
	int nParticlesClose = nTotalParticles/2;
	int nParticlesFar = nTotalParticles-nParticlesClose;

	double tEarliest = cosmology.z2t(ZmaxFull);
	double t0 = cosmology.z2t(0.);
	double tPoint = cosmology.z2t(ZmaxPoint);

	double weightFar = (tPoint-tEarliest)/(t0-tPoint)*nParticlesClose/nParticlesFar;

	Randomizer r;
	for(int nSources = nParticlesClose/nEnergyBins; nSources>0; nSources--)
	{//sampling close sources
		double t = tPoint + r.Rand()*(t0-tPoint);
		double z = cosmology.t2z(t);
		double weight = Evolution::EvolFactor(z, aEvolutionModel);
		InitSpec(nEnergyBins, particlesPointB0, aSpecAlpha, EminOutput, aEmax, z, weight);
		InitSpec(nEnergyBins, particlesPointB1, aSpecAlpha, EminOutput, aEmax, z, weight);
		InitSpec(nEnergyBins, particlesFull, aSpecAlpha, EminOutput, aEmax, z, weight);
	}
	for(int nSources = nParticlesFar/nEnergyBins; nSources>0; nSources--)
	{//sampling far away sources
		double t = tEarliest + r.Rand()*(tPoint-tEarliest);
		double z = cosmology.t2z(t);
		double weight = Evolution::EvolFactor(z, aEvolutionModel)*weightFar;
		InitSpec(nEnergyBins, particlesFull, aSpecAlpha, EminOutput, aEmax, z, weight);
	}

	std::cerr << "\n starting B0 calculation on " << time(0) - tStart << " sec" << std::endl;

	//calculate point source spectrum
	pePointB0.RunMultithread();
	resultPointB0.Flush();

	std::cerr << "\n starting B1 calculation on " << time(0) - tStart << " sec" << std::endl;

	//calculate point source spectrum
	pePointB1.RunMultithread();
	resultPointB1.Flush();

	std::cerr << "\n starting diffuse calculation on " << time(0) - tStart << " sec" << std::endl;

	//calculate full spectrum
	peFull.RunMultithread();
	resultFull.Flush();

	std::cerr << "\n starting diffuse to point ratio calculation on " << time(0) - tStart << " sec" << std::endl;

	//output diffuse to point fraction
	SmartPtr<TableFunction> fPointB0 = pointB0IntSpec->IntegralSpectrum();
	SmartPtr<TableFunction> fPointB1 = pointB1IntSpec->IntegralSpectrum();
	SmartPtr<TableFunction> fFull = fullIntSpec->IntegralSpectrum();
	double step = pow(10,0.05);
	double aEmaxeV = aEmax*1e6*units.Eunit;
	TableFunction* intFpoint[2] = {fPointB0, fPointB1};
	for(int i=0; i<2; i++)
	{
		std::ofstream ratioOut;
		ratioOut.open((ToString("B") + ToString(i) + "_" + modelStr + "_diff_to_point").c_str(),std::ios::out);
		TableFunction* fPoint = intFpoint[i];
		for(double E=fPoint->Xmin(); E<=fPoint->Xmax(); E*=step)
		{
			double pointF = fPoint->f(E);
			ASSERT(pointF>=0);
			if(pointF==0.)
				continue;
			double ratio = (fFull->f(E)-pointF)/pointF;
			ratioOut << E << "\t" << ratio << "\t" << aEmaxeV <<  "\t" << aEBLmodel <<
					"\t" << aSpecAlpha << "\t" << 0 << std::endl;
		}
		ratioOut.close();
	}
	std::cerr << "\n finished on " << time(0) - tStart << " sec" << std::endl;
}

double Evolution::EvolFactor(double z, int aEvolutionModel)
{
	switch (aEvolutionModel)
	{
	case EZDependenceOff:
		return 1.;
	case WaxmanBahcall://if m=0 corresponds to Waxman-Bahcall, hep-ph/9807282; Engel at al astro-ph/0101216 (oSFR in astro-ph/0605327v2)
		return WaxmanBahcallEvolution(z);
	case SFR03://if m=0 corresponds to Star Formation Rate from astro-ph/0309141
		return SfrEvolution(z);
	case SFR06://if m=0 corresponds to star formation rate from astro-ph/0607087 (lower figure 3)
		return SfrNewEvolution(z);
	case Weak://if m=0 corresponds to m=3 to z=1.8 and constant to z=3 going to zero there
		return weakEvolution(z);
	case Baseline://if m=0 corresponds to baseline evolution astro-ph-1512479
		return baselineEvolution(z);
	case Strong://if m=0 corresponds to fast evolution astro-ph-1512479
		return strongEvolution(z);
	case Engel9://if m=0 corresponds to Fig. 9 of http://lanl.arxiv.org/abs/astro-ph/0101216v2
		return Engel9Evolution(z);
	case SFR_Yuksel:
		return SFR_YukselEvolution(z);
	case GRB_Yuksel:
		return GRB_YukselEvolution(z);
	case AGN_Ahlers:
		return AGN_AhlersEvolution(z);
	default:
		Exception::Throw("invalid value of EvolutionModel parameter");
		break;
	};
	return 0.;
}

FermiTests::IntegralSpectrumOutput::IntegralSpectrumOutput(std::string aFile, double aEmin):
		fEmin(aEmin)
{
	fOut.open(aFile.c_str(),std::ios::out);
	if(!fOut.is_open())
		Exception::Throw("Failed to open " + aFile);
}

FermiTests::IntegralSpectrumOutput::~IntegralSpectrumOutput() {
	if(fOut.is_open())
		fOut.close();
}

void FermiTests::IntegralSpectrumOutput::AddParticle(const Particle& aParticle)
{
	if(aParticle.Energy>=fEmin && aParticle.Type == Photon)
		fParticles.insert(aParticle);
}

void FermiTests::IntegralSpectrumOutput::Flush(bool aIsFinal)
{
    if(!aIsFinal)
        return;//intermediate output is not supported
	if(!fOut.is_open())
	{
		std::cerr << "FermiTests::IntegralSpectrumOutput::Flush(): ignoring second call";
		return;
	}
	double E=-1;
	double sum=0.;
	double mult = units.Eunit*1e6;
	for(multiset<Particle,More>::iterator it = fParticles.begin(); it != fParticles.end(); it++)
	{
		const Particle& p = *it;
		if(p.Energy != E && sum > 0)
		{
			fOut << E*mult << "\t" << sum*mult << "\n";
			fE.insert(fE.begin(), E*mult);
			fIntFlux.insert(fIntFlux.begin(), sum*mult);
		}
		E = p.Energy;
		sum += p.Weight;
	}
	if(sum>0)
	{
		fOut << E*mult << "\t" << sum*mult << "\n";
		fE.insert(fE.begin(), E*mult);
		fIntFlux.insert(fIntFlux.begin(), sum*mult);
	}
	fOut.close();
}

TableFunction* FermiTests::IntegralSpectrumOutput::IntegralSpectrum() const
{
	ASSERT(fE.size()>0);//End must be called prior to this method call
	ASSERT(fE.size()>=2);//at least two data points should be present
	return new LinearFunc(fE, fIntFlux);
}

void TestList::RegisterTest(std::string aKey, TestProc aFunc) {
	fTests[aKey] = aFunc;
}

TestList::TestList() {
	TestList::RegisterTest("TurbulentMF",TurbulentMF::UnitTest);
	TestList::RegisterTest("MonochromaticMF",MonochromaticMF::UnitTest);
	TestList::RegisterTest("Deflection3D",Deflection3D::UnitTest);
	TestList::RegisterTest("NeutronDecay",NeutronDecay::UnitTestEconserv);
	TestList::RegisterTest("Stecker05EBL",Stecker2005EBL::UnitTest);
	TestList::RegisterTest("Stecker16EBL",Stecker16Background::UnitTest);
	TestList::RegisterTest("nucl_dis",PhotoDisintegration::UnitTest);
	TestList::RegisterTest("sampler",MathUtils::UnitTest);
}


EConservationTest::EConservationTest(double aMaxZ, double Emin_eV, u_long aRandSeed):
randomizer(aRandSeed),
result(Emin_eV*units.eV),
primary(Electron,0),
alphaThinning(0)
{
	if(!cosmology.IsInitialized())
		cosmology.Init();

	out = new TotalEnergyOutput();
	result.AddOutput(out);

	double cmbTemp = 2.73*units.K;
	backgr.AddComponent(new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp));
	backgr.AddComponent(new Backgrounds::Kneiske1001EBL());

	double logStepK = pow(10,0.05);
	double epsRel = 1e-3;

	stepZ = aMaxZ <0.2 ? aMaxZ /2 : 0.1;
	backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, aMaxZ, epsRel);
	pe = new PropagationEngine(particles, result, randomizer.CreateIndependent());
}

void EConservationTest::SetPrimary(const Particle& aPrimary, int aCount)
{
	ASSERT(nParticles>0);
	primary = aPrimary;
	for(int i=0; aCount > i; i++)
	{
		particles.AddPrimary(aPrimary);
	}
	nParticles = aCount;
}

void EConservationTest::Run()
{
	ASSERT(nParticles>0);//SetPrimary should be called prior to this method
	EnergyBasedThinning thinning(alphaThinning);
	pe->SetThinning(&thinning);
	pe->Run();

	double startE = nParticles* primary.Energy;
	double Ep = out->TotalE(primary.Type);
	double Esec = 0;
	for(int i=0; i<ParticleTypeEOF; i++)
		if(i!= primary.Type)
			Esec += out->TotalE((ParticleType)i);
	double deltaEp = startE-Ep;
	std::cout << "|deltaE_prim|/E_prim = " << deltaEp/startE << "\t\t(|deltaE_prim|-E_sec)/|deltaE_prim| = " << (deltaEp-Esec)/deltaEp <<  std::endl;
}


} /* namespace Test */
