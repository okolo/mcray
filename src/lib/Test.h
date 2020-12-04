/*
 * Test.h
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


#ifndef TEST_H_
#define TEST_H_

#include <string>
#include <map>
#include "Background.h"
#include "ParticleStack.h"
#include "PropagationEngine.h"

namespace Test {

using namespace mcray;
	typedef int (*TestProc)();

class TestList {
public:

	TestList();

	int Main(int argc, char** argv);
	static void FrameworkTest();
	static void BackgroundTest();
	static void BackgroundTest(const IBackground& aBackground, std::ostream& aOut);
	static void SplittedBackgroundTest(bool aSplitted);
	static void ThinningTest(double aAlpha);
	static void KachelriesTest(double aEmin_eV, double aEmax_eV, double aZmax,
			double aAlpha, bool aSplitted, int aNoParticles,
			bool aNoPP, bool aNoICS, bool aNoICScel, bool aNoRedshift,
			ParticleType aPrimary, const char* aOutputDir = 0);
	static void MonochromaticAndGausianBackgroundTest(bool aMonochromatic, bool aICS, bool aICScel, bool aPP);
	static void MultithreadTest(bool aMultithread);
	static void BichromaticTest(bool aGauss);
	static void SamplingTest();
	static void EngineTest(double aAlphaThinning);
	void RegisterTest(std::string aKey, TestProc aFunc);
private:
	std::map<std::string, TestProc> fTests;
};

class FermiTests {
	class IntegralSpectrumOutput : public TSmartReferencedObj<IOutput> {
	    struct More : public std::binary_function<Particle, Particle, bool>
	        {
	          bool
	          operator()(const Particle& __x, const Particle& __y) const
	          {
	        	  return __x.Energy > __y.Energy;
	          }
	        };
	    std::multiset<Particle,More> fParticles;//sort by energy in descending order
	public:
	    IntegralSpectrumOutput(std::string aFile, double aEmin);
		virtual ~IntegralSpectrumOutput();
		void AddParticle(const Particle& aParticle);
		void Flush(bool aIsFinal);
		TableFunction* IntegralSpectrum() const;//may only be called after End()
	private:
		std::ofstream	fOut;
		double fEmin;
		std::vector<double> fE;
		std::vector<double> fIntFlux;
	};
public:
	static void Main(int argc, char** argv);
private:
	static void InitCosmology();
	static void PrintUsage();
	static void AbsorptionTest(bool aInfiniteB, double aSpecAlpha, int aEBLmodel, double aEmax, double aDistanceMpc);
	static void DiffuseToPointTest(double aSpecAlpha, int aEBLmodel, double aEmax, int aEvolutionModel);
	static void InitSpec(int aParticlesCount, ParticleStack& aStack, double aSpecAlpha, double aEmin, double aEmax, double aZ, double aWeight, Result* aInitialOutput = 0);
	static Particle SamplePowerLawSpec(double aSpecAlpha, double aEmin, double aEmax, double aZ, double aWeight, double aRand);
	static IBackground* CreateEBL(int aEBLmodel);
	static void InitPropagEngine(PropagationEngine& aEngine, bool aInfiniteB, int aEBLmodel, IThinning* aThinning, double aZmax, double aEminCalc);
	static const double MinE;
};

class Evolution
{
public:
	enum TZDependence{
		EZDependenceOff = 0,//(1+z)^3
		EZDependencePowerLaw,//(1+z)^(3+m)
		EZDependencePowerTD,//t^{-4+p}

	//injection spectrum types below are proportional to (1+z)^m:

		WaxmanBahcall,//if m=0 corresponds to Waxman-Bahcall, hep-ph/9807282; Engel at al astro-ph/0101216 (oSFR in astro-ph/0605327v2)
		SFR03,//if m=0 corresponds to Star Formation Rate from astro-ph/0309141
		Weak,//if m=0 corresponds to m=3 to z=1.8 and constant to z=3 going to zero there
		Baseline,//if m=0 corresponds to baseline evolution astro-ph-1512479
		Strong,//if m=0 corresponds to fast evolution astro-ph-1512479
		Engel9,//if m=0 corresponds to Fig. 9 of http://lanl.arxiv.org/abs/astro-ph/0101216v2
		PowerCut,//(1+z)^m from Zmin up to PowerCutZ and const from PowerCutZ to Zmax
		SFR06,//star formation rate from astro-ph/0607087 (lower figure 3)
		SFR_Yuksel,//star formation rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 11)
		GRB_Yuksel,//GRB rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 12)
		AGN_Ahlers,//AGN rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 13)
		EZDependenceEOF
	};

	static double EvolFactor(double aZ, int aEvolutionModel);

	static inline double weakEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
	//m=3 to z=1.8 and constant to z=3 going to zero there
		return (aZ<1.8)?pow(1.+aZ,3.):((aZ<3.)?21.952:0.);
	}

	//baseline evolution astro-ph-0512479
	static inline double baselineEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
	//m=4 to z=1 and then constant to z=6 going to zero there.
		return (aZ<1.3)?pow(1.+aZ,3.1):((aZ<6)?13.2238005912547:0.);
	}

	//fast evolution astro-ph-0512479
	static inline double strongEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
	//m=4 to z=1 and then constant to z=6 going to zero there.
		return (aZ<1.0)?pow(1.+aZ,4.):((aZ<6)?16.:0.);
	}

	//AGN
	//Waxman-Bahcall, hep-ph/9807282; Engel at al astro-ph/0101216 (oSFR in astro-ph/0605327v2)
	static inline double WaxmanBahcallEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
		return (aZ<1.9)?pow(1.+aZ,3.):(24.389*((aZ<2.7)?1.:exp((2.7-aZ)/2.7)));
	}

	//Fig. 9 of http://lanl.arxiv.org/abs/astro-ph/0101216v2
	static inline double Engel9Evolution(double aZ){//corresponds to comoving frame only if m=0!!!
		return ((aZ<1.9)?pow(1.+aZ,4.):70.7281);
	}

	//Star Formation Rate from astro-ph/0309141
	static inline double SfrEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
		return (aZ<1.2)?pow(1.+aZ, 3.5):40.68052976*pow(1.+aZ, -1.2);//40.68052976==(1+1.2)^(3.5+1.2)
	}

	//star formation rate from astro-ph/0607087 (lower figure 3)
	static inline double SfrNewEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
		if(aZ<1.)
			return pow(1.+aZ, 4.65);
		else if(aZ<3.)
			return 7.7274906314*pow(1.+aZ, 1.7);
		else if(aZ<6.)
			return 7912.950406552335*pow(1.+aZ, -3.3);
		else if(aZ<8.)
			return 606880463687.06*pow(1.+aZ, -12.63);
		else
			return 0.;
	}

	//star formation rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 11)
	static inline double SFR_YukselEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
		if(aZ<1.)
			return pow(1.+aZ, 3.4);
		else if(aZ<4.)
			return 10.5560632861832*pow((1.+aZ)/2., -0.3);
		else
			return 8.01899573803639*pow((1.+aZ)/5., -3.5);
	}

	//GRB rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 12)
	static inline double GRB_YukselEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
		if(aZ<1.)
			return pow(1.+aZ, 4.8);
		else if(aZ<4.)
			return 27.857618025476*pow((1.+aZ)/2., 1.1);
		else
			return 76.3269641062938*pow((1.+aZ)/5., -2.1);
	}

	//AGN rate from http://lanl.arxiv.org/abs/1103.3574v2  (formula 13)
	static inline double AGN_AhlersEvolution(double aZ){//corresponds to comoving frame only if m=0!!!
		if(aZ<1.7)
			return pow(1.+aZ, 5);
		else if(aZ<2.7)
			return 143.48907;
		else
			return 143.48907*pow(10.,2.7-aZ);
	}
};

	class EConservationTest{
	public:
		EConservationTest(double aMaxZ, double Emin_eV=1e6, u_long aRandSeed = 0);
		void SetPrimary(const Particle& aPrimary, int aCount);
		void Run();
		Randomizer& GetRandomizer() { return randomizer;}
		PropagationEngine& Engine(){return *pe;}
		BackgroundIntegral* Backgr(){ return backgrI;}
		double alphaThinning;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	private:
		Randomizer randomizer;
		ParticleStack particles;
		Result result;
		SmartPtr<TotalEnergyOutput> out;
		CompoundBackground backgr;
		double stepZ;
		std::vector<PInteraction>     fInteractions;
		PBackgroundIntegral backgrI;
		SafePtr<PropagationEngine> pe;
		Particle primary;
		int nParticles;
	};
} /* namespace Test */
#endif /* TEST_H_ */
