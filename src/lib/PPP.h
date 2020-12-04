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

#ifndef PPP_H_
#define PPP_H_

#include "GammaPP.h"
#include "Interaction.h"

namespace Interactions {
using namespace mcray;

/// This class samples secondaries from pair production by protons
/// The energy loss by protons due to pair production is implemented
/// by class ProtonPPcel
class PPP: public RandomInteractionS, public EMInteraction {
	class SigmaTot : public Function
	{
	public:
		SigmaTot();
		double f(double s) const;
		double Xmin() const { return fThresholdS; }
		double Xmax() const { return fMaxS; }
	protected:
		const double fMp;
		const double fMe;
		const double fThresholdS;
		const double fMaxS;//introduced to avoid computational problems (with differential cross section)
							//for extremely large S (sqrt(S)>10 m_p,  where PPP is much less important than pion production anyway)
		const double fSigmaConst;
	};
	class DifSigma_r : public SigmaTot
	{
		class Kern : public Function
			{
			public:
				double f(double eps) const;
				void SetParams(double _k, double _e);
				double Xmin() const { return fEpsMin; }
				double Xmax() const { return k - 1.; }
			private:
				double fEpsMin;
				double k;
				double k2;//k^2
				double e;
				double e2;
				double e4;
			};
	public:
		double f(double r) const;
		void SetS(double aS);
		double SampleR(double aS, double aEpsRel, const Utils::MathUtils& aMath, double aRand);
		double Xmin() const { return f_rMin; }
		double Xmax() const { return f_rMax; }
	private:
		//Kern fKern;
		double fK;//photon energy in proton rest frame in units of electron mass
		double f_rMin;
		double f_rMax;
		Utils::MathUtils fMath;
	};

	class DifSigma_Eprime : public SigmaTot
	{// dSigma/dE_, where E_ electron gamma factor in proton rest frame
	public:
		DifSigma_Eprime();
		void SetK(double aK);//k-photon energy in proton rest frame in units of electron mass
		double f(double E) const;
		static double sigma(double aCosTheta, void* data);
		double Xmin() const { return 1.; }
		double Xmax() const { return fK-1.; }
	private:
		double fE;//electron gamma factor in proton rest frame
		double fK;//photon energy in proton rest frame in units of electron mass
		gsl_function fKern_Theta;
		Utils::MathUtils fMath;
		double fCoef;
	};

	class CosThetaSamplingEq : public Function
	{
		double fE;
		double fK;
	public:
		CosThetaSamplingEq(double aEprime, double aK):fE(aEprime),fK(aK){}
		double Xmin() const { return -1.; }
		double Xmax() const { return 1.; }
		double f(double aCosT) const;
	};

public:
	PPP(BackgroundIntegral* aBackground, double aScaleFactor = 1.);
	static void UnitTest();
	static void UnitTestSampling(double aE, double aZ);
	static void UnitTestEconserv(double aE, double aZ);
	static double DiffSigma(double aK, double aE, double aCosTheta);
	virtual ~PPP();
	RandomInteraction* Clone() const;
	virtual bool KillsPrimary() const
	    {
	    	return false;
	    }
	double Rate(const Particle& aParticle) const;
	bool SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const;
	void SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const;
	void SampleSecondariesMean(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const;
	void SampleSecondariesExact(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const;


private:
	SmartPtr<BackgroundIntegral>	fBackground;
	SigmaTot						fSigma;
	DifSigma_r						fDifSigma;
	DifSigma_Eprime					fDifSigma_Eprime;
	double f_ScaleFactor;//>=0 divide interaction rate by f_ScaleFactor and
	// multiply weight of secondary electrons and positrons by f_ScaleFactor
	double fMp;
	double fMe;
	Utils::MathUtils math;
	static double s_epsRel;//Desired relative error of integration
};

} /* namespace Interactions */

#endif /* PPP_H_ */
