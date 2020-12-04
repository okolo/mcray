/*
 * GammaPP.h
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

#ifndef GAMMAPP_H_
#define GAMMAPP_H_

#include "Interaction.h"
#include "Background.h"

namespace Interactions {
using namespace mcray;

class EMInteraction
{
public:
	static const double AlphaEM;//fine structure constant
	static const double SigmaCoef;//AlphaEM*AlphaEM*2*M_PI used by PP and ICS
};

class GammaPP : public RandomInteractionS, EMInteraction{
	class Sigma : public Function
	{
	public:
		Sigma();
		double f(double s) const;
		double Xmin() const { return fThresholdS; }
	private:
		double fThresholdS;
	};
public:
	GammaPP(BackgroundIntegral* aBackground);
	RandomInteraction* Clone() const;
	static void UnitTest();
	static void UnitTest(bool aPrintRates, bool aPropagate, bool aSplit, bool aBestFitEBL);
	static void UnitTestSampling(double aE, double aZ);
	virtual ~GammaPP();
	double Rate(const Particle& aParticle) const;
	bool SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const;
	void SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const;

private:
	double SampleSecondaryFracE(double aS, double aRand) const;
	struct SecondaryEnergySamplingEquationPars
	{
		double s;
		double Integral;
	};
	/*!
	@param[in]  s center of mass energy in units of electron mass squared. s>4
	@param[in]  r ratio of lower energy secondary lepton energy to primary photon energy  1/2*(1 - (1 - 4/s)^(1/2)) < r < 1/2
	*/
	static double DifSigmaIntegral(double r, double s);
	static double DifSigmaIntegralUnnorm(double r, double s);
	static double DifSigmaIntegralTest(double r, double s, bool aSeries, double& a);

	/*!
	@param[in]  r ratio of lower energy secondary lepton energy to primary photon energy
	@param[in]  aParams pointer to SecondaryEnergySamplingEquationPars struct
	@return DifSigmaIntegral(r, aParams->s) - aParams->Integral
	*/
	static double SecondaryEnergySamplingEquation(double r, void* aParams);

	SmartPtr<BackgroundIntegral>	fBackground;
	Sigma							fSigma;
	static bool						fDebugOutput;
};

} /* namespace Interactions */

#endif /* GAMMAPP_H_ */
