/*
 * ICS.h
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
#ifndef ICS_H_
#define ICS_H_
#include "Interaction.h"
#include "Background.h"
#include "GammaPP.h"

namespace Interactions {
using namespace mcray;

/*!
 * Inverse Compton Scattering (ICS) with secondary photon energy higher than minimal is calculated in this class.
 * ICS part with soft secondary photons is approximated by continuous energy loss in class IcsCEL
 */
class ICS : public RandomInteractionS, EMInteraction{
	class Sigma : public Function
	{
	public:
		Sigma(double rMin);
		double f(double s) const;
		double Xmin() const {return fXmin;}
	private:
		double f_rMin;//minimal fraction of energy lost in single interaction
		double fXmin;
		double fM2;
	};
public:
	ICS(BackgroundIntegral* aBackground, double aGammaEmin);
	RandomInteraction* Clone() const;
	virtual ~ICS();
	double Rate(const Particle& aParticle) const;
	static void UnitTest(bool aPrintRates, bool aPropagate, bool aSplit, double aEmin, double aEmax, double aZmax);
	static void UnitTestSampling(double aEmin, double aE, double aZ);
	static void UnitTestDistribution(double aEmin, double aE, double aZ);
	static void UnitTestMono(bool aMono, int aAc, bool aRatesOnly);
	bool SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const;
	void SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const;

private:
	double SampleSecondaryGammaFracE(double aS, double rMin, double aRand) const;
	struct SecondaryEnergySamplingEquationPars
	{
		double s;
		double Integral;
	};
	/*!
	@param[in]  rMin ratio of minimal secondary photon energy to primary electron energy
	@param[in]  s center of mass energy in units of electron mass squared
	*/
	static double PartialSigma(double rMin, double s);

	/*!
	@param[in]  r ratio of lower energy secondary lepton energy to primary photon energy
	@param[in]  aParams pointer to SecondaryEnergySamplingEquationPars struct
	@return DifSigmaIntegral(r, aParams->s) - aParams->Integral
	*/
	static double SecondaryEnergySamplingEquation(double r, void* aParams);

	SmartPtr<BackgroundIntegral>	fBackground;
	double							fGammaEmin;//minimal energy of secondary photon
	bool							fDeleteElectron;//used for unit test
};

/*!
 * Continuous energy loss rate due to Inverse Compton Scattering (ICS) with soft secondary photons
 * (energy lower than minimal value provided in constructor) is calculated in this class
 */
class IcsCEL : public CELInteraction, EMInteraction
{
	/*!
	 * 1/E * Integrate[E_gamma * dSigma/dE_gamma,{E_gamma,0,E*f_rGammaMax}]
	 * */
	class SigmaE : public Function
	{
	public:
		SigmaE(double rGammaMax);
		double f(double s) const;
		double Xmin() const;
	private:
		double f_rGammaMax;//maximal fraction of energy lost in single interaction
		double f_m2;
	};
public:
	static void UnitTest();
	IcsCEL(BackgroundIntegral* aBackground, double aMaxGammaE);
	virtual ~IcsCEL();
	double Rate(const Particle& aParticle) const;
	CELInteraction* Clone() const;
private:
	SmartPtr<BackgroundIntegral>	fBackground;
	double							fMaxGammaE;//maximal energy of secondary photons
};

} /* namespace Interactions */

#endif /* ICS_H_ */
