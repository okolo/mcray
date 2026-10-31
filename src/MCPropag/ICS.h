/*
 * ICS.h
 *
 *  Created on: Oct 22, 2011
 *      Author: ok
 */

#ifndef ICS_H_
#define ICS_H_
#include "Interaction.h"
#include "Background.h"
#include "GammaPP.h"

namespace Interactions {
using namespace MCPropag;

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
