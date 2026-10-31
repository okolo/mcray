/*
 * ProtonPP.h
 *
 *  Created on: Oct 1, 2013
 *      Author: ok
 */

#ifndef PROTONPP_H_
#define PROTONPP_H_

#include "Interaction.h"
#include "GammaPP.h"

namespace Interactions {

using namespace MCPropag;

class ProtonPPSigma : public Utils::Function, EMInteraction
{
public:
	ProtonPPSigma();
	double f(double s) const;
	double sigmaLab(double k) const;
	double Xmin() const { return fThresholdS; }
protected:
	double fThresholdS;
	double fMpMe_1;
	double fMp2;
	double fCoef;
};

	/// This is outdated implementation of pair production by protons
	/// Use PPP class instead unless you are sure what you are doing
class ProtonPP : public RandomInteractionS{
public:
	ProtonPP(BackgroundIntegral* aBackground);
	virtual double Rate(const Particle& aParticle) const;
	virtual ~ProtonPP();
	virtual bool SampleS(const Particle& aParticle, double& aS, Randomizer& aRandomizer) const;
	virtual void SampleSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, double aS, Randomizer& aRandomizer) const;
	virtual RandomInteraction* Clone() const;
	static void UnitTest();
private:
	double SampleE(double aEmin, double aEmax, Randomizer& aRandomizer) const;
	SmartPtr<BackgroundIntegral>	fBackground;
	ProtonPPSigma 			fSigma;
	static double 					maxPPPemaxError;
};

///Used from PPP class to speed-up simulations
class ProtonPPcel : public CELInteraction
{
	/*!
	 * 2/E * Integrate[E_e * dSigma/dE_e,{E_e,Emin,Emax}]
	 * */
	class SigmaE : public ProtonPPSigma
	{
	public:
		SigmaE();
		double f(double s) const;
	private:
		static const double a[];
		static const double b[];
		double fInelCoef;
	};
public:
	static void UnitTest();
	ProtonPPcel(BackgroundIntegral* aBackground);
	virtual ~ProtonPPcel();
	double Rate(const Particle& aParticle) const;
	virtual void GetSecondaries(const Particle& aParticle, std::vector<Particle>& aSecondaries,
								cosmo_time aDeltaT, Randomizer& aRandomizer) const;

	CELInteraction* Clone() const;
private:
	SmartPtr<BackgroundIntegral>	fBackground;
	SigmaE							fSigma;
};

} /* namespace Interactions */
#endif /* PROTONPP_H_ */
