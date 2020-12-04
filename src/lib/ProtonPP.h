/*
 * ProtonPP.h
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


#ifndef PROTONPP_H_
#define PROTONPP_H_

#include "Interaction.h"
#include "GammaPP.h"

namespace Interactions {

using namespace mcray;

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
