/*
 * Sophia.h
 *
 *  Created on: May 25, 2013
 *      Author: ok
 */

#ifndef SOPHIA_H_
#define SOPHIA_H_

#include "Particle.h"
#include "MathUtils.h"

namespace Interactions {
	using namespace MCPropag;

///stub for calls to SOPHIA
///the class is not thread-safe
	class SOPHIA {
		static int LastInit;
		static bool RandomSeedSet;

		static void Init(int aPrimary);

	public:
		static const int MaxNumberOfProducts;

		static int toSOPHIA(MCPropag::ParticleType aType);

		static MCPropag::ParticleType fromSOPHIA(int aType);

		//this method is not thread-safe!!!
		static void SamplePhotopion(MCPropag::ParticleType aPrimary, double aEnergyGeV, double aEpsilonGeV,
									double aThetaDeg, int &aNoSecondaries, double aSecEnergies[], int aSecTypes[]);

		//this method is not thread-safe!!!
		//sample energies of secondaries in units of primary energy in ultrarelativistic limit
		static void SamplePhotopionRel(MCPropag::ParticleType aPrimary, double aEpsPrimeGeV, int &aNoSecondaries,
									   double aSecEfrac[], int aSecTypes[]);

		//this method is not thread-safe
		static double CrossSection(int aPrimary, double aEpsPrimeGeV);

		static double Mass(MCPropag::ParticleType aParticle);

		static void SetRandomSeed(int aSeed);
	};

	class SophiaCS : public Utils::Function {
	public:
		SophiaCS(MCPropag::ParticleType aPrimary);

		virtual double f(double _x) const;

		virtual double Xmin() const;

	private:
		double fXmin;
		int fPrimary;
	};
}
#endif /* SOPHIA_H_ */
