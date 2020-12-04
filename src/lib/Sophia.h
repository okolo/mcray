/*
 * Sophia.h
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

#ifndef SOPHIA_H_
#define SOPHIA_H_

#include "Particle.h"
#include "MathUtils.h"

namespace Interactions {
	using namespace mcray;

///stub for calls to SOPHIA
///the class is not thread-safe
	class SOPHIA {
		static int LastInit;
		static bool RandomSeedSet;

		static void Init(int aPrimary);

	public:
		static const int MaxNumberOfProducts;

		static int toSOPHIA(mcray::ParticleType aType);

		static mcray::ParticleType fromSOPHIA(int aType);

		//this method is not thread-safe!!!
		static void SamplePhotopion(mcray::ParticleType aPrimary, double aEnergyGeV, double aEpsilonGeV,
									double aThetaDeg, int &aNoSecondaries, double aSecEnergies[], int aSecTypes[]);

		//this method is not thread-safe!!!
		//sample energies of secondaries in units of primary energy in ultrarelativistic limit
		static void SamplePhotopionRel(mcray::ParticleType aPrimary, double aEpsPrimeGeV, int &aNoSecondaries,
									   double aSecEfrac[], int aSecTypes[]);

		//this method is not thread-safe
		static double CrossSection(int aPrimary, double aEpsPrimeGeV);

		static double Mass(mcray::ParticleType aParticle);

		static void SetRandomSeed(int aSeed);
	};

	class SophiaCS : public Utils::Function {
	public:
		SophiaCS(mcray::ParticleType aPrimary);

		virtual double f(double _x) const;

		virtual double Xmin() const;

	private:
		double fXmin;
		int fPrimary;
	};
}
#endif /* SOPHIA_H_ */
