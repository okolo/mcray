/*
 * Sophia.cpp
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

#include "Sophia.h"

extern "C"
   {
		void sample_photopion_(int& aPrimary, double& aEnergyGeV, double& aEpsilonGeV, double& aThetaDeg, int& aNoSecondaries, double aSecEnergies[], int aSecTypes[]);

		void sample_photopion_rel_(int& aPrimary, double& aEpsPrimeGeV, int& aNoSecondaries, double aSecEfrac[], int aSecTypes[]);


		//initialize SOPHIA
		//L0 - primary particle type (13 = proton, 14 = neutron)
		void initial_(int& L0);

		//double crossection_(double x, double NDIR, int NL0);

		void crossec_(int& L0, double& eps_prime, double& sigma);

		void get_mass_(int& L0,double& m);

		void set_random_seed_(int& aSeed);
   }

namespace Interactions {
	using namespace mcray;
	using namespace Utils;

	int SOPHIA::toSOPHIA(ParticleType aType) {
		switch (aType) {
			case Proton:
				return 13;
			case Neutron:
				return 14;
			case Photon:
				return 1;
			case Positron:
				return 2;
			case Electron:
				return 3;
			case NeutrinoE:
				return 15;
			case NeutrinoAE:
				return 16;
			case NeutrinoM:
				return 17;
			case NeutrinoAM:
				return 18;
			default:
				Exception::Throw("unsupported particle type");
				break;
		}
		return -1;
	}

	ParticleType SOPHIA::fromSOPHIA(int aType) {
		switch (aType) {
			case 1:
				return Photon;
			case 2:
				return Positron;
			case 3:
				return Electron;
			case 13:
				return Proton;
			case 14:
				return Neutron;
			case 15:
				return NeutrinoE;
			case 16:
				return NeutrinoAE;
			case 17:
				return NeutrinoM;
			case 18:
				return NeutrinoAM;
		}
		return ParticleTypeEOF;
	}

	void SOPHIA::Init(int aPrimary) {
		ASSERT(aPrimary == 13 || aPrimary == 14);
		if (aPrimary != LastInit) {
			initial_(aPrimary);
			LastInit = aPrimary;
		}
	}

	void SOPHIA::SamplePhotopion(ParticleType aPrimary, double aEnergyGeV, double aEpsilonGeV, double aThetaDeg,
								 int &aNoSecondaries, double aSecEnergies[], int aSecTypes[]) {
		ASSERT(aPrimary == Proton || aPrimary == Neutron);
		int primary = toSOPHIA(aPrimary);
		Init(primary);
		sample_photopion_(primary, aEnergyGeV, aEpsilonGeV, aThetaDeg, aNoSecondaries, aSecEnergies, aSecTypes);
		for (int i = 0; i < aNoSecondaries; i++)
			aSecTypes[i] = fromSOPHIA(aSecTypes[i]);
	}

	double SOPHIA::CrossSection(int aPrimary, double aEpsPrimeGeV) {
		Init(aPrimary);
		double result = 0;
		crossec_(aPrimary, aEpsPrimeGeV, result);
		return result;
	}

	const int SOPHIA::MaxNumberOfProducts = 2000;
	int SOPHIA::LastInit = -1;

	SophiaCS::SophiaCS(ParticleType aPrimary) :
			fPrimary(SOPHIA::toSOPHIA(aPrimary)) {
		const double Smin = 1.1646;
		ASSERT(fPrimary == 13 || fPrimary == 14);
		double pm = fPrimary == 13 ? 0.93827 : 0.93957;
		fXmin = 0.5 * (Smin / pm - pm);
	}

	double SophiaCS::f(double _x) const {
		return SOPHIA::CrossSection(fPrimary, _x);
	}

	double SophiaCS::Xmin() const {
		return fXmin;
	}

	void SOPHIA::SamplePhotopionRel(mcray::ParticleType aPrimary, double aEpsPrimeGeV, int &aNoSecondaries,
									double aSecEfrac[], int aSecTypes[]) {
		ASSERT(aPrimary == Proton || aPrimary == Neutron);
		int primary = toSOPHIA(aPrimary);
		Init(primary);
		sample_photopion_rel_(primary, aEpsPrimeGeV, aNoSecondaries, aSecEfrac, aSecTypes);
		for (int i = 0; i < aNoSecondaries; i++)
			aSecTypes[i] = fromSOPHIA(aSecTypes[i]);
	}

	bool SOPHIA::RandomSeedSet = false;

	void SOPHIA::SetRandomSeed(int aSeed) {
		bool doJob = true;
#pragma omp critical (SOPHIA)
		{
			if (RandomSeedSet)
				doJob = false;//initialize only once
			else
				RandomSeedSet = true;
		}
		if(doJob)
			set_random_seed_(aSeed);
	}

	double SOPHIA::Mass(mcray::ParticleType aParticle) {
		int particle = toSOPHIA(aParticle);
		double mGeV = 0.;
		get_mass_(particle, mGeV);
		return mGeV * units.GeV;
	}
}
