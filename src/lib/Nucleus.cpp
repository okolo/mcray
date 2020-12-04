/*
 * Nucleus.cpp
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

//
// Photodisintegration cross sections as defined in F.W. Stecker, M.H. Salamon, Astrophys.J. 512 (1999) 521-526 (astro-ph/9808110v3)
//
#include "Nucleus.h"
#include "gsl/gsl_sf_erf.h"
#include "GammaPP.h"

#define AUTO -1

namespace Interactions {
	using namespace Utils;
	using namespace mcray;

	const double portionsHe4[] = {0.8, 0.2};
	static const TBranching branchingHe4 = {2, portionsHe4, AUTO, AUTO};

	const double portionsBe8[] = {0, 0, 0, 2.};
	// a trick to make Be9 decay to 2He4 + n only
	static const TBranching branchingBe8 = {4, portionsBe8, 0., 0.};

	const double portionsBe9[] = {0, 0, 0, 0, 2.};
	static const TBranching branchingBe9 = {5, portionsBe9, 0., 1.};

	const double portionsB10[] = {0.1, 0.3, 0.1, 0.1, 0.2, 0.2};
	static const TBranching branchingB10 = {6, portionsB10, AUTO, AUTO};//3.6

	const double portionsNa23[] = {0.10, 0.35, 0.10, 0.05, 0.15, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.018, 0.015,
								   0.012, 0.01};
	static const TBranching branchingNa23 = {15, portionsNa23, AUTO, AUTO};//4.349

	static const int EUnstable = -1;

#define NON 0
#define INF 1000
#define UNKNOWN_K 5
#define UNKNOWN_DELTA 5,25
#define FAST_DECAY {1, UNKNOWN_K, INF , 5,25}
#define NO_REACTION {2,15,0.,25,25}

// Photodisintegration cross sections as defined in Puget, J.L., F .W. Stecker & J.H. Bredekamp 1976, ApJ 205, 638
// with threshold energies from F.W. Stecker, M.H. Salamon, Astrophys.J. 512 (1999) 521-526 (astro-ph/9808110v3)
	static const TIsotope stableIsotopes[] = {
			{1,  0,  {NO_REACTION,                          NO_REACTION},    NON, NULL}, //  neutron
			{1,  1,  {NO_REACTION,                          NO_REACTION},    NON, NULL}, //  proton
			{2,  1,  {{2.2,  5,  0.97, 3,  15},             NO_REACTION},    NON, NULL}, //  H2
			{3,  2,  {{5.5,  13, 0.33, 18, 18}, /* NO_REACTION */{2,    15,   0.33, 13, 13}}, NON, NULL},//No threshold data for double nucleon decay of A03 in astro-ph/9808110v3
			// -> one might either use 2 MeV threshold as in Puget, J.L., F .W. Stecker & J.H. Bredekamp 1976, ApJ 205, 638
			// or just suppress the chanel. We choose the former option
			{4,  2,  {{19.8, 27, 0.47, 12, 12}, {26.1, 45,   0.11, 40, 40}}, 1.11, &branchingHe4},
			{5,  3,  {FAST_DECAY,                           NO_REACTION},    NON, NULL},
			{6,  3,  {{4.6, UNKNOWN_K, INF, UNKNOWN_DELTA}, FAST_DECAY},     NON, NULL},
			{7,  4,  {{7.3, UNKNOWN_K, INF, UNKNOWN_DELTA}, FAST_DECAY},     NON, NULL},
			{8,  4,  {NO_REACTION,                          NO_REACTION},    INF,  &branchingBe8},// a trick to make Be9 decay to 2He4 + n only
			{9,  4,  {{1.7,  26, 0.67, 20, 20}, {18.9, NON, NON, NON, NON}}, 1.00, &branchingBe9},
			{10, 5,  {{6.6,  25, 0.54, 11, 11}, {8.3,  25,   0.15, 11, 11}}, 1.03, &branchingB10},
			{11, 5,  {{11.2, 26, 0.85, 11, 11}, {18,   26,   0.15, 11, 11}}, 1.03, &branchingB10},
			{12, 6,  {{16,   23, 0.76, 6,  6},  {27.2, NON, NON, NON, NON}}, 1.06, &branchingB10},
			{13, 6,  {{4.9,  23, 0.71, 8,  8},  {20.9, 27,   0.05, 8,  8}},  1.06, &branchingB10},
			{14, 7,  {{7.6,  23, 0.46, 10, 10}, {12.5, 23,   0.37, 10, 10}}, 1.07, &branchingB10},
			{15, 7,  {{10.2, 23, 0.73, 10, 10}, {18.4, 23,   0.10, 10, 10}}, 1.07, &branchingB10},
			{16, 8,  {{12.1, 24, 0.83, 9,  9},  {22.3, 30,   0.04, 10, 10}}, 1.10, &branchingB10},
			{17, 8,  {{4.1,  24, 0.77, 9,  9},  {16.3, 29,   0.20, 10, 10}}, 1.10, &branchingB10},
			{18, 8,  {{8,    24, 0.67, 9,  9},  {21.8, 29,   0.20, 10, 10}}, 1.10, &branchingB10},
			{19, 9,  {{8,    23, 0.76, 14, 14}, {16,   29,   0.14, 14, 14}}, 1.10, &branchingB10},
			{20, 10, {{12.8, 22, 0.87, 12, 12}, {20.8, 26,   0.05, 8,  8}},  1.09, &branchingB10},
			{21, 10, {{6.8,  22, 0.84, 12, 12}, {19.6, 25,   0.08, 6,  6}},  1.09, &branchingB10},
			{22, 10, {{10.4, 22, 0.81, 12, 12}, {23.4, 21,   0.11, 4,  4}},  1.09, &branchingB10},
			{23, 11, {{8.8,  22, 0.83, 12, 12}, {19.2, 25,   0.12, 10, 10}}, 1.09, &branchingNa23},
			{24, 12, {{11.7, 19, 0.94, 11, 11}, {20.5, 29,   0.03, 6,  6}},  1.08, &branchingNa23},
			{25, 12, {{7.3,  23, 0.77, 9,  9},  {19,   28,   0.20, 7,  7}},  1.08, &branchingNa23},
			{26, 12, {{11.1, 18, 0.77, 8,  8},  {23.2, 26,   0.20, 8,  8}},  1.08, &branchingNa23},
			{27, 13, {{8.3,  21, 0.80, 8,  8},  {19.4, 29,   0.20, 12, 12}}, 1.05, &branchingNa23},
			{28, 14, {{11.6, 21, 1.01, 8,  8},  {19.9, 30,   0.02, 8,  8}},  1.04, &branchingNa23},
			{29, 14, {{8.5,  20, 0.83, 7,  7},  {20.1, 26,   0.20, 8,  8}},  1.04, &branchingNa23},
			{30, 14, {{10.6, 20, 0.83, 7,  7},  {22.9, 26,   0.20, 8,  8}},  1.04, &branchingNa23},
			{31, 15, {{7.3,  21, 0.85, 8,  8},  {17.9, 29,   0.20, 12, 12}}, 1.02, &branchingNa23},
			{32, 16, {{8.9,  22, 0.97, 12, 12}, {16.2, 30,   0.10, 12, 12}}, 1.00, &branchingNa23},
			{33, 16, {{8.6,  22, 0.82, 12, 12}, {17.5, 22,   0.25, 12, 12}}, 1.00, &branchingNa23},
			{34, 16, {{10.9, 22, 0.87, 12, 12}, {20.4, 22,   0.20, 12, 12}}, 1.00, &branchingNa23},
			{35, 17, {{6.4,  20, 0.87, 7,  7},  {17.3, 26,   0.22, 10, 10}}, 1.00, &branchingNa23},
			{36, 16, {{9.9,  22, 0.82, 12, 12}, {21.5, 22,   0.25, 12, 12}}, 1.00, &branchingNa23},
			{37, 17, {{8.4,  20, 0.81, 7,  7},  {18.3, 24,   0.28, 7,  7}},  1.00, &branchingNa23},
			{38, 18, {{10.2, 18, 0.86, 8,  8},  {18.6, 22,   0.24, 8,  8}},  0.98, &branchingNa23},
			{39, 19, {{6.4,  20, 0.73, 7,  7},  {16.6, 25,   0.38, 12, 12}}, 0.98, &branchingNa23},
			{40, 20, {{8.3,  20, 0.84, 6,  6},  {14.7, 26,   0.28, 10, 10}}, 0.96, &branchingNa23},
			{41, 20, {{7.8,  20, 0.92, 6,  6},  {17.7, 26,   0.20, 8,  8}},  0.96, &branchingNa23},
			{42, 20, {{10.3, 20, 1.02, 7,  7},  {18.1, 26,   0.10, 8,  8}},  0.96, &branchingNa23},
			{43, 20, {{7.9,  20, 0.97, 8,  8},  {18.2, 26,   0.15, 8,  8}},  0.96, &branchingNa23},
			{44, 20, {{11.1, 20, 0.92, 9,  9},  {21.6, 26,   0.20, 8,  8}},  0.96, &branchingNa23},
			{45, 21, {{6.9,  19, 0.97, 9,  9},  {18,   26,   0.15, 8,  8}},  0.95, &branchingNa23},
			{46, 22, {{10.3, 19, 1.03, 8,  8},  {17.2, 25,   0.10, 6,  6}},  0.95, &branchingNa23},
			{47, 22, {{8.9,  19, 1.03, 8,  8},  {18.7, 25,   0.10, 6,  6}},  0.95, &branchingNa23},
			{48, 22, {{9.9,  19, 1.03, 8,  8},  {24.2, 25,   0.10, 6,  6}},  0.95, &branchingNa23},
			{49, 22, {{8.1,  19, 1.03, 8,  8},  {19.6, 25,   0.10, 6,  6}},  0.95, &branchingNa23},
			{50, 22, {{10.9, 19, 1.03, 8,  8},  {21.8, 25,   0.10, 6,  6}},  0.95, &branchingNa23},
			{51, 23, {{8.1,  19, 1.02, 7,  7},  {19,   25,   0.11, 6,  6}},  0.95, &branchingNa23},
			{52, 24, {{10.5, 18, 1.08, 7,  7},  {18.6, 24,   0.05, 8,  8}},  0.95, &branchingNa23},
			{53, 24, {{7.9,  18, 1.03, 7,  7},  {18.4, 24,   0.10, 8,  8}},  0.95, &branchingNa23},
			{54, 24, {{9.7,  18, 0.93, 7,  7},  {20.9, 24,   0.20, 8,  8}},  0.95, &branchingNa23},
			{55, 25, {{8.1,  18, 0.93, 7,  7},  {17.8, 23.5, 0.20, 8,  8}},  0.95, &branchingNa23},
			{56, 26, {{10.2, 18, 0.98, 8,  8},  {18.3, 22,   0.15, 7,  7}},  0.95, &branchingNa23},
	};

	void CNucleus::getPhotopionMultipliers(int aA, double &aRn, double &aRp, double &aRN) {
		ASSERT(aA >= 0 && aA <= 56);
		const TIsotope &theIsotope = stableIsotopes[aA];
		aRN = pow((double) aA, -1. / 3.);//	= A^{2/3}(cross section)  / A(energy loss suppression)
		aRn = ((double) (aA - theIsotope.Z)) * aRN / ((double) aA);//  * N/A (fraction of neutrons)
		aRp = ((double) theIsotope.Z) * aRN / ((double) aA);//  * Z/A (fraction of protons)
	}

	double CNucleus::PPPmultiplierR(int aA) {
		const TIsotope &theIsotope = stableIsotopes[aA];
		double result = theIsotope.Z * theIsotope.Z;
		result /= (double) aA;
		return result;
	}

	int CNucleus::getZ(ParticleType aNucleus) {
		ASSERT(aNucleus >= StartNuclei - 2 && aNucleus < EndNuclei);
		int no = 2 + (int) aNucleus - StartNuclei;
		return stableIsotopes[no].Z;
	}

	int CNucleus::getA(ParticleType aNucleus) {
		ASSERT(aNucleus >= StartNuclei - 2 && aNucleus < EndNuclei);
		int no = 2 + (int) aNucleus - StartNuclei;
		return stableIsotopes[no].A;
	}

	CPhotoDisintegrationMap::CPhotoDisintegrationMap() {
		int i;
		fCanals.setIndexShift(-2);//
		const int NucleiCount = EndNuclei - StartNuclei;

		for (i = 0; i < NucleiCount + 2; i++) {
			fSecondaryCanals.add(new CPDMSecLine());
		}
		int curParticle = StartNuclei;
		for (i = 2; i < NucleiCount + 2; i++, curParticle++) {

			const TIsotope &isotope = stableIsotopes[i];
			CPDMLine *line = new CPDMLine;
			//if (IsEnabled(curParticle)) {
				if (isotope.pd[ESingle].sigma > 0.) {
					CPhotoDisintegrationCanal *c = new CPhotoDisintegrationCanal(isotope, ESingle);
					line->add(c);
					initSecMap(c);
				}
				if (isotope.pd[EDouble].sigma > 0.) {
					CPhotoDisintegrationCanal *c = new CPhotoDisintegrationCanal(isotope, EDouble);
					line->add(c);
					initSecMap(c);
				}
				if (isotope.pmn > 0) {
					CPhotoDisintegrationCanal *c = new CPhotoDisintegrationCanal(isotope, EMultiple);
					line->add(c);
					initSecMap(c);
				}
			//}
			fCanals.add(line);
		}
	}

	void CPhotoDisintegrationMap::initSecMap(CPhotoDisintegrationCanal *aCanal) {
		const TIsotope &isotope = aCanal->isotope();
		int maxDeltaA = aCanal->getMaxDeltaA();
		int A = isotope.A;
		ASSERT(maxDeltaA < A);
		double autoN = 0.;
		double autoP = 0.;
		for (int deltaA = 1; deltaA <= maxDeltaA; deltaA++) {
			const TIsotope &aProduct = stableIsotopes[A - deltaA];//in case A=2,3 secondary canal for proton (H1) may be duplicated, but it is ok, since summary rate is correct
			double rate = aCanal->getRate(deltaA);
			int deltaZ = isotope.Z - aProduct.Z;
			if (deltaZ < 0) {//processes with neutron emission followed by beta decay
				deltaZ = 0;
			}
			int deltaN = deltaA - deltaZ;
			autoN += deltaN * rate;
			autoP += deltaZ * rate;
			if (rate > 0.) {
				fSecondaryCanals[A - deltaA]->add(new TPhotoDisintegrationMapEntry(aCanal, rate));
			}
		}
		double rateP = aCanal->getRateP();
		double rateN = aCanal->getRateN();
		if (rateP == AUTO) {
			rateP = autoP;
		};
		if (rateN == AUTO) {
			rateN = autoN;
		};
		if (rateN > 0.) {
			fSecondaryCanals[0]->add(new TPhotoDisintegrationMapEntry(aCanal, rateN));
		}
		if (rateP > 0.) {
			fSecondaryCanals[1]->add(new TPhotoDisintegrationMapEntry(aCanal, rateP));
		}
	}
/*
//reset all entries
//used to recalculate R when background changes
	void CPhotoDisintegrationMap::recalculate(const CBackgroundIntegral *aBackground) {
		const int NucleiCount = EndNuclei - StartNuclei;
		for (int i = 2; i < NucleiCount + 2; i++) {
			CPDMLine &line = *(fCanals[i]);
			int jMax = line.size();
			for (int j = 0; j < jMax; j++)
				line[j]->recalculate(aBackground);
		}
		//printValidity(56);
	}

	void CPhotoDisintegrationMap::approximate(CPhotoDisintegrationMap &aMap1, CPhotoDisintegrationMap &aMap2,
											  double aPart) {
		const int NucleiCount = EndNuclei - StartNuclei;
		for (int i = 2; i < NucleiCount + 2; i++) {
			CPDMLine &line = *(fCanals[i]);
			CPDMLine &line1 = *(aMap1.fCanals[i]);
			CPDMLine &line2 = *(aMap2.fCanals[i]);
			int jMax = line.size();
			for (int j = 0; j < jMax; j++)
				line[j]->approximate(*(line1[j]), *(line2[j]), aPart);
		}
	}*/


// get list of canals contributing to nucleus A
// A=1 proton
// A=0 neutron
	const CPDMSecLine *CPhotoDisintegrationMap::getIncome(int A) const {
		return fSecondaryCanals[A];
	}

// get all suppression canals for nucleus A
	const CPDMLine *CPhotoDisintegrationMap::getOutcome(int A) const {
		return fCanals[A];
	}

	CPhotoDisintegrationCanal::CPhotoDisintegrationCanal(const TIsotope &aIsotope, TPhotoDisintegrationChanalType aType)
			:
			iIsotope(aIsotope),
			iType(aType),
			iMeanDeltaA(0) {
		double meanNucleonMassMeV = 0.5*(Particle::Mass(Proton)+Particle::Mass(Neutron));
		//const double alpha = 0.00729735307639648;//fine structure constant
		iSumTRK = 2. * M_PI * M_PI * EMInteraction::AlphaEM / meanNucleonMassMeV * units.Eunit;
		//double test = iSumTRK*Lunit*Lunit*Eunit*1e27;// (should be ~ 59.8 MeV * mb)
		iSumTRK *= ((double) ((iIsotope.A - iIsotope.Z) * iIsotope.Z)) / ((double) iIsotope.A);
		//iSumTRK *= 0.001;//test
		ASSERT(iSumTRK > 0);

		if (iType == EMultiple)
			iSigma = new TLinearSigma(iIsotope.pmn, iSumTRK);
		else
			iSigma = new TGaussianSigma(iIsotope.pd[iType], iSumTRK);
	}


	int    CPhotoDisintegrationCanal::getMaxDeltaA() {
		switch (iType) {
			case ESingle:
				return 1;
			case EDouble:
				return 2;
			case EMultiple:
				return iIsotope.branching->noOfModes;
		}
		return 0;//never goes here
	}

	double    CPhotoDisintegrationCanal::getMeanDeltaA() {
		switch (iType) {
			case ESingle:
				return 1;
			case EDouble:
				return 2;
			case EMultiple: {
				if (iMeanDeltaA) return iMeanDeltaA;
				int maxDeltaA = iIsotope.branching->noOfModes;
				for (int i = 1; i <= maxDeltaA; i++)
					iMeanDeltaA += iIsotope.branching->nucleiRates[i - 1] * ((double) i);
				return iMeanDeltaA;
			}
		}
		return 0;//never goes here
	}

	const double CPhotoDisintegrationCanal::AutoRate = AUTO;

	void CPhotoDisintegrationCanal::getEnergyLoss(std::vector<double> &aOutput) {
		/* //TODO implement if needed
		aOutput.copy(iR.ptr());
		aOutput *= (getMeanDeltaA() / ((double) iIsotope.A));
		 */
	}

	double    CPhotoDisintegrationCanal::getRate(int aDeltaA) {
		switch (iType) {
			case ESingle:
				return (aDeltaA == 1) ? 1. : 0.;
			case EDouble:
				return (aDeltaA == 2) ? 1. : 0.;
			case EMultiple:
				return (aDeltaA > 0 || aDeltaA <= iIsotope.branching->noOfModes) ? iIsotope.branching->nucleiRates[
						aDeltaA - 1] : 0.;
		}
		return 0.;//never goes here
	}

	double    CPhotoDisintegrationCanal::getRateP() {
		return (iType != EMultiple) ? AUTO : iIsotope.branching->pRate;
	}

	double    CPhotoDisintegrationCanal::getRateN() {
		return (iType != EMultiple) ? AUTO : iIsotope.branching->nRate;
	}

	/*
	void CPhotoDisintegrationCanal::recalculate(const CBackgroundIntegral *aBackground) {
		const int nn = iR.size();//Ranges().nE();
		Function *sigma = NULL;
		if (iType == EMultiple)
			sigma = new TLinearSigma(iIsotope.pmn);
		else
			sigma = new TGaussianSigma(iIsotope.pd[iType]);

		double minK = sigma->Xmin();
		double maxK = sigma->Xmax();

		for (int i = 0; i < nn; i++) {
			double R = 0;
			//TODO: implement if needed
			//double R = iSumTRK * aBackground->calculateR(Ranges().midNucleonGamma()[i], sigma, minK, maxK);
			ASSERT_VALID_NO(R);
			iR[i] = R;
		}
	}*/

	const double CPhotoDisintegrationCanal::TGaussianSigma::iDefaultMaxK = 30.;
	//30 MeV
	const double CPhotoDisintegrationCanal::TLinearSigma::iMaxK = 150.;//150MeV

//temporary
//double gsl_sf_erf(double x){return x;}

	CPhotoDisintegrationCanal::TGaussianSigma::TGaussianSigma(const TGaussCrossSection &aCS, double aMult, bool aCompatibilityMode) : iCS(aCS) {
		if(aCompatibilityMode) {//this version was originally used in propagation program
			iMaxK = iCS.k +
					0.5 *
					iCS.deltaPlus;//this adjustment of iMaxK is done to take into account wide deltas for light nuclei
			if (iMaxK <
				iDefaultMaxK)// the way it is done need to be checked, the only unclear reference is page 646 of ApJ 205 (Puget at al 1976)
				iMaxK = iDefaultMaxK;
		}
		else
			iMaxK = iDefaultMaxK;


		double leftIntegral =
				iCS.deltaMinus * sqrt(M_PI / 8.) * gsl_sf_erf(sqrt(2.) * (iCS.k - iCS.kTh) / iCS.deltaMinus);
		double rightIntegral = iCS.deltaPlus * sqrt(M_PI / 8.) * gsl_sf_erf(sqrt(2.) * (iMaxK - iCS.k) / iCS.deltaPlus);
		iNorm = iCS.sigma / (leftIntegral + rightIntegral)*aMult;
	}


	double CPhotoDisintegrationCanal::TGaussianSigma::f(double _x) const {
		double power = (_x * units.Eunit - iCS.k);
		power /= ((power < 0.) ? iCS.deltaMinus : iCS.deltaPlus);
		return iNorm * exp(-2. * power * power);
	}

	double CPhotoDisintegrationCanal::TGaussianSigma::Xmin() const {
		return iCS.kTh / units.Eunit;
	}

	double CPhotoDisintegrationCanal::TGaussianSigma::Xmax() const {
		return TGaussianSigma::iMaxK / units.Eunit;
	}

	CPhotoDisintegrationCanal::TLinearSigma::TLinearSigma(double aSigmaIntegral, double aMult) {
		iNorm = aMult*aSigmaIntegral / (iMaxK - TGaussianSigma::iDefaultMaxK);
	}

	double CPhotoDisintegrationCanal::TLinearSigma::Xmin() const {
		return CPhotoDisintegrationCanal::TGaussianSigma::iDefaultMaxK / units.Eunit;
	}

	double CPhotoDisintegrationCanal::TLinearSigma::Xmax() const {
		return TLinearSigma::iMaxK / units.Eunit;
	}
}//end of namespace Interactions