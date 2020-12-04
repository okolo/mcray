/*
 * Nucleus.h
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


#if !defined(NUCLEUS_H__INCLUDED_)
#define NUCLEUS_H__INCLUDED_


#include "Particle.h"
#include "MathUtils.h"

namespace Interactions {
	using namespace Utils;
	using namespace mcray;

	class CBackgroundIntegral;

	struct TGaussCrossSection {
		double kTh;
		//threshold energy (MeV)
		double k; //middle energy (MeV)
		double sigma; //total sigma integrated
		double deltaMinus; //energy range (MeV)
		double deltaPlus; //energy range (MeV)
	};


	struct TBranching {
		int noOfModes;
		const double *nucleiRates;
		double pRate;
		double nRate;
	};

	struct TIsotope {
		int A; //mass
		int Z; //electric charge

		// photodisintegration with single and double nucleon emission
		TGaussCrossSection pd[2];

		// photodisintegration with multiple nucleon emission
		double pmn;

		// branching rates for multiple nucleon emission
		const TBranching *branching;
	};

	enum TPhotoDisintegrationChanalType {
		ESingle = 0, //order and numbers are important (see for ex. CPhotoDisintegrationCanal::recalculate)
				EDouble,
		EMultiple
	};

	class CPhotoDisintegrationCanal {
	public:
		CPhotoDisintegrationCanal(const TIsotope &aIsotope, TPhotoDisintegrationChanalType aType);

		Function* Sigma(){return iSigma;}

		//forces recalculation of R next time R(int) is called
		//used to recalculate R when background changes
		//void recalculate(const CBackgroundIntegral *aBackground);

		//void approximate(CPhotoDisintegrationCanal &aCanal1, CPhotoDisintegrationCanal &aCanal2, double aPart);

		int getMaxDeltaA();

		double getMeanDeltaA();

		//needed for energy loss outputs
		double getRate(int aDeltaA);

		double getRateP();

		double getRateN();

		void getEnergyLoss(std::vector<double> &aOutput);

		inline ParticleType getParticle() const { return (ParticleType)(StartNuclei + iIsotope.A - 2); };

		inline TPhotoDisintegrationChanalType type() { return iType; };

		inline const TIsotope &isotope() { return iIsotope; };

		static const double AutoRate;


		class TGaussianSigma : public Function {
		public:
			TGaussianSigma(const TGaussCrossSection &aCS, double aMult, bool aCompatibilityMode=false);

			double f(double _x) const;

			//sigma devided by TRK sum
			double Xmin() const;

			double Xmax() const;

		private:
			double iNorm;
			const TGaussCrossSection &iCS;
			double iMaxK;
		public:
			static const double iDefaultMaxK;//30MeV
		};

		class TLinearSigma : public Function {
		public:
			TLinearSigma(double aSigmaIntegral, double aMult);

			double f(double _x) const { return iNorm; };

			//sigma devided by TRK sum
			double Xmin() const;

			double Xmax() const;

		private:
			double iNorm;
		public:
			static const double iMaxK;//150MeV
		};

		const TIsotope &iIsotope;
		TPhotoDisintegrationChanalType iType;
		double iSumTRK;
		double iMeanDeltaA;
		SafePtr<Function> iSigma;
	};

	struct TPhotoDisintegrationMapEntry {
		CPhotoDisintegrationCanal *iCanal;
		double iRate;

		TPhotoDisintegrationMapEntry(CPhotoDisintegrationCanal *aCanal, double aRate) : iCanal(aCanal),
																						iRate(aRate) { };
	};

	typedef AutoDeletePtrArray <TPhotoDisintegrationMapEntry> CPDMSecLine;
	typedef AutoDeletePtrArray <CPhotoDisintegrationCanal> CPDMLine;//canals for single isotope


	class CPhotoDisintegrationMap {
	public:
		CPhotoDisintegrationMap();

		//reset all entries
		//used to recalculate R when background changes
		//void recalculate(const CBackgroundIntegral *aBackground);

		//void approximate(CPhotoDisintegrationMap &aMap1, CPhotoDisintegrationMap &aMap2, double aPart);

		// get list of canals contributing to nucleus A
		// A=1 proton
		// A=0 neutron
		const CPDMSecLine *getIncome(int A) const;

		// get all suppression canals for nucleus A
		const CPDMLine *getOutcome(int A) const;

		void printValidity(int aA);

//debug
		void printEnergyLoss(int aA, std::string &aFileName);

	private:

		void initSecMap(CPhotoDisintegrationCanal *aCanal);

		typedef AutoDeletePtrArray <CPDMLine> CPDMTable;
		typedef AutoDeletePtrArray <CPDMSecLine> CPDMSecTable;

		CPDMTable fCanals;
		CPDMSecTable fSecondaryCanals;
	};

	class CNucleus {
	public:
		static double PPPmultiplierR(int aA);

		static int getZ(ParticleType aNucleus);

		static int getA(ParticleType aNucleus);

		static void getPhotopionMultipliers(int aA, double &aRn, double &aRp, double &aRN);

	};
}

#endif // !defined(NUCLEUS_H__INCLUDED_)
