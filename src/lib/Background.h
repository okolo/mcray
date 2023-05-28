/*
 * Background.h
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

#ifndef BACKGROUND_H_
#define BACKGROUND_H_

#include "TableFunction.h"
#include <limits>
#include "Particle.h"
#include "Randomizer.h"

namespace mcray {

using namespace Utils;

class IBackground {
public:
	/*!
	 Calculate density E*dn(E)/dE at given redshift in internal units. Physical density should be returned, not a comoving one
	 @param[in]  aE  background photon energy (in internal units)
	 @param[in]  aZ  red shift
	 @return background photon density E*dn(E)/dE in internal units
	 */
	virtual double n(double aE, double aZ) const = 0;
	virtual double MinE() const = 0;
	virtual double MaxE() const = 0;
	virtual double MinZ() const = 0;
	virtual double MaxZ() const = 0;
	virtual std::string Name() const = 0;
	virtual ~IBackground(){};
};

class BackgroundUtils
{
public:
	static double CalcIntegralDensity(const IBackground& aBackground, double aRelError=1e-3);
	static void Print(const IBackground* aBackground, double aZ, std::ostream& aOut, int nIntervals = 100);
	static void UnitTest(const IBackground& aBackground, int aPointsZ=5, int aPointsE = 100);
private:
	static double DensityCalculator (double x, void * params);
};

class TableBackground : public IBackground{
public:
	TableBackground(std::string aName);

	virtual ~TableBackground();

	double n(double aE, double aZ) const;
	double MinE() const;
	double MaxE() const;
	double MinZ() const;
	double MaxZ() const;
	std::string Name() const { return fName; }
protected:
	/*!
	@param[in]  aDir directory containing data files
	@param[in]  aFileList list of file names convertable to double (represent value of z). Last element in the array should be NULL, the files should be ordered increasingly.
	@param[in]  aIsLogscaleY if true left and right values of functions in fData will be set to -Infinity, otherwise they will be set to 0
	*/
	void Init(std::string aDir, const char** aFileList, bool aIsLogscaleY);
	//convert energy in internal units to table x scale
	virtual double EtoRecordX(double aE, double aZ) const = 0;
	//retrieve energy in internal units from the data record
	virtual double recordToE(double aX, double aZ) const = 0;
	////retrieve concentration E*dn/dE in internal units from the data record
	virtual double recordToN(double aX, double aY, double aZ) const = 0;
private:
	std::vector<GSLTableFunc*>	fData;
	Vector						fZ;
	Vector						fZindex;
	GSLTableFunc*				fZtoIndex;
	double						fMinE;
	double						fMaxE;
	std::string					fName;
};

/// This class was originally designed to read slightly modified Elmag background tables
/// (see script ElmagTest/Tables/makeTablesFormcray.sh)
/// The table file is supposed to use the following format
/// 	 z1				z2				z3				...			zMax
/// E1	 dN/dE(z1,E1)	dN/dE(z2,E1)	dN/dE(z3,E1)	...			dN/dE(zMax,E1)
/// E2	 dN/dE(z1,E2)	dN/dE(z2,E2)	dN/dE(z3,E2)	...			dN/dE(zMax,E2)
/// ................
/// Emax dN/dE(z1,Emax)	dN/dE(z2,Emax)	dN/dE(z3,Emax)	...			dN/dE(zMax,Emax)
///
/// where [E]=eV, [dN/dE] = 1/cm^3/eV
/// Linear interpolation in z, log(E) and log (N) is performed as in Elmag
/// By default for all z<z_min dN/dE(z) == dN/dE(zMin)
class MatrixBackground : public IBackground{
public:
	MatrixBackground(std::string aName, std::string aTableFile, bool aIsComoving, bool aExtendToZero);

	virtual ~MatrixBackground();

	double n(double aE, double aZ) const;
	double MinE() const { return fMinE; }
	double MaxE() const { return fMaxE; }
	double MinZ() const { return fExtendToZero ? 0 : fZ[0]; } //For all z<z_min dN/dE(z) == dN/dE(zMin) if aExtendToZero == true
	double MaxZ() const { return *(fZ.end()-1); }
	std::string Name() const { return fName; }
private:
	std::vector<Vector*>		fLogEdN_dE;
	Vector						fZ;
	Vector						fZindex;
	GSLTableFunc*				fZtoIndex;
	Vector						fLogE;
	Vector						fLogEindex;
	GSLTableFunc*				fLogEtoIndex;
	double						fMinE;
	double						fMaxE;
	std::string					fName;
	bool						fExtendToZero;
	bool						fIsComoving;
	char*						fBuffer;
};

class PlankBackground : public IBackground{
public:
	/*!
	@param[in]  aT0 background temperature (internal units)
	@param[in]  aEmin0 minimal background energy (internal units)
	@param[in]  aEmax maximal background energy (internal units)
	@param[in]  aMinZ minimal background red shift
	@param[in]  aMaxZ maximal background red shift
	*/
	PlankBackground(double aT0, double aEmin, double aEmax, double aMinZ=0, double aMaxZ=10000, bool aFixedTemperature = false);
	double n(double aE, double aZ) const;
	double MinE() const { return fMinE; }
	double MaxE() const { return fMaxE; }
	double MinZ() const { return fMinZ; }
	double MaxZ() const { return fMaxZ; }
	std::string Name() const { return fName; }
	virtual ~PlankBackground(){};
private:
	double fT0;
	double fMinE;
	double fMaxE;
	double fMinZ;
	double fMaxZ;
	std::string fName;
	bool	fFixedTemperature;
};

class GaussianBackground : public IBackground{
public:
	/*!
	@param[in]  aE background photons energy central value (internal units)
	@param[in]  aSigma Gaussian distribution sigma (internal units)
	@param[in]  aWidth width at which the distribution is cut (in units of sigma)
	@param[in]  aConcentration integral concentration of photons (internal units)
	@param[in]  aMinZ minimal background red shift
	@param[in]  aMaxZ maximal background red shift
	*/
	GaussianBackground(double aE, double aSigma, double aWidth, double aConcentration, double aMinZ=0, double aMaxZ=10000);
	double n(double aE, double aZ) const;
	double MinE() const { return fMinE; }
	double MaxE() const { return fMaxE; }
	double MinZ() const { return fMinZ; }
	double MaxZ() const { return fMaxZ; }
	std::string Name() const { return fName; }
	virtual ~GaussianBackground(){};
private:
	double fE0;
	double fSigma;
	double fMinE;
	double fMaxE;
	double fMinZ;
	double fMaxZ;
	double fNorm;
	double fWidth;
	std::string fName;
};

class CompoundBackground : public IBackground{
public:
	//argument ownership is transfered to CompoundBackground object
	void AddComponent(IBackground* aBackground, double aWeight=1.);
	double n(double aE, double aZ) const;
	double MinE() const;
	double MaxE() const;
	double MinZ() const;
	double MaxZ() const;
	std::string Name() const;
	virtual ~CompoundBackground();
private:
	std::vector<IBackground*> fBackgrounds;
	std::vector<double>		  fWeights;
};

class BackgroundIntegral : public SmartReferencedObj
{
public:
    virtual Function* GetRateDistributionS(const Function& aSigmaS, const Particle& aParticle, double& aSmin, double& aSmax) = 0;
	///calculate rates using sigma(s) where s - total energy squared in center of mass frame
	virtual double GetRateAndSampleS(const Function& aSigmaS, const Particle& aParticle, Randomizer& aRandomizer, double& aS, double aAbsError=0) = 0;//todo: make use of aAbsError parameter in propagation engine
	virtual double GetRateS(const Function &aSigmaS, const Particle &aParticle, double aAbsError = 0) = 0;//todo: make use of aAbsError parameter in propagation engine

	///calculate rates using sigma(k) where k - background photon energy in the massive particle rest frame
	virtual double GetRateAndSampleK(const Function& aSigmaK, const Particle& aParticle, Randomizer& aRandomizer, double& aK, double aAbsError = 0) = 0;
	virtual double GetRateK(const Function &aSigmaK, const Particle &aParticle, double aAbsError = 0) = 0;
	virtual BackgroundIntegral* Clone() const = 0;
};

class MonochromaticBackgroundIntegral : public BackgroundIntegral
{
	class RateSKern : public Function
	{
	public:
		RateSKern(const Function& aSigmaS, const Particle& aParticle);

		//(s-m^2) sigma(s)
		double f(double _x) const;
	private:
		const Function&	fSigma;
		double				fM2;
	};
	class RateKKern : public Function
	{
	public:
		RateKKern(const Function& aSigmaK):fSigma(aSigmaK){};
		//k sigma(k)
		double f(double _x) const;
	private:
		const Function&	fSigma;
	};
public:
	/*!
	@param[in]  aConcentration Z-dependence of background particles concentration (internal units)
	@param[in]  aEnergy Z-dependence of background particles energy (internal units)
	@param[in]  aLogStepE energy binning used to calculate rates and sample S
	@param[in]  aEpsRel maximal relative error criterion (used in rate calculation)
	*/
	MonochromaticBackgroundIntegral(Function& aConcentration, Function& aEnergy, double aLogStepE, double aEpsRel):
		fConcentration(aConcentration),
		fEnergy(aEnergy),
		fLogStepE(aLogStepE),
		fEpsRel(aEpsRel){};
	MonochromaticBackgroundIntegral(const MonochromaticBackgroundIntegral& aBackgroundIntegral);
    Function* GetRateDistributionS(const Function& aSigmaS, const Particle& aParticle, double& aSmin, double& aSmax);

    double GetRateAndSampleS(const Function& aSigmaS, const Particle& aParticle, Randomizer& aRandomizer, double& aS, double aAbsError=0){
        return GetRateAndSampleS(aSigmaS, aParticle, aRandomizer.Rand(), aS, aAbsError);
    }

	double GetRateS(const Function& aSigma, const Particle& aParticle, double aAbsError=0);

	double GetRateAndSampleK(const Function& aSigmaK, const Particle& aParticle, Randomizer& aRandomizer, double& aK, double aAbsError = 0);
	double GetRateK(const Function &aSigmaK, const Particle &aParticle, double aAbsError = 0);

	BackgroundIntegral* Clone() const;
	virtual ~MonochromaticBackgroundIntegral(){};
private:
    double GetRateAndSampleS(const Function& aSigmaS, const Particle& aParticle, double aRand, double& aS, double aAbsError=0);
	MathUtils		math;
	bool			fEnableDebugOutput;
	Function&		fConcentration;
	Function&		fEnergy;
	double			fLogStepE;
	const double 	fEpsRel;
};

/// calculates background integral function I(k,z)=Integrate[x^-2 * dn(x,z)/dx, {x, k, Infinity}]
/// calculates Interaction rate Omega=1/8E^2/beta * Integrate[(s-m^2) sigma(s) I((s-m^2)/2/E/(1+beta),z), {s,s_min,Infinity}]
/// The class is not thread safe. Copy should be created for each thread using copy constructor;
class ContinuousBackgroundIntegral : public BackgroundIntegral{
	class Kern : public Function
	{
	public:
		Kern(IBackground&	aBackground) : fBackground(aBackground){}
		//x^-2 * dn(x,z)/dx === x^-3 x dn/dx
		virtual double f(double _x) const;
		double			Z;
	private:
		const IBackground&	fBackground;
	};
	class RateSKern : public Function
	{
	public:
		RateSKern(const Function&	aBackrgoundIntegral, const Function& aSigma, const Particle& aParticle);

		//(s-m^2) sigma(s) I((s-m^2)/2/E/(1+beta),z)
		double f(double _x) const;
	private:
		const Function&		fBackrgoundIntegral;
	protected:
		const Function&		fSigma;
		double				fM2;
		// 1/2/E/(1+beta)
		double				fMult;
		// 1/2/E/(1-beta) or -1 if beta=1
		double				fMaxEmult;
	};
	class RateSKernZ : public RateSKern
	{
	public:
		RateSKernZ(const Function&	aBackrgoundIntegralZ1, const Function&	aBackrgoundIntegralZ2, double aZ1, double aZ2, const Function& aSigma, const Particle& aParticle);

		//(s-m^2) sigma(s) I((s-m^2)/2/E/(1+beta),z)
		double f(double _x) const;
	private:
		const Function&		fBackrgoundIntegralZ1;
		const Function&		fBackrgoundIntegralZ2;

		double				fCoef1;
		double				fCoef2;
	};
	class RateKKern : public Function
	{
	public:
		RateKKern(const Function&	aBackrgoundIntegral, const Function& aSigma, const Particle& aParticle);
		double f(double _x) const;
	private:
		const Function&		fBackrgoundIntegral;
	protected:
		const Function&		fSigma;

		// m/E/(1+beta)
		double				fMult;
		// m/E/(1-beta) or -1 if beta=1
		double				fMaxEmult;
	};
	class RateKKernZ : public RateKKern
	{
	public:
		RateKKernZ(const Function&	aBackrgoundIntegralZ1, const Function&	aBackrgoundIntegralZ2, double aZ1, double aZ2, const Function& aSigma, const Particle& aParticle);
		double f(double _x) const;
	private:
		const Function&		fBackrgoundIntegralZ1;
		const Function&		fBackrgoundIntegralZ2;

		double				fCoef1;
		double				fCoef2;
	};
public:
	ContinuousBackgroundIntegral(IBackground& aBackground, double aStepZ, double aLogStepK, double aZmaxLimit, double aEpsRel = 0);
	BackgroundIntegral* Clone() const;
	double GetRateAndSampleS(const Function& aSigmaS, const Particle& aParticle, Randomizer& aRandomizer, double& aS, double aAbsError = 0){
        return GetRateAndSampleS(aSigmaS, aParticle, aRandomizer.Rand(), aS, aAbsError);
    }
    Function* GetRateDistributionS(const Function& aSigmaS, const Particle& aParticle, double& aSmin, double& aSmax);
	double GetRateS(const Function &aSigmaS, const Particle &aParticle, double aAbsError = 0);

	double GetRateAndSampleK(const Function& aSigmaK, const Particle& aParticle, Randomizer& aRandomizer, double& aK, double aAbsError = 0){
		return GetRateAndSampleK(aSigmaK, aParticle, aRandomizer.Rand(), aK, aAbsError);
	}
	double GetRateK(const Function &aSigmaK, const Particle &aParticle, double aAbsError = 0);

	virtual ~ContinuousBackgroundIntegral();
	static void UnitTest();
private:
    double GetRateAndSampleS(const Function& aSigmaS, const Particle& aParticle, double aRand, double& aS, double aAbsError = 0);
	double GetRateAndSampleK(const Function& aSigmaK, const Particle& aParticle, double aRand, double&aK, double aAbsError = 0);
	ContinuousBackgroundIntegral(const ContinuousBackgroundIntegral& aBackgroundIntegral);
	//used in unit test
	static double PlankBackgroundIntegral(double aT, double aE);
	//used in unit test
	static double PlankBackgroundTestRate(double aT, double aE, double aS);

	double I(double k, double z);

	double			fStepZ;
	double			fLogStepK;
	double			fZmaxLimit;
	double			fMaxK;
	double			fMinK;
	const double 	fEpsRel;
	std::vector<Function*> fIntegrals;//integrals I for different values of z
	std::vector<Vector*> fData;
	Vector			fBinningK;
	MathUtils		math;
	bool			fEnableDebugOutput;
};

class HighRedshiftBackgrExtension : public IBackground
{
public:
	HighRedshiftBackgrExtension(IBackground* aBackground, double aDeltaZconst, double aPowerLow, double aDeltaZexp, double aZmax=1000);
	virtual double n(double aE, double aZ) const;
	virtual double MinE() const;
	virtual double MaxE() const;
	virtual double MinZ() const;
	virtual double MaxZ() const;
	virtual std::string Name() const;
	virtual ~HighRedshiftBackgrExtension();
private:
	IBackground* fBackground;
	double fZconst;
	double fPowerLow;
	double fDeltaZexp;
	double fZmax;
	double fInnerZmax;
};

class CuttedBackground : public IBackground
{
public:
	//takes ownership of aBackgr
	CuttedBackground(IBackground* aBackgr, double aMinE=0, double aMaxE=1e300, double aMinZ=0, double aMaxZ=1e300);
	virtual double n(double aE, double aZ) const;
	virtual double MinE() const;
	virtual double MaxE() const;
	virtual double MinZ() const;
	virtual double MaxZ() const;
	virtual std::string Name() const;
private:
	SafePtr<IBackground> fBackgr;
	double fMinE;
	double fMaxE;
	double fMinZ;
	double fMaxZ;
};


typedef SmartPtr<BackgroundIntegral> PBackgroundIntegral;



} /* namespace mcray */
#endif /* BACKGROUND_H_ */
