/*
 * Background.cpp
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


#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <gsl/gsl_sf_erf.h>

#include "Background.h"
#include "Utils.h"
#include "TableFunction.h"
#include "PrecisionTests.h"


namespace mcray {

using namespace std;

double BackgroundUtils::CalcIntegralDensity(const IBackground& aBackground, double aRelError)
{
	ASSERT_VALID_NO(aBackground.MinE());
	ASSERT_VALID_NO(aBackground.MaxE());
	double x1 = log(aBackground.MinE());
	double x2 = log(aBackground.MaxE());
	MathUtils math;
	gsl_function f;
	f.function = DensityCalculator;
	f.params = (void*)(&aBackground);
	return math.Integration_qag(f, x1, x2, 0, aRelError, 1000);
}

void BackgroundUtils::Print(const IBackground* aBackground, double aZ, std::ostream& aOut, int nIntervals)
{
	class Spec : public Function
	{
		const IBackground* fBackground;
		double fZ;
		virtual double f(double aE_eV) const
		{
			return fBackground->n(aE_eV*units.eV, fZ)*units.cm3;
		}
		virtual double Xmin() const {return fBackground->MinE()/units.eV; }
		virtual double Xmax() const {return fBackground->MaxE()/units.eV; }
	public:
		Spec(const IBackground* aBackground, double aZ):fBackground(aBackground), fZ(aZ){}
	};
	Spec spec(aBackground, aZ);
	spec.Print(aOut, nIntervals, true);
}

double BackgroundUtils::DensityCalculator (double x, void * params)
{
	IBackground* pBackgr = (IBackground*)params;
	return pBackgr->n(exp(x),0);
}

void BackgroundUtils::UnitTest(const IBackground& aBackground, int aPointsZ, int aPointsE)
{
	double zMax=aBackground.MaxZ();
	if(zMax==0)
		aPointsZ = 1;
	for(int iZ=0; iZ<aPointsZ; iZ++){
		double z = zMax/aPointsZ*iZ;
		ofstream out((aBackground.Name() + "z" + ToString(z)).c_str());
		Print(&aBackground, z, out, aPointsE);
	}
	std::cerr << "density " << CalcIntegralDensity(aBackground)*units.cm3 << std::endl;
}

TableBackground::TableBackground(std::string aName):
		fZtoIndex(0),
		fMinE(DBL_MAX),
		fMaxE(0),
		fName(aName)
{

}

void TableBackground::Init(std::string aDir, const char** aFileList, bool aIsLogscaleY)
{
	//the files should be ordered increasingly?
	//x scale should by monotonic with increasing order
	aDir = TABLES_DIR + aDir;
	double lastZ = -1;
	for(int iZ=0; aFileList[iZ]; iZ++)
	{
		double z = -1.;
		z = atof(aFileList[iZ]);
		if(z<0)
			throw "Invalid table background file name";
		if(z<lastZ)
			throw "Invalid table background file order";

		std::ifstream file;
		file.open((aDir + DIR_DELIMITER_STR + aFileList[iZ]).c_str(), std::ios::in);
		if(!file.is_open())
			Exception::Throw("Failed to open " + aDir + DIR_DELIMITER_STR + aFileList[iZ]);
		TableReader* reader = new TableReader(file, 2);
		Vector& x = reader->getColumn(0);
		Vector& y = reader->getColumn(1);
		int iMax = x.size();
		for(int iX=0; iX<iMax; iX++)
		{//TODO add support of arbitrary ordered data
			double E = recordToE(x[iX],z);
			if(E<fMinE)
				fMinE=E;
			else if(E>fMaxE)
				fMaxE=E;
		}
		double leftValue = aIsLogscaleY?(y[0]-1000):0;
		double rightValue = aIsLogscaleY?(y[y.size()-1]-1000):0;

		GSLTableFunc* f = new GSLTableFunc(reader,gsl_interp_linear,leftValue,rightValue);
		fData.push_back(f);
		fZ.push_back(z);
		fZindex.push_back((double)iZ);
	}
	fZtoIndex = new GSLTableFunc(fZ, fZindex, gsl_interp_linear);
}

TableBackground::~TableBackground() {
	delete fZtoIndex;
	for(int i = fData.size()-1; i>=0; i--)
		delete fData[i];
}

MatrixBackground::MatrixBackground(std::string aName, std::string aTableFile, bool aIsComoving, bool aExtendToZero):
		fZtoIndex(0),
		fLogEtoIndex(0),
		fName(aName),
		fExtendToZero(aExtendToZero),
		fIsComoving(aIsComoving),
		fBuffer(0)
{
	const size_t bufSize = 1048576;//1M
	fBuffer = new char[bufSize];
	std::istream_iterator<double> eos;
	aTableFile = TABLES_DIR + aTableFile;
	std::ifstream file(aTableFile.c_str());

	if(!file.good())
		Exception::Throw("Invalid background table file " + aTableFile);
	file.getline(fBuffer, bufSize);
	std::istringstream ist(fBuffer);
	std::istream_iterator<double> iit (ist);
	for(; iit!=eos; iit++)
	{
		fZ.push_back(*iit);
		fZindex.push_back(fZ.size()-1);
	}
	double Epriv = 0;
	while(!file.eof())
	{
		file.getline(fBuffer, bufSize);
		std::istringstream ist(fBuffer);
		std::istream_iterator<double> iit (ist);
		if(iit==eos)
			break;
		double E = *iit;
		if(E<=0 || E<Epriv)
			Exception::Throw("Invalid background table file " + aTableFile + ": invalid energy value in line " +
				ToString(fLogEdN_dE.size()+1));
		Epriv = E;
		fLogE.push_back(log(E*1e-6/units.Eunit));
		fLogEindex.push_back(fLogE.size()-1);
		std::vector<double>* pN = new std::vector<double>();
		fLogEdN_dE.push_back(pN);
		for(iit++; iit!=eos; iit++)
		{
			double dN_dE = *iit;
			if(dN_dE<0)
				Exception::Throw("Invalid background table file " + aTableFile + ": line " +
									ToString(fLogEdN_dE.size()+1) + " contains " + ToString(pN->size()+1) +
									" negative value");
			pN->push_back(dN_dE>0?log(E*dN_dE*units.Vunit):-700.);
		}
		if(pN->size()!=fZ.size())
			Exception::Throw("Invalid background table file " + aTableFile + ": line " +
					ToString(fLogEdN_dE.size()+1) + " contains " + ToString(pN->size()+1) +
					" records while " + ToString(fZ.size()+1) + " expected");
	}
	fZtoIndex = new GSLTableFunc(fZ, fZindex, gsl_interp_linear, aExtendToZero ? 0. : -1., -1.);
	fLogEtoIndex = new GSLTableFunc(fLogE, fLogEindex, gsl_interp_linear, -1., -1.);
	fMinE = exp(fLogE[0]);
	fMaxE = exp(*(fLogE.end()-1));
	delete fBuffer;
	fBuffer = 0;
}

MatrixBackground::~MatrixBackground()
{
	delete fBuffer;
	delete fZtoIndex;
	delete fLogEtoIndex;
	for(std::vector<Vector*>::iterator  it = fLogEdN_dE.begin(); it != fLogEdN_dE.end(); it++)
		delete (*it);
}

double MatrixBackground::n(double aE, double aZ) const
{
	double logE = log(aE);
	double indexE = fLogEtoIndex->f(logE);
	if(indexE<0.)
		return 0.;
	double indexZ = fZtoIndex->f(aZ);
	if(indexZ<0.)
		return 0.;
	int iE = indexE;
	double fracE = indexE-iE;
	int iZ = indexZ;
	double fracZ = indexZ-iZ;
	double result = fLogEdN_dE[iE]->at(iZ)*(1.-fracZ)*(1.-fracE);
	if(fracZ>0)
		result += (fLogEdN_dE[iE]->at(iZ+1)*fracZ*(1.-fracE));
	if(fracE>0)
	{
		result += (fLogEdN_dE[iE+1]->at(iZ)*(1.-fracZ)*fracE);
		if(fracZ>0)
			result += (fLogEdN_dE[iE+1]->at(iZ+1)*fracZ*fracE);
	}
	result = exp(result);
	if(fIsComoving)
	{
		double z1=1.+aZ;
		result *= (z1*z1*z1);
	}
	return result;
}

PlankBackground::PlankBackground(double aT0, double aEmin, double aEmax, double aMinZ, double aMaxZ, bool aFixedTemperature):
	fT0(aT0),
	fMinE(aEmin),
	fMaxE(aEmax),
	fMinZ(aMinZ),
	fMaxZ(aMaxZ),
	fFixedTemperature(aFixedTemperature)
{
	ASSERT(aEmin<aEmax);
	ASSERT(aMinZ<aMaxZ);
	fName = "PlankT=" + ToString(aT0*units.Eunit*Units::phTemperature_mult) + "K";
}

double PlankBackground::n(double aE, double aZ) const
{
	double gamma=aE/fT0;
	if(!fFixedTemperature)
		gamma /= (1.0+aZ);
	if(gamma<700)
	{
		double ex_1 = (gamma > 1e-2) ? exp(gamma)-1.0 : gamma*(1.+ 0.5*gamma*(1. + 0.333333333333333*gamma*(1. + 0.25*gamma)));
		return aE*aE*aE/9.869604401089/ex_1;//that is E^3/Pi^2/(...)
	}
	return 0;
}

	/*!
	@param[in]  aE background photons energy central value (internal units)
	@param[in]  aSigma Gaussian distribution sigma (internal units)
	@param[in]  aWidth width at which the distribution is cut (in units of sigma)
	@param[in]  aConcentration integral concentration of photons (internal units)
	@param[in]  aMinZ minimal background red shift
	@param[in]  aMaxZ maximal background red shift
	*/
GaussianBackground::GaussianBackground(double aE, double aSigma, double aWidth,
		double aConcentration, double aMinZ, double aMaxZ):
		fE0(aE), fSigma(aSigma), fMinZ(aMinZ), fMaxZ(aMaxZ), fWidth(aWidth)
{
	fMinE = fE0 - 0.5*fSigma*aWidth;
	fMaxE = fE0 + 0.5*fSigma*aWidth;
	fNorm = 1./fSigma/sqrt(2.*M_PI);
	double mult = 1./sqrt(2.)/fSigma;
	double conc = 0.5*(gsl_sf_erf(mult*(fMaxE-fE0))-gsl_sf_erf(mult*(fMinE-fE0)));
	fNorm *= aConcentration/conc;
	fName = "Gaussian ( E0 = " + ToString(fE0*units.Eunit*1e6) + " eV; sigma/E0 = " +
			ToString(fSigma/fE0)  + "; width/sigma = " + ToString(aWidth)  +
			"; n = " + ToString(aConcentration/units.Vunit) + "cm^-3 )";
}

double GaussianBackground::n(double aE, double aZ) const
{
	double ratio = fabs((aE-fE0)/fSigma);
	if(ratio>fWidth)
		return 0.;
	return fNorm*exp(-0.5*ratio*ratio)*aE;
}

double TableBackground::n(double aE, double aZ) const
{
	ASSERT(fMaxE>fMinE);//make sure Init() method was called prior
	int iZ = (int)fZtoIndex->f(aZ);
	if(iZ<0)
		return 0;
	double x = EtoRecordX(aE,aZ);
	double leftY = fData[iZ]->f(x);
	double y=0.;
	if(iZ==((int)fZ.size())-1)
		y = leftY;
	else
	{
		double rightY = fData[iZ+1]->f(x);
		double zL = fZ[iZ];
		double zR = fZ[iZ+1];
		y = leftY + (rightY-leftY)/(zR-zL)*(aZ-zL);
	}
	return recordToN(x, y, aZ);
}

double TableBackground::MinE() const
{
	ASSERT(fMaxE>fMinE);//make sure Init() method was called prior
	return fMinE;
}

double TableBackground::MaxE() const
{
	ASSERT(fMaxE>fMinE);//make sure Init() method was called prior
	return fMaxE;
}

double TableBackground::MinZ() const
{
	ASSERT(fMaxE>fMinE);//make sure Init() method was called prior
	return fZ[0];
}

double TableBackground::MaxZ() const
{
	ASSERT(fMaxE>fMinE);//make sure Init() method was called prior
	return fZ[fZ.size()-1];
}

void CompoundBackground::AddComponent(IBackground* aBackground, double aWeight)
{
	fBackgrounds.push_back(aBackground);
	fWeights.push_back(aWeight);
}

double CompoundBackground::n(double aE, double aZ) const
{
	double sum = 0.;
	vector<double>::const_iterator wit = fWeights.begin();
	for(vector<IBackground*>::const_iterator it = fBackgrounds.begin(); it != fBackgrounds.end(); it++, wit++)
	{
		const IBackground* b = *it;
		double weight = *wit;
		if(aE>=b->MinE() && aE<=b->MaxE() && aZ>=b->MinZ() && aZ<=b->MaxZ())
			sum += weight*b->n(aE,aZ);
	}
	return sum;
}

std::string CompoundBackground::Name() const
{
	std::string result = "";
	for(vector<IBackground*>::const_iterator it = fBackgrounds.begin(); it != fBackgrounds.end(); it++)
	{
		const IBackground* b = *it;
		if(result.length()>0)
			result = result + " + " + b->Name();
		else
			result = b->Name();
	}
	return result;
}

double CompoundBackground::MinE() const
{
	double result = DBL_MAX;
	for(vector<IBackground*>::const_iterator it = fBackgrounds.begin(); it != fBackgrounds.end(); it++)
	{
		const IBackground* b = *it;
		double min = b->MinE();
		if(result>min)
			result = min;
	}
	return result;
}

double CompoundBackground::MaxE() const
{
	double result = 0;
	for(vector<IBackground*>::const_iterator it = fBackgrounds.begin(); it != fBackgrounds.end(); it++)
	{
		const IBackground* b = *it;
		double max = b->MaxE();
		if(result<max)
			result = max;
	}
	return result;
}

double CompoundBackground::MinZ() const
{
	double result = DBL_MAX;
	for(vector<IBackground*>::const_iterator it = fBackgrounds.begin(); it != fBackgrounds.end(); it++)
	{
		const IBackground* b = *it;
		double min = b->MinZ();
		if(result>min)
			result = min;
	}
	return result;
}

double CompoundBackground::MaxZ() const
{
	double result = 0;
	for(vector<IBackground*>::const_iterator it = fBackgrounds.begin(); it != fBackgrounds.end(); it++)
	{
		const IBackground* b = *it;
		double max = b->MaxZ();
		if(result<max)
			result = max;
	}
	return result;
}

CompoundBackground::~CompoundBackground()
{
	for(vector<IBackground*>::iterator it = fBackgrounds.begin(); it != fBackgrounds.end(); it++)
		delete *it;
}

BackgroundIntegral* ContinuousBackgroundIntegral::Clone() const
{
	return new ContinuousBackgroundIntegral(*this);
}

ContinuousBackgroundIntegral::ContinuousBackgroundIntegral(const ContinuousBackgroundIntegral& aBackgroundIntegral):
		fStepZ(aBackgroundIntegral.fStepZ),
		fLogStepK(aBackgroundIntegral.fLogStepK),
		fZmaxLimit(aBackgroundIntegral.fZmaxLimit),
		fMaxK(aBackgroundIntegral.fMaxK),
		fMinK(aBackgroundIntegral.fMinK),
		fEpsRel(aBackgroundIntegral.fEpsRel),
		fBinningK(aBackgroundIntegral.fBinningK),
		fEnableDebugOutput(aBackgroundIntegral.fEnableDebugOutput)
{
	int size = aBackgroundIntegral.fIntegrals.size();
	for(int iZ=0; iZ<size; iZ++)
	{
		fIntegrals.push_back(aBackgroundIntegral.fIntegrals[iZ]->Clone());
	}
}

ContinuousBackgroundIntegral::ContinuousBackgroundIntegral(IBackground& aBackground, double aStepZ, double aLogStepK, double aZmaxLimit, double aEpsRel):
fStepZ(aStepZ),
fLogStepK(aLogStepK),
fZmaxLimit(aZmaxLimit),
fEpsRel(aEpsRel),
fEnableDebugOutput(false)
{
	ASSERT(fEpsRel<=0.1);
	if(aLogStepK<=1)
		Exception::Throw("invalid BackgroundIntegral::BackgroundIntegral argument aLogStepE");
	fMaxK = aBackground.MaxE();
	fMinK = aBackground.MinE();
	if(aBackground.MaxZ()<fZmaxLimit)
		fZmaxLimit=aBackground.MaxZ();
	int nIntervalsZ = (int)(fZmaxLimit/fStepZ + 1.0);
	ASSERT(nIntervalsZ>=0);
	fStepZ = (nIntervalsZ>0)?(fZmaxLimit/((double)nIntervalsZ)):0;
	int nIntervalsE = (int)(log10(fMaxK/fMinK)/log10(fLogStepK) + 0.5);
	if(nIntervalsE<10)
		nIntervalsE = 2;
	fLogStepK = pow(fMaxK/fMinK, 1./nIntervalsE);
	if(fLogStepK<=1)
			Exception::Throw("invalid BackgroundIntegral::BackgroundIntegral argument aBackground (zero energy range)");

	{double k=fMinK;
	for(int iK=0; iK<=nIntervalsE; iK++, k*=fLogStepK)
	{//build energy binning
		fBinningK.push_back(k);
	}}
	Kern kern(aBackground);
	{for(int iZ=0; iZ<=nIntervalsZ; iZ++)
	{
		fData.push_back(new Vector(nIntervalsE+1,0.));
		double z=iZ*fStepZ;
		kern.Z = z;
		double k=fMaxK/fLogStepK;
		double sum = 0.;
		(*fData[iZ])[nIntervalsE]=0.;//I(k_max)=0
		for(int iK=nIntervalsE-1; iK>=0; iK--,k/=fLogStepK)
		{
			double I = math.Integration_qag(kern, k, k*fLogStepK, 0, fEpsRel, (int)(1./fEpsRel));
			ASSERT_VALID_NO(I);
			sum += I;
			(*fData[iZ])[iK]=sum;
		}
	}}
	{for(int iZ=0; iZ<=nIntervalsZ; iZ++)
	{
		GSLTableFunc* func = new GSLTableFunc(fBinningK, *fData[iZ], new LogScale(), new LogScale(), gsl_interp_linear);
		func->SetAutoLimits();//for k<k_min I(k)=I(k_min); for k>k_max I(k)=I(k_max)=0
		fIntegrals.push_back(func);
	}}
}

double ContinuousBackgroundIntegral::Kern::f(double aX) const
{
	// TODO think about switching to log scale integration
	return fBackground.n(aX, Z)/(aX*aX*aX);
}

ContinuousBackgroundIntegral::RateSKern::RateSKern(const Function&	aBackrgoundIntegral, const Function& aSigma, const Particle& aParticle) :
	fBackrgoundIntegral(aBackrgoundIntegral),fSigma(aSigma)
{
	fM2 = aParticle.Mass();
	fM2*=fM2;
	fMaxEmult = 2.0*aParticle.Energy*(1.-aParticle.beta());
	if(fMaxEmult>0)
		fMaxEmult = 1./fMaxEmult;
	else
		fMaxEmult = -1.;

	fMult = 0.5/aParticle.Energy/(1.+aParticle.beta());
}
//(s-m^2) sigma(s) I((s-m^2)/2/E/(1+beta),z)
double ContinuousBackgroundIntegral::RateSKern::f(double s) const
{
	double bi = fBackrgoundIntegral((s-fM2)*fMult);
	if(fMaxEmult>0)
	{
		bi -= fBackrgoundIntegral((s-fM2)*fMaxEmult);
	}
	return fSigma(s)*(s-fM2)*bi;
}

ContinuousBackgroundIntegral::RateKKern::RateKKern(const Function&	aBackrgoundIntegral, const Function& aSigma, const Particle& aParticle) :
		fBackrgoundIntegral(aBackrgoundIntegral),fSigma(aSigma)
{
	double m = aParticle.Mass();

	fMaxEmult = aParticle.Energy*(1.-aParticle.beta());
	if(fMaxEmult>0)
		fMaxEmult = m/fMaxEmult;
	else
		fMaxEmult = -1.;

	fMult = m/aParticle.Energy/(1.+aParticle.beta());
}

double ContinuousBackgroundIntegral::RateKKern::f(double k) const
{
	double bi = fBackrgoundIntegral(k*fMult);
	if(fMaxEmult>0)
	{
		bi -= fBackrgoundIntegral(k*fMaxEmult);
	}
	return fSigma(k)*k*bi;
}

ContinuousBackgroundIntegral::RateSKernZ::RateSKernZ(const Function&	aBackrgoundIntegralZ1, const Function&	aBackrgoundIntegralZ2, double aZ1, double aZ2, const Function& aSigma, const Particle& aParticle):
		RateSKern(aBackrgoundIntegralZ1, aSigma, aParticle),
		fBackrgoundIntegralZ1(aBackrgoundIntegralZ1),
		fBackrgoundIntegralZ2(aBackrgoundIntegralZ2)
{
	//initializing linear interpolation coefficients
	fCoef2 = (aParticle.Time.z()-aZ1)/(aZ2-aZ1);
	fCoef1 = 1.-fCoef2;
}

//(s-m^2) sigma(s) I((s-m^2)/2/E/(1+beta),z)
double ContinuousBackgroundIntegral::RateSKernZ::f(double s) const
{
	double k = (s-fM2)*fMult;
	double bi1 = fBackrgoundIntegralZ1(k);
	double bi2 = fBackrgoundIntegralZ2(k);
	if(fMaxEmult>0)
	{
		k = (s-fM2)*fMaxEmult;
		bi1 -= fBackrgoundIntegralZ1(k);
		bi2 -= fBackrgoundIntegralZ2(k);
	}
	double bi = fCoef1*bi1 + fCoef2*bi2;
	double result = fSigma(s)*(s-fM2)*bi;
	//std::cout.precision(15);
	//std::cout << s << "\t" << result << std::endl;
	return result;
}

ContinuousBackgroundIntegral::RateKKernZ::RateKKernZ(const Function&	aBackrgoundIntegralZ1, const Function&	aBackrgoundIntegralZ2, double aZ1, double aZ2, const Function& aSigma, const Particle& aParticle):
		RateKKern(aBackrgoundIntegralZ1, aSigma, aParticle),
		fBackrgoundIntegralZ1(aBackrgoundIntegralZ1),
		fBackrgoundIntegralZ2(aBackrgoundIntegralZ2)
{
	//initializing linear interpolation coefficients
	fCoef2 = (aParticle.Time.z()-aZ1)/(aZ2-aZ1);
	fCoef1 = 1.-fCoef2;
}

double ContinuousBackgroundIntegral::RateKKernZ::f(double k) const
{
	double eps = k*fMult;
	double bi1 = fBackrgoundIntegralZ1(eps);
	double bi2 = fBackrgoundIntegralZ2(eps);
	if(fMaxEmult>0)
	{
		eps = k*fMaxEmult;
		bi1 -= fBackrgoundIntegralZ1(eps);
		bi2 -= fBackrgoundIntegralZ2(eps);
	}
	double bi = fCoef1*bi1 + fCoef2*bi2;
	double result = fSigma(k)*k*bi;
	return result;
}

Function* ContinuousBackgroundIntegral::GetRateDistributionS(const Function& aSigma, const Particle& aParticle, double& aSmin, double& aSmax)
{
    double SminSigma = aSigma.Xmin();
    ASSERT_VALID_NO(SminSigma);
    double beta = aParticle.beta();
    double m=aParticle.Mass();
    aSmax = m*m+2.*aParticle.Energy*fMaxK*(1.+beta);
    double SminKinematic = m*m+2.*aParticle.Energy*fMinK*(1.-beta);
    aSmin = SminKinematic<SminSigma ? SminSigma : SminKinematic;
    if(aSigma.Xmax()<aSmax)
        aSmax=aSigma.Xmax();
    if(aSmax<=aSmin)
        return NULL;

    double z = aParticle.Time.z();
    int iZ= fStepZ>0 ? ((int)(z/fStepZ)):0;
    //int maxIntervals = (int)(0.1/fEpsRel + 10.5);

    //int nStepsS = (int)(log10(Smax/Smin)/log10(fLogStepK)+0.5);
    //if(nStepsS<5)
    //	nStepsS = 5;
    //double stepS = pow(Smax/Smin,1./nStepsS);
    double totalRate = 0.;//rate times 8*E*beta

    bool zDependenceExists = (iZ<(int)fIntegrals.size()-1);

    if(zDependenceExists)
        return new RateSKernZ(*fIntegrals[iZ], *fIntegrals[iZ + 1], iZ * fStepZ, (iZ + 1) * fStepZ, aSigma, aParticle);
    else
        return new RateSKern(*fIntegrals[iZ], aSigma, aParticle);
}

double ContinuousBackgroundIntegral::GetRateAndSampleS(const Function& aSigma, const Particle& aParticle, double aRand, double& aS, double aAbsError)
{
    double Smin, Smax, totalRate;
    SafePtr<Utils::Function> kern = GetRateDistributionS(aSigma, aParticle, Smin, Smax);
	if(MathUtils::SampleLogDistribution(*kern, aRand, aS, totalRate, Smin, Smax, fEpsRel))
	{
		return totalRate/8./aParticle.Energy/aParticle.Energy/aParticle.beta();
	}
	aS=0;
	return 0;
}

double ContinuousBackgroundIntegral::GetRateAndSampleK(const Function& aSigmaK, const Particle& aParticle, double aRand, double&aK, double aAbsError)
{
	ASSERT(aRand>0 || aRand<1);
	double KminSigma = aSigmaK.Xmin();
	ASSERT_VALID_NO(KminSigma);
	double beta = aParticle.beta();
	double m=aParticle.Mass();
	double Kmax = aParticle.Energy * fMaxK * (1. + beta)/m;
	double KminKinematic = aParticle.Energy * fMinK * (1. - beta)/m;
	double Kmin = KminKinematic < KminSigma ? KminSigma : KminKinematic;
	if(aSigmaK.Xmax() < Kmax)
		Kmax = aSigmaK.Xmax();
	if(Kmax <= Kmin)
		return 0.;

	double z = aParticle.Time.z();
	int iZ= fStepZ>0 ? ((int)(z/fStepZ)):0;

	double totalRate = 0.;//rate times 8*E*beta

	bool zDependenceExists = (iZ<(int)fIntegrals.size()-1);
	SafePtr<RateKKern> kern;
	if(zDependenceExists)
		kern = new RateKKernZ(*fIntegrals[iZ], *fIntegrals[iZ + 1], iZ * fStepZ, (iZ + 1) * fStepZ,
											 aSigmaK, aParticle);
	else
		kern = new RateKKern(*fIntegrals[iZ], aSigmaK, aParticle);

	if(MathUtils::SampleLogDistribution(*kern, aRand, aK, totalRate, Kmin, Kmax, fEpsRel))
	{
		double gamma=aParticle.Energy/m;
		return totalRate/(2.0*gamma*gamma*aParticle.beta());
	}
	aK =0;
	return 0;
}

double ContinuousBackgroundIntegral::GetRateS(const Function &aSigma, const Particle &aParticle, double aAbsError)
{
	double SminSigma = aSigma.Xmin();
	ASSERT_VALID_NO(SminSigma);
	double beta = aParticle.beta();
	double m=aParticle.Mass();
	double Smax = m*m+2.*aParticle.Energy*fMaxK*(1.+beta);
	double SminKinematic = m*m+2.*aParticle.Energy*fMinK*(1.-beta);
	double Smin = SminKinematic<SminSigma ? SminSigma : SminKinematic;

	if(Smax<=Smin)
		return 0.;
	if(aSigma.Xmax()<Smax)
		Smax=aSigma.Xmax();

	double z = aParticle.Time.z();
	int iZ= fStepZ>0 ? ((int)(z/fStepZ)):0;
	int maxIntervals = (int)(log10(Smax/Smin)/log10(fLogStepK)/fEpsRel*1e-2 + 100);

	RateSKern leftKern(*fIntegrals[iZ], aSigma, aParticle);

	double leftRate = math.Integration_qag(leftKern,Smin,Smax,0,fEpsRel,maxIntervals);
	ASSERT_VALID_NO(leftRate);

	double result = 0;
	if(iZ==(int)fIntegrals.size()-1)
		result = leftRate;
	else
	{
		RateSKern rightKern(*fIntegrals[iZ + 1], aSigma, aParticle);
		double rightRate = math.Integration_qag(rightKern,Smin,Smax,aAbsError,fEpsRel,maxIntervals);
		ASSERT_VALID_NO(rightRate);
		result = leftRate + (rightRate-leftRate)/fStepZ*(z-iZ*fStepZ);
		ASSERT_VALID_NO(result);
	}
	return result/8./aParticle.Energy/aParticle.Energy/aParticle.beta();
}

double ContinuousBackgroundIntegral::GetRateK(const Function &aSigmaK, const Particle &aParticle, double aAbsError)
{
	double KminSigma = aSigmaK.Xmin();
	ASSERT_VALID_NO(KminSigma);
	double beta = aParticle.beta();
	double m=aParticle.Mass();
	double Kmax = aParticle.Energy * fMaxK * (1. + beta)/m;
	double KminKinematic = aParticle.Energy * fMinK * (1. - beta)/m;
	double Kmin = KminKinematic < KminSigma ? KminSigma : KminKinematic;

	if(Kmax <= Kmin)
		return 0.;
	if(aSigmaK.Xmax() < Kmax)
		Kmax = aSigmaK.Xmax();

	double z = aParticle.Time.z();
	int iZ= fStepZ>0 ? ((int)(z/fStepZ)):0;
	int maxIntervals = (int)(log10(Kmax / Kmin) / log10(fLogStepK) / fEpsRel * 1e-2 + 100);

	RateKKern leftKern(*fIntegrals[iZ], aSigmaK, aParticle);

	double leftRate = math.Integration_qag(leftKern, Kmin, Kmax, 0, fEpsRel, maxIntervals);
	ASSERT_VALID_NO(leftRate);

	double result = 0;
	if(iZ==(int)fIntegrals.size()-1)
		result = leftRate;
	else
	{
		RateKKern rightKern(*fIntegrals[iZ + 1], aSigmaK, aParticle);
		double rightRate = math.Integration_qag(rightKern, Kmin, Kmax, aAbsError, fEpsRel, maxIntervals);
		ASSERT_VALID_NO(rightRate);
		result = leftRate + (rightRate-leftRate)/fStepZ*(z-iZ*fStepZ);
		ASSERT_VALID_NO(result);
	}
	double gamma=aParticle.Energy/m;
	return result/(2.0*gamma*gamma*aParticle.beta());
}

ContinuousBackgroundIntegral::~ContinuousBackgroundIntegral()
{
	//ASSERT(fData.size()==fIntegrals.size());//fData may be empty if this object was created with Clone() method
	for(int iZ=fData.size()-1; iZ>=0; iZ--)
	{
		delete fData[iZ];
	}
	for(int iZ=fIntegrals.size()-1; iZ>=0; iZ--)
	{
		delete fIntegrals[iZ];
	}
	fData.clear();
	fIntegrals.clear();
}

void ContinuousBackgroundIntegral::UnitTest()
{
	double zMax = 1;
	double kStep = pow(10,0.05);
	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;

	////// Background integral calculation test

	double kMin = 1e-3*cmbTemp;
	double kMax = 1e3*cmbTemp;
	PlankBackground plank(cmbTemp, kMin, kMax, 0, zMax);
	ContinuousBackgroundIntegral bi(plank, 0.1, kStep, zMax, 1e-3);
	Function& I = *(bi.fIntegrals[0]);
	double kMaxI = PlankBackgroundIntegral(cmbTemp, kMax);
	std::cout << "#Z=0" << "\n";
	for(double k=0.5*kMin; k<2*kMax; k*=kStep)
	{
		std::cout << k/cmbTemp << "\t" << I(k) << "\t" << PlankBackgroundIntegral(cmbTemp, k) - kMaxI << '\n';
	}
	std::cout << "\n\n#table z=0\n";
	I.Print(std::cout,1./cmbTemp,1.);

	std::cout << "\n\n#Z=" << zMax << "\n";
	Function& I1 = *(bi.fIntegrals[bi.fIntegrals.size()-1]);
	kMaxI = PlankBackgroundIntegral((1.+zMax)*cmbTemp, kMax);
	for(double k=0.5*kMin; k<2*kMax; k*=kStep)
	{
		std::cout << k/cmbTemp << "\t" << I1(k) << "\t" << PlankBackgroundIntegral((1.+zMax)*cmbTemp, k) - kMaxI << '\n';
	}

	///// Total rate and sampling test z=0

	class TestSigma : public Function
	{
	public:
		virtual double f(double s) const { return 1./s;}
		virtual double Xmin() const {return sMin;}
		virtual double Xmax() const {return sMax;}
		double sMin;
		double sMax;
	};

	kMin = 1e-6*cmbTemp;
	kMax = 1e-5*cmbTemp;
	PlankBackground plankLow(cmbTemp, kMin, kMax, 0, zMax);
	ContinuousBackgroundIntegral biLow(plankLow, 0.1, kStep, zMax, 1e-3);
	TestSigma ts;
	Particle p(Photon, 0.);
	p.Energy = 10;
	ts.sMin = 2*kMin*4*p.Energy;
	ts.sMax = 0.5*kMax*4*p.Energy;
	double rate = biLow.GetRateS(ts, p);
	kMaxI = PlankBackgroundIntegral(cmbTemp, kMax);
	double thRate = PlankBackgroundTestRate(cmbTemp,p.Energy,ts.sMax) -
			PlankBackgroundTestRate(cmbTemp,p.Energy,ts.sMin) - kMaxI*(ts.sMax-ts.sMin)/8./p.Energy/p.Energy;
	double S;
	biLow.fEnableDebugOutput = true;
	double rateSample = biLow.GetRateAndSampleS(ts, p, 0.5, S);
	biLow.fEnableDebugOutput = false;

	cout << "#Total rate test (z=0): rate/sample rate/thRate : " << rate << " / " << rateSample <<  " / " << thRate << "\n";
	cout << "#<rand> <S> <Smin> <Smax>\n";
	for(double rand=0.; rand<=1.; rand+=0.1)
	{
		biLow.GetRateAndSampleS(ts, p, rand, S);
		cout << rand << "\t" << S << "\t" << ts.sMin << "\t" << ts.sMax << "\n";
	}

	///// Total rate and sampling test z=0.55*zMax

	double z = 0.55*zMax;
	p.Time.setZ(z);
	cmbTemp *= (1.+z);
	rate = biLow.GetRateS(ts, p);

	kMaxI = PlankBackgroundIntegral(cmbTemp, kMax);
	thRate = PlankBackgroundTestRate(cmbTemp,p.Energy,ts.sMax) -
			PlankBackgroundTestRate(cmbTemp,p.Energy,ts.sMin) - kMaxI*(ts.sMax-ts.sMin)/8./p.Energy/p.Energy;

	biLow.fEnableDebugOutput = true;
	rateSample = biLow.GetRateAndSampleS(ts, p, 0.5, S);
	biLow.fEnableDebugOutput = false;

	cout << "#Total rate test(z=" << z << ") : rate/sample rate/thRate : " << rate << " / " << rateSample <<  " / " << thRate << "\n";
	cout << "#<rand> <S> <Smin> <Smax>\n";
	for(double rand=0.; rand<=1.; rand+=0.1)
	{
		biLow.GetRateAndSampleS(ts, p, rand, S);
		cout << rand << "\t" << S << "\t" << ts.sMin << "\t" << ts.sMax << "\n";
	}
}

double ContinuousBackgroundIntegral::PlankBackgroundTestRate(double aT, double aE, double aS)
{
	double x = 0.25*aS/aE/aT;
	return 0.5*aT*aT/aE/M_PI/M_PI*(0.25*x*x + x - x*log(x));
}

double ContinuousBackgroundIntegral::PlankBackgroundIntegral(double aT, double aK)
{
	double gamma = aK/aT;
	double result = 0;
	if(gamma<1e-5)
		result = 0.5*gamma-gamma*gamma/24.0+gamma*gamma*gamma*gamma/2880.-log(gamma);
	else if(gamma<700)
		result = gamma-log(exp(gamma)-1);
	else
		result = 0;
	return aT/M_PI/M_PI*result;
}

BackgroundIntegral* MonochromaticBackgroundIntegral::Clone() const
{
	return new MonochromaticBackgroundIntegral(*fConcentration.Clone(), *fEnergy.Clone(), fLogStepE, fEpsRel);
}

MonochromaticBackgroundIntegral::RateSKern::RateSKern(const Function& aSigma, const Particle& aParticle) :
	fSigma(aSigma)
{
	fM2 = aParticle.Mass();
	fM2*=fM2;
}
//(s-m^2) sigma(s) I((s-m^2)/2/E/(1+beta),z)
double MonochromaticBackgroundIntegral::RateSKern::f(double s) const
{
	return fSigma(s)*(s-fM2);
}

double MonochromaticBackgroundIntegral::RateKKern::f(double k) const
{
	return fSigma(k)*k;
}

Function* MonochromaticBackgroundIntegral::GetRateDistributionS(const Function& aSigma, const Particle& aParticle, double& aSmin, double& aSmax)
{
    double SminSigma = aSigma.Xmin();
    ASSERT_VALID_NO(SminSigma);
    double beta = aParticle.beta();
    double m=aParticle.Mass();

    double z = aParticle.Time.z();
    double n=fConcentration.f(z);
    if(n<=0.)
        return NULL;
    double k=fEnergy.f(z);

    aSmax = m*m+2.*aParticle.Energy*k*(1.+beta);
    double SminKinematic = m*m+2.*aParticle.Energy*k*(1.-beta);
    aSmin = SminKinematic<SminSigma ? SminSigma : SminKinematic;
    if(aSigma.Xmax()<aSmax)
        aSmax=aSigma.Xmax();
    if(aSmax<=aSmin)
        return NULL;

    return new RateSKern(aSigma, aParticle);
}

double MonochromaticBackgroundIntegral::GetRateAndSampleS(const Function& aSigma, const Particle& aParticle, double aRand, double& aS, double aAbsError)
{
	ASSERT(aRand>0 || aRand<1);
	aS=0;
    double z = aParticle.Time.z();
    double n=fConcentration.f(z);
    if(n<=0.)
        return NULL;
    double k=fEnergy.f(z);
    double Smin, Smax;
	double totalRate = 0.;//rate times 8*E*beta
	SafePtr<Function> kern = GetRateDistributionS(aSigma, aParticle, Smin, Smax);
    if (kern == NULL)
        return 0;

	if(MathUtils::SampleLogDistribution(*kern, aRand, aS, totalRate, Smin, Smax, fEpsRel))
	{
		return totalRate*n/k/k/8./aParticle.Energy/aParticle.Energy/aParticle.beta();
	}
	aS=0;
	return 0;
}

double MonochromaticBackgroundIntegral::GetRateS(const Function& aSigma, const Particle& aParticle, double aAbsError)
{
	double z = aParticle.Time.z();
	double n=fConcentration.f(z);
	if(n<=0.)
		return 0.;
	double k=fEnergy.f(z);
	double SminSigma = aSigma.Xmin();
	ASSERT_VALID_NO(SminSigma);
	double beta = aParticle.beta();
	double m=aParticle.Mass();
	double Smax = m*m+2.*aParticle.Energy*k*(1.+beta);
	double SminKinematic = m*m+2.*aParticle.Energy*k*(1.-beta);
	double Smin = SminKinematic<SminSigma ? SminSigma : SminKinematic;

	if(Smax<=Smin)
		return 0.;
	if(aSigma.Xmax()<Smax)
		Smax=aSigma.Xmax();

	int maxIntervals = (int)(log10(Smax/Smin)/log10(fLogStepE)/fEpsRel*1e-2 + 100);

	RateSKern kern(aSigma, aParticle);

	double rate = math.Integration_qag(kern,Smin,Smax,aAbsError,fEpsRel,maxIntervals);
	ASSERT_VALID_NO(rate);

	return rate*n/k/k/8./aParticle.Energy/aParticle.Energy/aParticle.beta();
}

HighRedshiftBackgrExtension::HighRedshiftBackgrExtension(IBackground* aBackground, double aDeltaZconst, double aPowerLow, double aDeltaZexp, double aZmax):
		fBackground(aBackground),
		fPowerLow(aPowerLow),
		fDeltaZexp(aDeltaZexp),
		fZmax(aZmax),
		fInnerZmax(fBackground->MaxZ())
{
	if(aDeltaZconst<0)
		Exception::Throw("HighRedshiftBackgrExtension: invalid DeltaZconst");
	fZconst=fInnerZmax+aDeltaZconst;
}

double HighRedshiftBackgrExtension::n(double aE, double aZ) const
{
	if(aZ<=fInnerZmax)
		return fBackground->n(aE, aZ);
	if(aZ>fZmax)
		return 0.;
	double expansionFactor = (1.+aZ)/(1+fInnerZmax);
	expansionFactor *= (expansionFactor*expansionFactor);// comoving volumes ratio

	if(aZ<fZconst)
		return expansionFactor*fBackground->n(aE, 0.9999*fInnerZmax);

	double powerLowFactor = pow((1.+aZ)/(1.+fZconst),fPowerLow);
	double expFactor = (fDeltaZexp>0)?exp((fZconst-aZ)/fDeltaZexp):1;

	return expansionFactor*powerLowFactor*expFactor*fBackground->n(aE, 0.9999*fInnerZmax);
}

double HighRedshiftBackgrExtension::MinE() const
{
	return fBackground->MinE();
}

double HighRedshiftBackgrExtension::MaxE() const
{
	return fBackground->MaxE();
}

double HighRedshiftBackgrExtension::MinZ() const
{
	return fBackground->MinZ();
}

double HighRedshiftBackgrExtension::MaxZ() const
{
	return fZmax;
}

std::string HighRedshiftBackgrExtension::Name() const
{
	std::string name = "extended ";
	name += fBackground->Name();
	return name;
}

HighRedshiftBackgrExtension::~HighRedshiftBackgrExtension()
{
	delete fBackground;
}

CuttedBackground::CuttedBackground(IBackground* aBackgr, double aMinE, double aMaxE, double aMinZ, double aMaxZ):
	fBackgr(aBackgr),fMinE(aMinE),fMaxE(aMaxE),fMinZ(aMinZ),fMaxZ(aMaxZ)
{
	if(aBackgr->MinZ()>aMinZ)
		fMinZ = aBackgr->MinZ();
	if(aBackgr->MinE()>aMinE)
		fMinE = aBackgr->MinE();
	if(aBackgr->MaxZ()<aMaxZ)
		fMaxZ = aBackgr->MaxZ();
	if(aBackgr->MaxE()<aMaxE)
		fMaxE = aBackgr->MaxE();
}

double CuttedBackground::n(double aE, double aZ) const {
	return (aE<=fMaxE && aE>=fMinE && aZ<=fMaxZ && aZ>=fMinZ) ? fBackgr->n(aE,aZ) : 0.;
}

double CuttedBackground::MinE() const {
	return fMinE;
}

double CuttedBackground::MaxE() const {
	return fMaxE;
}

double CuttedBackground::MinZ() const {
	return fMinZ;
}

double CuttedBackground::MaxZ() const {
	return fMaxZ;
}

std::string CuttedBackground::Name() const {
	return "Cutted " + fBackgr->Name();
}
double MonochromaticBackgroundIntegral::GetRateAndSampleK(const Function &aSigmaK, const Particle &aParticle,
														  Randomizer &aRandomizer, double &aK, double aAbsError) {
	double aRand = aRandomizer.Rand();
	aK=0;
	double z = aParticle.Time.z();
	double n=fConcentration.f(z);
	if(n<=0.)
		return 0.;
	double eps =fEnergy.f(z);

	double KminSigma = aSigmaK.Xmin();
	ASSERT_VALID_NO(KminSigma);
	double beta = aParticle.beta();
	double m=aParticle.Mass();

	double Kmax = aParticle.Energy * eps * (1. + beta)/m;
	double KminKinematic = aParticle.Energy * eps * (1. - beta)/m;
	double Kmin = KminKinematic < KminSigma ? KminSigma : KminKinematic;
	if(aSigmaK.Xmax() < Kmax)
		Kmax =aSigmaK.Xmax();
	if(Kmax <= Kmin)
		return 0.;

	double totalRate = 0.;
	RateKKern kern(aSigmaK);

	if(MathUtils::SampleLogDistribution(kern, aRand, aK, totalRate, Kmin, Kmax, fEpsRel))
	{
		double gamma = aParticle.Energy/m;
		return 0.5*totalRate * n /(eps * eps * gamma * gamma * aParticle.beta());
	}
	aK=0;
	return 0;
}

double MonochromaticBackgroundIntegral::GetRateK(const Function &aSigmaK, const Particle &aParticle,
												 double aAbsError) {
	double z = aParticle.Time.z();
	double n=fConcentration.f(z);
	if(n<=0.)
		return 0.;
	double eps = fEnergy.f(z);
	double KminSigma = aSigmaK.Xmin();
	ASSERT_VALID_NO(KminSigma);
	double KmaxSigma = aSigmaK.Xmin();
	double beta = aParticle.beta();
	double m=aParticle.Mass();
	double KmaxKinematic = aParticle.Energy * eps * (1. + beta)/m;
	double KminKinematic = aParticle.Energy * eps * (1. - beta)/m;

	double Kmin = KminKinematic<KminSigma ? KminSigma : KminKinematic;
	double Kmax = KmaxSigma<KmaxKinematic ? KmaxSigma : KmaxKinematic;

	if(Kmax<=Kmin)
		return 0.;

	int maxIntervals = (int)(log10(Kmax/Kmin)/log10(fLogStepE)/fEpsRel*1e-2 + 100);

	RateKKern kern(aSigmaK);

	double rate = math.Integration_qag(kern,Kmin,Kmax,aAbsError,fEpsRel,maxIntervals);
	ASSERT_VALID_NO(rate);
	double gamma = aParticle.Energy/m;

	return 0.5*rate * n /(eps * eps * gamma * gamma * aParticle.beta());
}
} /* namespace mcray */
