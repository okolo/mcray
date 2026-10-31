/*
 * MathUtils.cpp
 *
 *  Created on: Oct 21, 2011
 *      Author: ok
 */

#include "MathUtils.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "TableFunction.h"

namespace Utils {

double MathUtils::SolveEquation(
			gsl_function F,	double x_lo, double x_hi,
			const double relError, const int max_iter,
			const gsl_root_fsolver_type *T)
{
    double r = 0;
	int status;
	int iter = 0;

	gsl_root_fsolver *solver = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (solver, &F, x_lo, x_hi);
	do
	{
	 	iter++;
	    status = gsl_root_fsolver_iterate (solver);
	    r = gsl_root_fsolver_root (solver);
	    x_lo = gsl_root_fsolver_x_lower (solver);
	    x_hi = gsl_root_fsolver_x_upper (solver);
	    status = gsl_root_test_interval (x_lo, x_hi, 0, relError);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	ASSERT(status == GSL_SUCCESS);
	gsl_root_fsolver_free (solver);
	if (status != GSL_SUCCESS)
	{
		if(iter >= max_iter)
			Exception::Throw("Failed to solve equation: maximal number of iterations achieved");
	  	else
	  		Exception::Throw("Failed to solve equation: GSL status " + ToString(status));
	}
	return r;
}

double MathUtils::SolveEquation(
		double (*aEquation) (double, void*),
		double x_lo, double x_hi, void* aEquationPars,
		const double relError, const int max_iter,
		const gsl_root_fsolver_type *T)
{
	gsl_function F;
	F.function = aEquation;
	F.params = aEquationPars;
	return SolveEquation(F,	x_lo, x_hi, relError, max_iter, T);
}

MathUtils::MathUtils():
gslQAGintegrator(0)
{

}

MathUtils::~MathUtils()
{
	if(gslQAGintegrator)
		gsl_integration_workspace_free (gslQAGintegrator);
}

template<typename X> bool MathUtils::SampleLogscaleDistribution(const Function& aDistrib, double aRand, X& aOutputX, X& aOutputIntegral, int nStepsS, X xMin, X xMax, double aRelError)
{
	std::vector<X> sArray,ratesArray;
	double distrXmin = aDistrib.Xmin();
	if(xMin < distrXmin)
		xMin = distrXmin;
	double distrXmax = aDistrib.Xmax();
	if(xMax > distrXmax)
		xMax = distrXmax;
	ASSERT(xMin>0. && xMin<xMax);
	double stepS = pow(xMax/xMin,1./nStepsS);

	sArray.push_back(0);
	ratesArray.push_back(0);

	X taleAccLimit;
	RelAccuracy<X>(taleAccLimit);
	taleAccLimit*=10;
	int maxIntervals = (int)(0.1/aRelError + 10.5);
	aOutputIntegral=0.;
	double deltaLogS = log(stepS);
	double s=xMin;
	for(int iS=1; iS<=nStepsS; iS++)
	{
		double s2=s*stepS;
		X rate = Integration_qag(aDistrib,s,s2,1e-300,aRelError,maxIntervals);
		ASSERT_VALID_NO(rate);
		if(rate==0. && aOutputIntegral==0.)
		{//move Smax
			sArray[0] = deltaLogS*iS;
		}
		else
		{
			if(rate==0 || (aOutputIntegral>0 && rate/aOutputIntegral < taleAccLimit))
				rate = aOutputIntegral*taleAccLimit;//avoiding zero derivative (inverse function must be defined everywhere)
			aOutputIntegral += rate;
			sArray.push_back(deltaLogS*iS);
			ratesArray.push_back(aOutputIntegral);
		}
		s=s2;
	}
	if(aOutputIntegral>0)
	{//sampling s
		ASSERT(sArray.size()>1);
		double randRate = aRand*aOutputIntegral;
		for(int i = ratesArray.size()-1;i>=0; i--)
			ratesArray[i]-=randRate;
		GSLTableFunc func(sArray, ratesArray, 0, 0, gsl_interp_cspline);
		func.SetAutoLimits();
		double logError = aRelError/(nStepsS*deltaLogS);
		if(logError>aRelError)
			logError = aRelError;
		aOutputX = SolveEquation(func, 0, nStepsS*deltaLogS, logError);
		aOutputX = xMin*exp(aOutputX);
		ASSERT(aOutputX>=xMin && aOutputX<=xMax);
		return true;
	}
	return false;
}

IFunctionCallHandlerX<double>* MathUtils::fLogger = 0;

template<typename X> void MathUtils::RelAccuracy(X& aOutput)
{
	NOT_IMPLEMENTED
}

template<> void MathUtils::RelAccuracy<double>(double& aOutput)
{
	aOutput = 1e-15;
}

template<> void MathUtils::RelAccuracy<long double>(long double& aOutput)
{
	aOutput = 1e-18L;
}

template bool MathUtils::SampleLogscaleDistribution<double>(const Function& aDistrib, double aRand, double& aOutputX, double& aOutputIntegral, int nStepsS, double xMin, double xMax, double aRelError);
//template bool MathUtils::SampleLogscaleDistribution<long double>(const Function& aDistrib, double aRand, long double& aOutput, int nStepsS, long double xMin, long double xMax, double aRelError);

double MathUtils::Integration_qag (
		gsl_function aFunction,
		double aXmin,
		double aXmax,
		double epsabs,
		double epsrel,
		size_t limit,
		int key)
{
	if(gslQAGintegrator==0)
		gslQAGintegrator = gsl_integration_workspace_alloc (limit);
	else if(gslQAGintegrator->limit < limit)
	{
		gsl_integration_workspace_free (gslQAGintegrator);
		gslQAGintegrator = gsl_integration_workspace_alloc (limit);
	}
	double result, abserr;
	try{
		if(epsabs==0)
			epsabs = std::numeric_limits<double>::min();
		int failed = gsl_integration_qag (&aFunction, aXmin, aXmax, epsabs, epsrel, limit, key, gslQAGintegrator, &result, &abserr);
		if(failed)
		{
			ASSERT(0);
			Exception::Throw("Integration failed with code " + ToString(failed));
		}
	}catch(Exception* ex)
	{
#ifdef _DEBUG
		GslProxyFunction f(aFunction,aXmin,aXmax);
		std::cerr << "\n\n#Integration_qag debug output:" << std::endl;
		bool logscale = aXmin>0 && aXmax/aXmin > 100;
		f.Print(std::cerr, 100, logscale, aXmin, aXmax);
		std::cerr << "\n\n#end of Integration_qag debug output" << std::endl;
#endif
		throw ex;
	}
	return result;
}

} /* namespace Utils */
