/*
 * MathUtils.h
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

#ifndef MATHUTILS_H_
#define MATHUTILS_H_

#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include "Utils.h"
#include <limits>
#include "Randomizer.h"

namespace Utils {


template<typename X = double >
class FunctionX
{
public:
	virtual X f(X _x) const = 0;
	virtual X Xmin() const {return -std::numeric_limits<X>::max();}
	virtual X Xmax() const {return std::numeric_limits<X>::max();}
	inline X operator()(X _x) const {return f(_x);};
	virtual ~FunctionX(){};

	virtual FunctionX<X>* Clone() const { NOT_IMPLEMENTED; return 0; }

	inline operator gsl_function() const
	{
		gsl_function result = {GslProxyFunc, (void*)this};
		return result;
	}
	void Print(std::ostream& aOut, int nIntervals, bool aLogScale, double aXmin=-DBL_MAX, double aXmax=DBL_MAX) const
	{
		if(aXmin<Xmin())
			aXmin = Xmin();
		if(aXmax>Xmax())
			aXmax = Xmax();
		ASSERT(aXmin<=aXmax);
		ASSERT(aXmax<DBL_MAX);
		ASSERT(aXmin>-DBL_MAX);
		ASSERT(nIntervals>=1);
		ASSERT((!aLogScale) || aXmin>0);
		double step = aLogScale?pow(aXmax/aXmin,1./nIntervals) : (aXmax-aXmin)/nIntervals;
		double x = aXmin;
		for(int i=0; i<=nIntervals; i++, (x = aLogScale ? x*step : x+step) )
			aOut << x << "\t" << f(x) << "\n";
	}
private:
	static double GslProxyFunc (double x, void * params)
	{
		const FunctionX<X>* f = (const FunctionX<X>*)params;
		return f->f(x);
	}
};

template<typename X = double >
class Function2X
{
public:
	virtual X f(X x, X y) const = 0;
	virtual X MinArg(int aArgNo) const {return std::numeric_limits<X>::min();}
	virtual X MaxArg(int aArgNo) const  {return std::numeric_limits<X>::max();}

	virtual inline X operator()(X _x, X _y) const {return f(_x,_y);};
	virtual ~Function2X(){};
};

template<typename X = double >
class IFunctionCallHandlerX
{
public:
	virtual void OnCall(const FunctionX<X>& aFunc, X aX, X aY) const = 0;
};

template<typename X = double >
class DebugFunctionX : public FunctionX<X>
{
public:
	DebugFunctionX(const FunctionX<X>& aOrigFunc, const IFunctionCallHandlerX<X>& aDebugger):fOrigFunc(aOrigFunc),fDebugger(aDebugger){};
	virtual X f(X _x) const
	{
		X y = fOrigFunc.f(_x);
		fDebugger.OnCall(fOrigFunc, _x, y);
		return y;
	}
	virtual ~DebugFunctionX(){};
	virtual FunctionX<X>* Clone() const { return new DebugFunctionX<X>(fOrigFunc,fDebugger); }
private:
	const FunctionX<X>&				fOrigFunc;
	const IFunctionCallHandlerX<X>&	fDebugger;
};

template<typename X = double >
class MaxFunctionX : public FunctionX<X>
{
public:
    MaxFunctionX(const FunctionX<X>& aOrigFunc, X& aMaxArg, X& aMaxValue):
    fOrigFunc(aOrigFunc),
    fMaxArg(aMaxArg),
    fMaxVal(aMaxValue)
    {};
    virtual X f(X _x) const
    {
        X y = fOrigFunc.f(_x);
        if (y > fMaxVal){
            fMaxArg = _x;
            fMaxVal = y;
        }
        return y;
    }
    virtual X Xmin() const {return fOrigFunc.Xmin();}
    virtual X Xmax() const {return fOrigFunc.Xmax();}
    virtual ~MaxFunctionX(){};
    virtual FunctionX<X>* Clone() const { return new MaxFunctionX<X>(fOrigFunc,fMaxArg,fMaxVal); }
private:
    const FunctionX<X>&				fOrigFunc;
    X&	                            fMaxArg;
    X&	                            fMaxVal;
};

template<typename X = double >
class FunctionCallLoggerX : public IFunctionCallHandlerX<X>
{
public:
	FunctionCallLoggerX(std::ostream& aOutput) : fOutput(aOutput){}
	virtual void OnCall(const FunctionX<X>& aFunc, X aX, X aY) const
	{
		((std::ostream&)fOutput) << aX << "\t" << aY << "\n";
	}
	virtual ~FunctionCallLoggerX()
	{
		fOutput << std::endl;
	}
private:
	std::ostream& fOutput;
};

typedef FunctionX<double> Function;
typedef Function2X<double> Function2;
typedef MaxFunctionX<double> MaxFunction;

template<typename X = double >
class ParamlessFunctionX : public FunctionX<X>
{
public:
	ParamlessFunctionX(X (*aFunction)(X)):fFunction(aFunction){}
	X f(X _x) const { return fFunction(_x); }
private:
	X (*fFunction)(X);
};

typedef ParamlessFunctionX<double> ParamlessFunction;

	template<typename X = double >
	class ISampler {
	public:
		virtual X sample(const FunctionX<X> &f, X aTmin, X aTmax, X aRand,
				 X &aTotRate, X aRelErr, X aAbsErr = 1e300) = 0;
		void UnitTest();
	};

/// Mathematical functions
/// Static methods are thread-safe
/// Non-static methods are not guaranteed to be thread-safe
	///TODO: make implementation using boost odeint and dense output and get rid of nr library
class MathUtils {
public:
	MathUtils();
	~MathUtils();
    static inline double RelDifference(double aVal1, double aVal2){
        double meanNorm = 0.5*fabs(aVal1)+fabs(aVal2);
        return meanNorm==0. ? 0.:fabs(aVal1-aVal2)/meanNorm;
    }

	static double SolveEquation(
			double (*aEquation) (double, void*),
			double x_lo, double x_hi, void* aEquationPars = 0,
			const double relError=1e-3, const int max_iter=100,
			const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection);

	static double SolveEquation(
			gsl_function f,	double x_lo, double x_hi,
			const double relError=1e-3, const int max_iter=100,
			const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection);

	double Integration_qag (
			gsl_function aFunction,
			double aXmin,
			double aXmax,
			double epsabs,
			double epsrel,
			size_t limit,
			int key=GSL_INTEG_GAUSS15);

	template<typename X> bool SampleLogscaleDistribution(const Function& aDistrib, double aRand, X& aOutputX, X& aOutputIntegral, int nStepsS, X xMin, X xMax, double aRelError);
	static bool _SampleDistribution(const Function& aDistrib, double aRand, double& aOutputX, double& aOutputIntegral, double xMin, double xMax, double aRelError);
    static bool _SampleLogDistribution(const Function& aDistrib, double aRand, double& aOutputX, double& aOutputIntegral, double xMin, double xMax, double aRelError);

    bool SampleDistribution(const Function& aDistrib, mcray::Randomizer& aRandomizer,  double& aOutputX, double& aOutputIntegral, double xMin, double xMax, double aRelError, int max_recurse=3);
    bool SampleLogDistribution(const Function& aDistrib, mcray::Randomizer& aRandomizer, double& aOutputX, double& aOutputIntegral, double xMin, double xMax, double aRelError, int max_recurse=3);


//bool SampleLogDistributionNegative(const Function& aDistrib, Randomizer& aRandomizer, double& aOutputX, double& aOutputIntegral,
//                                       double xMin, double xMax, double aRelError, size_t limit=1000, int key=GSL_INTEG_GAUSS15);


    template<typename X> static void RelAccuracy(X& aOutput);
	static void SetLogger(IFunctionCallHandlerX<double>* aLogger) { fLogger = aLogger; }
	static int UnitTest();
private:
	static double GslProxySampleLogscaleDistributionFunc (double x, void * params)
	{
		x=exp(x);
		const Function* f = (const Function*)params;
		return x*f->f(x);
	}
	gsl_integration_workspace* gslQAGintegrator;
	static IFunctionCallHandlerX<double>* fLogger;
};

class GslProxyFunction : public Function
{
public:
	GslProxyFunction(gsl_function aFunction, double aXmin=-DBL_MAX, double aXmax=DBL_MAX) :
		fFunction(aFunction),
		fXmin(aXmin),
		fXmax(aXmax){};
	virtual double f(double _x) const
	{
		return fFunction.function(_x,fFunction.params);
	}
	virtual double Xmin() const {return fXmin;}
	virtual double Xmax() const {return fXmax;}
	Function* Clone() const
	{//There is no way to clone fFunction.params
		ASSERT(0);
		Exception::Throw("GslProxyFunction doesn't support Clone");
		return 0;
	}
private:
	gsl_function fFunction;
	double fXmin;
	double fXmax;
};

template<typename X = double >
class ConstFunctionX : public FunctionX<X>
{
public:
	ConstFunctionX(X aValue, X aXmin=-std::numeric_limits<X>::max(), X aXmax=std::numeric_limits<X>::max()) :
		fValue(aValue),
		fXmin(aXmin),
		fXmax(aXmax){};
	virtual X f(X _x) const
	{
		return fValue;
	}
	virtual X Xmin() const {return fXmin;}
	virtual X Xmax() const {return fXmax;}
	FunctionX<X>* Clone() const
	{
		return new ConstFunctionX<X>(fValue, fXmin, fXmax);
	}
private:
	X fValue;
	X fXmin;
	X fXmax;
};

typedef ConstFunctionX<double> ConstFunction;

} /* namespace Utils */
#endif /* MATHUTILS_H_ */
