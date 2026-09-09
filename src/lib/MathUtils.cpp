/*
 * MathUtils.cpp
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

#define DISABLE_NR // don't use Numerical Recipes API

#ifdef DISABLE_NR

//#include <boost/numeric/odeint/config.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>

#else

#include "nr/nr3.h"
#include "nr/stepper.h"
#include "nr/odeint.h"
#include "nr/stepperdopr5.h"

#endif //#ifndef DISABLE_NR

#include "Utils.h"
#include "MathUtils.h"
#include <gsl/gsl_errno.h>
#include "gsl/gsl_sf_erf.h"
#include <gsl/gsl_math.h>
#include "TableFunction.h"


namespace Utils {
#ifdef DISABLE_NR
	using namespace boost::numeric::odeint;

	template< class Obj , class Mem >
	class ode_wrapper
	{
		Obj* m_pObj;
		Mem m_mem;

	public:

		ode_wrapper( Obj* obj , Mem mem ) : m_pObj(obj ) , m_mem(mem ) { }

		template< class State , class Deriv , class Time >
		void operator()( const State &x , Deriv &dxdt , Time t )
		{
			((*m_pObj).*m_mem)(x , dxdt , t );
		}
	};

	template< class Obj , class Mem >
	ode_wrapper< Obj , Mem > make_ode_wrapper( Obj* obj , Mem mem )
	{
		return ode_wrapper< Obj , Mem >( obj , mem );
	}


	template< class Obj , class Mem >
	class observer_wrapper
	{
		Obj* m_pObj;
		Mem m_mem;

	public:

		observer_wrapper( Obj* obj , Mem mem ) : m_pObj(obj ) , m_mem(mem ) { }

		template< class State , class Time >
		void operator()( const State &x , Time t )
		{
			((*m_pObj).*m_mem)(x , t );
		}
	};

	template< class Obj , class Mem >
	observer_wrapper< Obj , Mem > make_observer_wrapper( Obj* obj , Mem mem )
	{
		return observer_wrapper< Obj , Mem >( obj , mem );
	}

	template< class X=double >
	class Sampler4 : public ISampler<X> {
		//typedef runge_kutta_dopri5<X> dopri5_type;
		//typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
		//typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
		//typedef bulirsch_stoer_dense_out<X> dopri5_type;
	public:
		Sampler4(){
			fIntermediateVals=new std::vector<X>;
			fIntermediateTs=new std::vector<X>;
			fIntermediateTs->reserve(256);
			fIntermediateVals->reserve(256);
		}

		~Sampler4(){
			delete fIntermediateVals;
			delete fIntermediateTs;
		}

		void System(const X &x , X &dxdt , X t ){
			dxdt=fFunc->f(t);
		}

		void Output(const X &x , X t ){
			if(fMaxX>0 && x>=fMaxX){
				X lastT = *(fIntermediateTs->end()-1);
				X lastX = *(fIntermediateVals->end()-1);
				if((t-lastT)/(t+lastT)<fRelErr && t-lastT<fAbsErr){
					X result = lastT + (t-lastT)/(x-lastX)*(fMaxX-lastX);
					throw result;
				}
				X totRate = x-lastX;
				throw sample(*fFunc, lastT, t, (fMaxX-lastX)/(x-lastX), totRate, fRelErr, fAbsErr);
			}
			else {
				fIntermediateVals->push_back(x);
				fIntermediateTs->push_back(t);
			}
		}

		/// If total rate is known it should be passed via aTotRate argument
		/// Otherwize aTotRate should be set to 0 (it will be calculated by function)
		X sample(const FunctionX<X>& f, X aTmin, X aTmax, X aRand,
				 X & aTotRate, X aRelErr, X aAbsErr = 1e300){
			fMaxX = aTotRate*aRand;
			fRelErr = aRelErr;
			fAbsErr = aAbsErr;

			size_t nPoints = (size_t)((aTmax-aTmin)/aRelErr+0.5)/2+1;
			X step= (aTmax-aTmin)/nPoints;

			fIntermediateTs->resize(0);//this will not decrease capacity of vectors
			fIntermediateVals->resize(0);

			X curT=aTmin;
			fFunc = &f;
			aTotRate = 0.;

			X dt = (aTmax - aTmin) * aRelErr;

			try {
				//dense_output_dopri5_type dopri5 = make_dense_output( aAbsErr , aRelErr , dopri5_type() );
				bulirsch_stoer_dense_out< X > stepper( aAbsErr , aRelErr);
				integrate_adaptive( stepper , make_ode_wrapper(this,&Sampler4::System) , aTotRate , aTmin , aTmax , dt , make_observer_wrapper(this, &Sampler4::Output) );
				//integrate_const
				//integrate_adaptive(dopri5, make_ode_wrapper(this,&Sampler3::System), aTotRate, aTmin, aTmax, dt,
								   //make_observer_wrapper(this, &Sampler3::Output));
			}catch (X aValue){
				return aValue;
			}
			if(aTotRate<=0){
				return aTmax;
			}
			std::vector<X>& intermediateVals = *fIntermediateVals;
			std::vector<X>& intermediateTs = *fIntermediateTs;
			X searchVal = aTotRate * aRand;
			size_t i2= intermediateVals.size();
			size_t i1=0;
			for(size_t i=(i2+i1)/2; i2-i1>1; i=(i2+i1)/2)
			{
				if(intermediateVals[i] < searchVal)
					i1 = i;
				else
					i2 = i;
			}
			X t1 = intermediateTs[i1];
			X t2 = intermediateTs[i2];
			X x1 = intermediateVals[i1];
			X x2 = intermediateVals[i2];
			if((t2-t1)/(t2+t1)<fRelErr && t2-t1<fAbsErr){
				return t1 + (t2-t1)/(x2-x1)*(searchVal-x1);
			}
			X totRate=x2-x1;
			return sample(*fFunc, t1, t2, (searchVal-x1)/(x2-x1), totRate, fRelErr, fAbsErr);
		}
	private:
		const FunctionX<X>* fFunc;
		std::vector<X>* fIntermediateVals;
		std::vector<X>* fIntermediateTs;
		X      fMaxX;
		X      fRelErr;
		X      fAbsErr;
	};

	template< class X=double >
	class Sampler3 : public ISampler<X> {
		typedef runge_kutta_dopri5<X> dopri5_type;
		typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
		typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
	public:
		Sampler3(){
			fIntermediateVals=new std::vector<X>;
			fIntermediateTs=new std::vector<X>;
			fIntermediateTs->reserve(256);
			fIntermediateVals->reserve(256);
		}

		~Sampler3(){
			delete fIntermediateVals;
			delete fIntermediateTs;
		}

		void System(const X &x , X &dxdt , X t ){
			dxdt=fFunc->f(t);
		}

		void Output(const X &x , X t ){
			if(fMaxX>0 && x>=fMaxX){
				X lastT = *(fIntermediateTs->end()-1);
				X lastX = *(fIntermediateVals->end()-1);
				if((t-lastT)/(t+lastT)<fRelErr && t-lastT<fAbsErr){
					X result = lastT + (t-lastT)/(x-lastX)*(fMaxX-lastX);
					throw result;
				}
				X totRate = x-lastX;
				throw sample(*fFunc, lastT, t, (fMaxX-lastX)/(x-lastX), totRate, fRelErr, fAbsErr);
			}
			else {
				fIntermediateVals->push_back(x);
				fIntermediateTs->push_back(t);
			}
		}

		/// If total rate is known it should be passed via aTotRate argument
		/// Otherwize aTotRate should be set to 0 (it will be calculated by function)
		X sample(const FunctionX<X>& f, X aTmin, X aTmax, X aRand,
				 X & aTotRate, X aRelErr, X aAbsErr = 1e300){
			fMaxX = aTotRate*aRand;
			fRelErr = aRelErr;
			fAbsErr = aAbsErr;

			size_t nPoints = (size_t)((aTmax-aTmin)/aRelErr+0.5)/2+1;
			X step= (aTmax-aTmin)/nPoints;

			fIntermediateTs->resize(0);//this will not decrease capacity of vectors
			fIntermediateVals->resize(0);

			X curT=aTmin;
			fFunc = &f;
			aTotRate = 0.;
			dense_output_dopri5_type dopri5 = make_dense_output( aAbsErr , aRelErr , dopri5_type() );
			X dt = (aTmax - aTmin) * aRelErr;

			try {
				integrate_adaptive(dopri5, make_ode_wrapper(this,&Sampler3::System), aTotRate, aTmin, aTmax, dt,
								   make_observer_wrapper(this, &Sampler3::Output));
			}catch (X aValue){
				return aValue;
			}
			if(aTotRate<=0){
				return aTmax;
			}
			std::vector<X>& intermediateVals = *fIntermediateVals;
			std::vector<X>& intermediateTs = *fIntermediateTs;
			X searchVal = aTotRate * aRand;
			size_t i2= intermediateVals.size();
			size_t i1=0;
			for(size_t i=(i2+i1)/2; i2-i1>1; i=(i2+i1)/2)
			{
				if(intermediateVals[i] < searchVal)
					i1 = i;
				else
					i2 = i;
			}
			X t1 = intermediateTs[i1];
			X t2 = intermediateTs[i2];
			X x1 = intermediateVals[i1];
			X x2 = intermediateVals[i2];
			if((t2-t1)/(t2+t1)<fRelErr && t2-t1<fAbsErr){
				return t1 + (t2-t1)/(x2-x1)*(searchVal-x1);
			}
			X totRate=x2-x1;
			return sample(*fFunc, t1, t2, (searchVal-x1)/(x2-x1), totRate, fRelErr, fAbsErr);
		}
	private:
		const FunctionX<X>* fFunc;
		std::vector<X>* fIntermediateVals;
		std::vector<X>* fIntermediateTs;
		X      fMaxX;
		X      fRelErr;
		X      fAbsErr;
	};

	template<typename X = double >
	class Sampler : public ISampler<X> {
		typedef runge_kutta_dopri5<X> dopri5_type;
		typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
		typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
	public:
		Sampler(){
			fIntermediateVals=0;
			fIntermediateTs=0;
		}
		void operator()(const X &x , X &dxdt , const X t ){
			dxdt=fFunc->f(t);
		}
		void operator()(const X &x , const X t ){
			(*fIntermediateVals)[fLastStep++]=x;
			double a=(*fIntermediateVals)[fLastStep-1];
			a=a+1;
		}

		X sample(const FunctionX<X>& f, X aTmin, X aTmax, X aRand,
				 X & aTotRate, X aRelErr, X aAbsErr = 1e300){
			// create a vector with observation time points
			size_t nPoints = (size_t)((aTmax-aTmin)/aRelErr/2+0.5);
			X step= (aTmax-aTmin)/nPoints;
			std::vector<X> intermediateVals(nPoints + 1);
			std::vector<X> intermediateTs(nPoints + 1);
			fIntermediateVals=&intermediateVals;
			fIntermediateTs=&intermediateTs;
			X curT=aTmin;
			for( size_t i=0 ; i<=nPoints ; ++i, curT+=step )
				intermediateTs[i] = curT;

			fFunc = &f;
			aTotRate = 0.;
			dense_output_dopri5_type dopri5 = make_dense_output( aAbsErr , aRelErr , dopri5_type() );
			X dt = (aTmax - aTmin) * aRelErr;

			fLastStep=0;
			integrate_times(dopri5 , (*this) , aTotRate , intermediateTs, dt , (*this) );
			if(aTotRate==0.){
				//aTotRate=0.;
				return aTmax;
			}

			X searchVal = aTotRate * aRand;
			size_t i2= intermediateVals.size();
			size_t i1=0;
			for(size_t i=(i2+i1)/2; i2-i1>1; i=(i2+i1)/2)
			{
				if(intermediateVals[i] < searchVal)
					i1 = i;
				else
					i2 = i;
			}
			X t1 = intermediateTs[i1];
			X t2 = intermediateTs[i2];
			X x1 = intermediateVals[i1];
			X x2 = intermediateVals[i2];
			return t1 + (t2 - t1) / (x2 - x1) * (searchVal - x1);
		}

		const FunctionX<X>* fFunc;
		std::vector<X>* fIntermediateVals;
		std::vector<X>* fIntermediateTs;

		size_t      fLastStep;
	};

	template<typename X = double >
	class LogSampler : public ISampler<X> {
		typedef runge_kutta_dopri5<X> dopri5_type;
		typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
		typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
	public:
		LogSampler(){
			fIntermediateVals=0;
			fIntermediateTs=0;
		}
		void operator()(const X &x , X &dxdt , const X t ){
			dxdt=fFunc->f(t);
		}
		void operator()(const X &x , const X t ){
			(*fIntermediateVals)[fLastStep++]=x;
			double a=(*fIntermediateVals)[fLastStep-1];
			a=a+1;
		}

		X sample(const FunctionX<X>& f, X aTmin, X aTmax, X aRand,
				 X & aTotRate, X aRelErr, X aAbsErr = 1e300){
			// create a vector with observation time points
			size_t nPoints = (size_t)(log(aTmax/aTmin)/log(1.0+aRelErr)/2+0.5);
			X step= pow(aTmax/aTmin, 1./nPoints);
			std::vector<X> intermediateVals(nPoints + 1);
			std::vector<X> intermediateTs(nPoints + 1);
			fIntermediateVals=&intermediateVals;
			fIntermediateTs=&intermediateTs;
			X curT=aTmin;
			for( size_t i=0 ; i<=nPoints ; ++i, curT*=step )
				intermediateTs[i] = curT;

			fFunc = &f;
			aTotRate = 0.;
			dense_output_dopri5_type dopri5 = make_dense_output( aAbsErr , aRelErr , dopri5_type() );
			X dt = (aTmax - aTmin) * aRelErr;

			fLastStep=0;
			integrate_times(dopri5 , (*this) , aTotRate , intermediateTs, dt , (*this) );
			if(aTotRate==0.){
				//aTotRate=0.;
				return aTmax;
			}

			X searchVal = aTotRate * aRand;
			size_t i2= intermediateVals.size();
			size_t i1=0;
			for(size_t i=(i2+i1)/2; i2-i1>1; i=(i2+i1)/2)
			{
				if(intermediateVals[i] < searchVal)
					i1 = i;
				else
					i2 = i;
			}
			X t1 = intermediateTs[i1];
			X t2 = intermediateTs[i2];
			X x1 = intermediateVals[i1];
			X x2 = intermediateVals[i2];
			return t1 + (t2 - t1) / (x2 - x1) * (searchVal - x1);
		}

		const FunctionX<X>* fFunc;
		std::vector<X>* fIntermediateVals;
		std::vector<X>* fIntermediateTs;

		size_t      fLastStep;
	};

	template<typename X = double >
	class Sampler2 : public ISampler<X> {
		typedef runge_kutta_dopri5<X> dopri5_type;
		typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
		typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;
	public:
		Sampler2(){
			fIntermediateVals=0;
			fIntermediateTs=0;
		}
		~Sampler2()
		{
			fIntermediateVals=0;
			fIntermediateTs=0;
		}
		void operator()(const X &x , X &dxdt , const double t ){
			dxdt=fFunc->f(t);
		}
		void operator()(const X &x , const X t ){
			if(fMaxX>0 && x>=fMaxX){
				X lastT = *(fIntermediateTs->end()-1);
				X lastX = *(fIntermediateVals->end()-1);
				if((t-lastT)/(t+lastT)<fRelErr && t-lastT<fAbsErr){
					X result = lastT + (t-lastT)/(x-lastX)*(fMaxX-lastX);
					throw result;
				}
				X totRate = x-lastX;
				throw sample(*fFunc, lastT, t, (fMaxX-lastX)/(x-lastX), totRate, fRelErr, fAbsErr);
			}
			else {
				fIntermediateVals->push_back(x);
				fIntermediateTs->push_back(t);
			}
		}

		/// If total rate is known it should be passed via aTotRate argument
		/// Otherwize aTotRate should be set to 0 (it will be calculated by function)
		X sample(const FunctionX<X>& f, X aTmin, X aTmax, X aRand,
				 X & aTotRate, X aRelErr, X aAbsErr = 1e300){
			fMaxX = aTotRate*aRand;
			//fTmin = aTmin;
			//fTmax = aTmax;
			//fTotRate = aTotRate;
			fRelErr = aRelErr;
			fAbsErr = aAbsErr;

			size_t nPoints = 128;//(size_t)((aTmax-aTmin)/aRelErr+0.5)/2+1;
			//X step= (aTmax-aTmin)/nPoints;
			std::vector<X> intermediateVals;
			std::vector<X> intermediateTs;
			intermediateTs.reserve(nPoints);
			intermediateVals.reserve(nPoints);
			fIntermediateVals=&intermediateVals;
			fIntermediateTs=&intermediateTs;
			X curT=aTmin;
			fFunc = &f;
			aTotRate = 0.;

			///TODO: try different steppers

			dense_output_dopri5_type dopri5 = make_dense_output( aAbsErr , aRelErr , dopri5_type() );
			X dt = (aTmax - aTmin) * aRelErr;

			try {
				integrate_adaptive(dopri5, (*this), aTotRate, aTmin, aTmax, dt,
								   (*this));//this will create two copies of (*this)
				/* //this one is a bit slower than dense_output_dopri5_type
				typedef runge_kutta_cash_karp54< X > error_stepper_type;
				integrate_adaptive( make_controlled< error_stepper_type >( aAbsErr , aRelErr ) ,
									(*this) , aTotRate , aTmin, aTmax, dt,
									(*this));*/

			}catch (X aValue){
				return aValue;
			}
			if(aTotRate<=0){
				return aTmax;
			}

			X searchVal = aTotRate * aRand;
			size_t i2= intermediateVals.size();
			size_t i1=0;
			for(size_t i=(i2+i1)/2; i2-i1>1; i=(i2+i1)/2)
			{
				if(intermediateVals[i] < searchVal)
					i1 = i;
				else
					i2 = i;
			}
			X t1 = intermediateTs[i1];
			X t2 = intermediateTs[i2];
			X x1 = intermediateVals[i1];
			X x2 = intermediateVals[i2];
			if((t2-t1)/(t2+t1)<fRelErr && t2-t1<fAbsErr){
				return t1 + (t2-t1)/(x2-x1)*(searchVal-x1);
			}
			X totRate=x2-x1;
			return sample(*fFunc, t1, t2, (searchVal-x1)/(x2-x1), totRate, fRelErr, fAbsErr);
		}

		const FunctionX<X>* fFunc;
		std::vector<X>* fIntermediateVals;
		std::vector<X>* fIntermediateTs;
		X      fMaxX;
		X      fRelErr;
		X      fAbsErr;
	};
#endif //#ifdef DISABLE_NR

	template<typename X = double >
	class NRSampler : public ISampler<X> {
		X sample(const FunctionX<X>& f, X aTmin, X aTmax, X aRand,
				 X & aTotRate, X aRelErr, X aAbsErr = 1e300){
			X x;
            MathUtils::_SampleLogDistribution(f,aRand,x,aTotRate,aTmin,aTmax,aRelErr);
		}
	};

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

	template<class X> class GaussDisr : public FunctionX<X>{
	public:
		X mean;
		X sigma;
		GaussDisr(X aMean, X aSigma):mean(aMean),sigma(aSigma){}
		virtual X f(X _x) const{
			X diff = (_x-mean)/sigma;
			return exp(-diff*diff*0.5)/sigma*0.398942280401433;
		}
	};

	template<class X> void ISampler<X>::UnitTest()
	{

		X mean = 1e10;
		X sigma = 1e9;
		X minX=1e9;
		X maxX=1e11;
		X minI = 0.5*(1+gsl_sf_erf((minX-mean)/sigma/sqrt(2.0)));
		X maxI = 0.5*(1+gsl_sf_erf((maxX-mean)/sigma/sqrt(2.0)));
		X totI=maxI-minI;

		GaussDisr<X> gd(mean,sigma);
		X relError = 1e-6;
		int nSteps = 10000;
//gnuplot command: 0.5*(1+erf((x-mean)/sigma/sqrt(2.0))) w l

		for(int i=1; i<nSteps; i++){
			X rand = 1./nSteps*i;
			X totRate = 0;
			X curX = sample(gd, minX, maxX, rand, totRate, relError);
			X exactFrac = (0.5*(1+gsl_sf_erf((curX-mean)/sigma/sqrt(2.0)))-minI)/totI;
			std::cout << curX << "\t" << rand << "\t"  << exactFrac << std::endl;
		}
	}

	template void ISampler<double>::UnitTest();
	template void ISampler<long double>::UnitTest();

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

#ifndef DISABLE_NR
class MathUtilsODE
{
public:
	MathUtilsODE(const Function& aDistrib):fDistrib(aDistrib)
	{}
	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
		double val = fDistrib(x);
		if(val!=0)
			dydx[0]= val;
		else
			dydx[0]= 0.;
		ASSERT(val > -1.6e308 && val < 1.6e308);
	}
private:
	const Function& fDistrib;
};
#endif

IFunctionCallHandlerX<double>* MathUtils::fLogger = 0;

bool MathUtils::_SampleDistribution(const Function& aDistrib, double aRand, double& aOutputX, double& aOutputIntegral, double xMin, double xMax, double aRelError)
{
#ifdef DISABLE_NR
    Sampler<double> dil;//slower 25 sec
    //Sampler2<double> dil;// 20 sec (todo: fix memory leaks)
    //Sampler3<double> dil;// 20 sec (todo: fix memory leaks)
    aOutputIntegral=0.;

    aOutputX=dil.sample(aDistrib, xMin, xMax, aRand, aOutputIntegral, aRelError);

    ASSERT_VALID_NO(aOutputX);
    return true;

#else
	const Function* distrib = &aDistrib;
	SafePtr<DebugFunctionX<double> > func;
	if(fLogger)
	{
		func = new DebugFunctionX<double>(aDistrib, *fLogger);
		distrib = func;
	}

	double distrXmin = distrib->Xmin();
	if(xMin < distrXmin)
		xMin = distrXmin;
	double distrXmax = distrib->Xmax();
	if(xMax > distrXmax)
		xMax = distrXmax;
	ASSERT(xMax>xMin);

	double initialStep = aRelError*(xMax-xMin);

	MathUtilsODE d(*distrib);

	const Doub atol=0.;//1.0e-3;
	const Doub hmin=0.0;//minimal step (can be zero)
	VecDoub ystart(1);
	ystart[0]=0.;
	Output out(-1); //output is saved at every integration step
	Odeint<StepperDopr5<MathUtilsODE> > ode(ystart,xMin,xMax,atol,aRelError,initialStep,hmin,out,d);
	ode.integrate();
	aOutputIntegral = ystart[0];
	if(aOutputIntegral<=0)
		return false;
	aRand *= aOutputIntegral;
	int i1=0;
	int i2=out.count-1;
	double* ysave = out.ysave[0];
	for(int i=(i1+i2)/2; i2-i1>1; i=(i1+i2)/2)
	{
		double y = ysave[i];
		if(y<aRand)
			i1 = i;
		else
			i2 = i;
	}
	double y1=ysave[i1];
	double y2=ysave[i2];
	double x1=out.xsave[i1];
	double x2=out.xsave[i2];
	aRand -= y1;
	aOutputX = x1 + (x2-x1)/(y2-y1)*aRand;//make a linear estimate of X
	if(fabs((x2-x1)/aOutputX)<=aRelError)
		return true;
	ystart[0]=0.;
	initialStep = aRelError*(aOutputX>0 ? aOutputX : (x2-x1));
	Output out2(2+(int)fabs((x2-x1)/initialStep));

	Odeint<StepperDopr5<MathUtilsODE> > ode2(ystart,x1,x2,atol,aRelError,initialStep,hmin,out2,d);
	ode2.integrate();
	i1=0;
	i2=out2.count-1;
	ysave = out2.ysave[0];
	for(int i=(i1+i2)/2; i2-i1>1; i=(i1+i2)/2)
	{
		double y = ysave[i];
		if(y<aRand)
			i1 = i;
		else
			i2 = i;
	}
	y1=ysave[i1];
	y2=ysave[i2];
	x1=out2.xsave[i1];
	x2=out2.xsave[i2];
	aOutputX = x1 + (x2-x1)/(y2-y1)*(aRand-y1);//make a linear estimate of X
	return true;
#endif
}

#ifndef DISABLE_NR
class MathUtilsLogODE
{
public:
	MathUtilsLogODE(const Function& aDistrib):fDistrib(aDistrib)
	{}
	void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
		double xx = exp(x);
		double val = xx*fDistrib(xx);
		if(val!=0)
			dydx[0]= val;
		else
			dydx[0]= 0.;
		ASSERT(val > -1.6e308 && val < 1.6e308);
	}
private:
	const Function& fDistrib;
};
#endif

bool MathUtils::_SampleLogDistribution(const Function& aDistrib, double aRand, double& aOutputX, double& aOutputIntegral, double xMin, double xMax, double aRelError)
{
    ASSERT(aRelError>0 && aRelError<=0.1);
    ASSERT(xMin>0 && xMax>xMin);
#ifdef DISABLE_NR
    LogSampler<double> dil;//slower 25 sec
    //Sampler2<double> dil;// 20 sec (todo: fix memory leaks)
    //Sampler3<double> dil;// 20 sec (todo: fix memory leaks)
    aOutputIntegral=0.;

    aOutputX=dil.sample(aDistrib, xMin, xMax, aRand, aOutputIntegral, aRelError);

    ASSERT_VALID_NO(aOutputX);
    return true;
#else
	const Function* distrib = &aDistrib;
	SafePtr<DebugFunctionX<double> > func;
	if(fLogger)
	{
		func = new DebugFunctionX<double>(aDistrib, *fLogger);
		distrib = func;
	}

	double distrXmin = distrib->Xmin();
	if(xMin < distrXmin)
		xMin = distrXmin;
	double distrXmax = distrib->Xmax();
	if(xMax > distrXmax)
		xMax = distrXmax;
	ASSERT(xMax>xMin && xMin>0);
	xMin=log(xMin);
	xMax=log(xMax);

	double initialStep = 0.5*(xMax-xMin);
//	if(aRelError<initialStep)
//		initialStep=aRelError;

	MathUtilsLogODE d(*distrib);

	const Doub atol=0;
	const Doub hmin=0.0;//minimal step (can be zero)
	VecDoub ystart(1);
	ystart[0]=0.;
	Output out(-1); //output is saved at every integration step
	Odeint<StepperDopr5<MathUtilsLogODE> > ode(ystart,xMin,xMax,atol,aRelError,initialStep,hmin,out,d);
	ode.integrate();
	aOutputIntegral = ystart[0];
	if(aOutputIntegral<=0)
		return false;
	double yRand = aRand*aOutputIntegral;
	int i1=0;
	int i2=out.count-1;
	double* ysave = out.ysave[0];
	for(int i=(i1+i2)/2; i2-i1>1; i=(i1+i2)/2)
	{
		double y = ysave[i];
		if(y<yRand)
			i1 = i;
		else
			i2 = i;
	}
	double y1=ysave[i1];
	double y2=ysave[i2];
	double x1=out.xsave[i1];
	double x2=out.xsave[i2];
	double yFrac = (yRand-y1)/(y2-y1);
	if((x2-x1)<aRelError)
	{
		aOutputX = x1 + (x2-x1)*yFrac;//make a linear estimate of X in log scale
		ASSERT(aOutputX>=xMin && aOutputX<=xMax);
	}
	else
	{
		ystart[0]=0.;
		initialStep = 0.5*(x2-x1);//aRelError;
		Output out2(2+(int)fabs((x2-x1)/initialStep));
		Odeint<StepperDopr5<MathUtilsLogODE> > ode2(ystart,x1,x2,atol,aRelError,initialStep,hmin,out2,d);
		ode2.integrate();
		yRand = ystart[0]*yFrac;
		i1=0;
		i2=out2.count-1;
		ysave = out2.ysave[0];
		for(int i=(i1+i2)/2; i2-i1>1; i=(i1+i2)/2)
		{
			double y = ysave[i];
			if(y<yRand)
				i1 = i;
			else
				i2 = i;
		}
		y1=ysave[i1];
		y2=ysave[i2];
		x1=out2.xsave[i1];
		x2=out2.xsave[i2];
		aOutputX = x1 + (x2-x1)/(y2-y1)*(yRand-y1);//make a linear estimate of X in log scale
		ASSERT(aOutputX>=xMin && aOutputX<=xMax);
	}
	aOutputX = exp(aOutputX);
	ASSERT_VALID_NO(aOutputX);
	return true;
#endif
}

// diagnostics
//struct SampleStat{
//     long long tot_calls;
//     long long failed_calls;
//     long long level;
//};
//static SampleStat stats_SampleDistribution;
//static SampleStat stats_SampleLogDistribution;
//
//void print_sampling_stat(SampleStat s){
//#pragma omp critical (stderr)
//    std::cerr << "sampling failures %: " << 100.0*s.failed_calls/s.tot_calls << std::endl;
//}

    bool MathUtils::SampleDistribution(const Function& aDistrib, mcray::Randomizer& aRandomizer,  double& aOutputX, double& aOutputIntegral, double xMin, double xMax,
                                       double aRelError, int max_recurse){
        size_t limit = 1000;
        int sampling_limit = 1000;
        ASSERT(aRelError > 0 && aRelError <= 0.1);

        double distrXmin = aDistrib.Xmin();
        if(xMin < distrXmin)
            xMin = distrXmin;
        double distrXmax = aDistrib.Xmax();
        if(xMax > distrXmax)
            xMax = distrXmax;
        ASSERT(xMax>xMin);

        double max_arg = 0;
        double max_val = -1e-300; //will be used for sampling
        MaxFunction pdf(aDistrib, max_arg, max_val);

        aOutputIntegral = Integration_qag (pdf,xMin,xMax,0,aRelError, limit);
// Rejection Sampling
        for(int attempts_left = sampling_limit; attempts_left>0; attempts_left--){
            double x = xMin + aRandomizer.Rand() * (xMax-xMin);
            double threshold = aRandomizer.Rand() * max_val;
            double f = aDistrib(x);
            if (f >= threshold) {
                aOutputX = x;
                return true;
            }
        }

        if (max_recurse <= 0){
            // perhaps the distribution is very narrow return argmax of PDF
            std::cerr << "SampleDistribution maximal number of iteration reached, returning argmax of PDF" << std::endl;
            aOutputX = max_arg;
            return true;
        }
        // perhaps the distribution is very narrow try to sample around its maximum
        double new_xmin = max_arg - 0.05*(xMax-xMin);
        new_xmin = (new_xmin<xMin)?xMin:new_xmin;
        double new_xmax = max_arg + 0.05*(xMax-xMin);
        new_xmax = (new_xmax>xMax)?xMax:new_xmax;
        double out_int;
        return SampleDistribution(aDistrib, aRandomizer,  aOutputX, out_int, new_xmin, new_xmax, aRelError, max_recurse-1);
}

class LogscaleDistr : public Function
{
public:
    LogscaleDistr(const Function& aOrigFunc, double& aMaxArg, double& aMaxValue):
            fOrigFunc(aOrigFunc),
            fMaxArg(aMaxArg),
            fMaxVal(aMaxValue),
            fXmin(log(aOrigFunc.Xmin())),
            fXmax(log(aOrigFunc.Xmax()))
    {};
    virtual double f(double _x) const
    {
        double orig_x = exp(_x);
        double y = orig_x * fOrigFunc.f(orig_x);
        if (y > fMaxVal){
            fMaxArg = _x;
            fMaxVal = y;
        }
        return y;
    }
    virtual double Xmin() const {return fXmin;}
    virtual double Xmax() const {return fXmax;}
    virtual ~LogscaleDistr(){};
    virtual Function* Clone() const { return new LogscaleDistr(fOrigFunc,fMaxArg,fMaxVal); }
private:
    const Function&				fOrigFunc;
    double&	                    fMaxArg;
    double&	                    fMaxVal;
    double                      fXmin;
    double                      fXmax;
};

bool MathUtils::SampleLogDistribution(const Function& aDistrib, mcray::Randomizer& aRandomizer, double& aOutputX, double& aOutputIntegral, double xMin, double xMax, double aRelError, int max_recurse){
    size_t limit = 1000;
    int sampling_limit = 1000;
    ASSERT(aRelError > 0 && aRelError <= 0.1);

    double distrXmin = aDistrib.Xmin();
    if(xMin < distrXmin)
        xMin = distrXmin;
    double distrXmax = aDistrib.Xmax();
    if(xMax > distrXmax)
        xMax = distrXmax;
    ASSERT(xMax>xMin && xMin>0);
    xMin=log(xMin);
    xMax=log(xMax);

    double max_arg = 0;
    double max_val = -1e-300; //will be used for sampling
    LogscaleDistr pdf(aDistrib, max_arg, max_val);

    aOutputIntegral = Integration_qag (pdf,xMin,xMax,0,aRelError, limit);

    // Rejection Sampling
    for(int attempts_left = sampling_limit; attempts_left>0; attempts_left--){
        double x = xMin + aRandomizer.Rand()*(xMax-xMin);
        double threshold = aRandomizer.Rand() * max_val;
        double f = pdf(x);
        if (f >= threshold) {
            aOutputX = exp(x);
            return true;
        }
    }

    if (max_recurse <= 0){
        // perhaps the distribution is very narrow return argmax of PDF
        std::cerr << "SampleLogDistribution maximal number of iteration reached, returning argmax of PDF" << std::endl;
        aOutputX = exp(max_arg);
        return true;
    }
    // perhaps the distribution is very narrow try to sample around its maximum
    double new_xmin = max_arg - 0.05*(xMax-xMin);
    new_xmin = (new_xmin<xMin)?xMin:new_xmin;
    double new_xmax = max_arg + 0.05*(xMax-xMin);
    new_xmax = (new_xmax>xMax)?xMax:new_xmax;
    new_xmin = exp(new_xmin);
    new_xmax = exp(new_xmax);
    double out_int;
    return SampleLogDistribution(aDistrib, aRandomizer,  aOutputX, out_int, new_xmin, new_xmax, aRelError, max_recurse-1);
}

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

	int MathUtils::UnitTest(){
		//LogSampler<double> dil;//slower 25 sec
		//Sampler2<double> dil;// 20 sec (todo: fix memory leaks)
		//Sampler3<double> dil;// 20 sec (todo: fix memory leaks)
		//Sampler4<double> dil;
		NRSampler<double> dil;
		dil.UnitTest();
        return 0;
	}

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
