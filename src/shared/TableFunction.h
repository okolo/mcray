/*
 * TableFunction.h
 *
 *  Created on: Sep 29, 2011
 *      Author: ok
 */

#ifndef TABLEFUNCTION_H_
#define TABLEFUNCTION_H_

#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "MathUtils.h"
#include "TableReader.h"


namespace Utils {

enum ELimits
{
	LimitNone = 0,
	LimitLeft = 1,
	LimitRight = 2,
	LimitAll = 3
};

template<typename X = double >
class IScale
{
public:
	virtual X ToScale(X aX) = 0;
	virtual X FromScale(X aXprime) = 0;
	virtual IScale<X>* Clone() const = 0;
	virtual ~IScale(){}
};

template<typename X = double >
class LogScaleX : public IScale<X>
{
public:
	X ToScale(X aX);
	X FromScale(X aXprime);
	IScale<X>* Clone() const;
};

typedef LogScaleX<double> LogScale;

template<typename X = double >
class TableFunctionX : public SmartReferencedObj, virtual public FunctionX<X>{
public:
	TableFunctionX(const std::vector<X>& _x, const std::vector<X>& _y, X _leftVal=0., X _rightVal=0.);

	//takes ownership of aXscale and aYscale
	TableFunctionX(const std::vector<X>& _x, const std::vector<X>& _y, IScale<X>* aXscale, IScale<X>* aYscale, X _leftVal=0., X _rightVal=0.);

	//stores smart ref to aReader
	TableFunctionX(TableReaderX<X>* aReader, X _leftVal=0., X _rightVal=0.);

	//takes ownership of aXscale and aYscale, deletes aReader if it has no zero reference count
	TableFunctionX(TableReaderX<X>* aReader, IScale<X>* aXscale, IScale<X>* aYscale, X _leftVal=0., X _rightVal=0.);

	virtual X f(X _x) const;
	virtual X f_scaled(X _x) const = 0;
	virtual X Xmin() const;
	virtual X Xmax() const;

	void SetAutoLimits(ELimits aLimits = LimitAll);
	void Print(std::ostream& aOut, X aXcoef=1., X aYcoef=1.);
private:
	//common part of constructors
	void Init();
protected:
	int FindLeftX(X xValue) const;
	//used for cloning
	TableFunctionX(const TableFunctionX<X>& aTableFunction);

	std::vector<X>			fCloneX;
	std::vector<X>			fCloneY;
	const std::vector<X>&	m_x;
	const std::vector<X>&	m_y;
	X						m_leftVal;
	X						m_rightVal;
	SafePtr<IScale<X> >		fXscale;
	SafePtr<IScale<X> >		fYscale;
	SmartPtr<TableReaderX<X> >	fTable;
	X						fXmin;
	X						fXmax;
};

typedef TableFunctionX<double> TableFunction;

template<typename X = double >
class LinearFuncX  : public TableFunctionX<X>
{
public:
	LinearFuncX(const std::vector<X>& _x,const  std::vector<X>& _y, X _leftVal=0., X _rightVal=0.):
	  TableFunctionX<X>(_x, _y, _leftVal, _rightVal)  {};

	//stores smart ref to aReader
	LinearFuncX(TableReaderX<X>* aReader, X _leftVal=0., X _rightVal=0.):
		TableFunctionX<X>(aReader,_leftVal,_rightVal) {};

	//takes ownership of aXscale and aYscale
	LinearFuncX(const std::vector<X>& _x,const std::vector<X>& _y, IScale<X>* aXscale, IScale<X>* aYscale, X _leftVal=0., X _rightVal=0.):
		TableFunctionX<X>(_x, _y, aXscale, aYscale, _leftVal, _rightVal) {}

	//takes ownership of aXscale and aYscale, deletes aReader if it has no zero reference count
	LinearFuncX(TableReaderX<X>* aReader, IScale<X>* aXscale, IScale<X>* aYscale, X _leftVal=0., X _rightVal=0.):
		TableFunctionX<X>(aReader, aXscale, aYscale, _leftVal, _rightVal) {}

	X f_scaled(X _x) const;
	FunctionX<X>* Clone() const { return new LinearFuncX<X>(*this); }
protected:
	LinearFuncX(const TableFunctionX<X>& aTableFunction):
		TableFunctionX<X>(aTableFunction) {};
};

typedef LinearFuncX<double> LinearFunc;

class GSLTableFunc  : public TableFunction
{
public:
	GSLTableFunc(const Vector& _x, const Vector& _y, const gsl_interp_type * aInterpType=gsl_interp_linear,
			double _leftVal=0., double _rightVal=0.):
	  TableFunction(_x, _y, _leftVal, _rightVal),fAcc(0),fSpline(0) {Init(aInterpType);}

	//stores smart ref to aReader
	GSLTableFunc(TableReader* aReader, const gsl_interp_type * aInterpType=gsl_interp_linear, double _leftVal=0., double _rightVal=0.):
		TableFunction(aReader,_leftVal,_rightVal),fAcc(0),fSpline(0) {Init(aInterpType);}

	//takes ownership of aXscale and aYscale
	//_x and _y are values are modified with scale functions given by aXscale and aYscale parameters
	GSLTableFunc(const Vector& _x, const Vector& _y, IScale<double>* aXscale, IScale<double>* aYscale, const gsl_interp_type * aInterpType=gsl_interp_linear, double _leftVal=0., double _rightVal=0.):
		TableFunction(_x, _y, aXscale, aYscale, _leftVal, _rightVal) {Init(aInterpType);}

	virtual ~GSLTableFunc();

	double f_scaled(double _x) const;
	Function* Clone() const { return new GSLTableFunc(*this); }
protected:
	GSLTableFunc(const GSLTableFunc& aTableFunction):
		TableFunction(aTableFunction) {Init(aTableFunction.fSpline->interp->type);};
private:
	void Init(const gsl_interp_type * aInterpType);
    gsl_interp_accel *fAcc;
    gsl_spline *fSpline;
    const gsl_interp_type *fInterpType;
};

} /* namespace Utils */

#endif /* TABLEFUNCTION_H_ */
