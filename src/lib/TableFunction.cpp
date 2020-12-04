/*
 * TableFunction.cpp
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


#include "PrecisionTests.h"
#include "TableFunction.h"

namespace Utils {


template<typename X> X LogScaleX<X>::ToScale(X aX)
{
	ASSERT(aX>=0);
	return aX>0?log(aX):-1000;
}

template<typename X> X LogScaleX<X>::FromScale(X aXprime)
{
	return exp(aXprime);
}

template<typename X> IScale<X>* LogScaleX<X>::Clone() const
{
	return new LogScaleX();
}

template class LogScaleX<double>;
template class LogScaleX<long double>;

template<typename X> void TableFunctionX<X>::Init()
{
	ASSERT(m_x.size()==m_y.size() && m_x.size()>1);
	fXmin = m_x[0];
	fXmax = m_x[m_x.size()-1];
}

template<typename X> TableFunctionX<X>::TableFunctionX(const std::vector<X>& _x, const std::vector<X>& _y, X _leftVal, X _rightVal):
m_x(_x),
m_y(_y),
m_leftVal(_leftVal),
m_rightVal(_rightVal)
{
	Init();
}

template<typename X> TableFunctionX<X>::TableFunctionX(const std::vector<X>& _x, const std::vector<X>& _y, IScale<X>* aXscale, IScale<X>* aYscale, X _leftVal, X _rightVal):
fCloneX(_x),
fCloneY(_y),
m_x(fCloneX),
m_y(fCloneY),
m_leftVal(_leftVal),
m_rightVal(_rightVal),
fXscale(aXscale),
fYscale(aYscale)
{
	Init();
	int iMax = _x.size();
	for(int i=0; i<iMax; i++)
	{
		if(aXscale)
			fCloneX[i] = aXscale->ToScale(_x[i]);
		if(aYscale)
			fCloneY[i] = aYscale->ToScale(_y[i]);
	}
}

template<typename X> TableFunctionX<X>::TableFunctionX(TableReaderX<X>* aReader, X _leftVal, X _rightVal):
		m_x(aReader->getColumn(0)),
		m_y(aReader->getColumn(1)),
		m_leftVal(_leftVal),
		m_rightVal(_rightVal),
		fTable(aReader)
{
	Init();
}


template<typename X> TableFunctionX<X>::TableFunctionX(const TableFunctionX<X>& aTableFunction):
fCloneX(aTableFunction.m_x),
fCloneY(aTableFunction.m_y),
m_x(fCloneX),
m_y(fCloneY),
m_leftVal(aTableFunction.m_leftVal),
m_rightVal(aTableFunction.m_rightVal),
fXscale(aTableFunction.fXscale==0 ? 0 : aTableFunction.fXscale->Clone()),
fYscale(aTableFunction.fYscale==0 ? 0 : aTableFunction.fYscale->Clone()),
fXmin(aTableFunction.fXmin),
fXmax(aTableFunction.fXmax)
{

}

template<typename X> int TableFunctionX<X>::FindLeftX(const std::vector<X>&aSet, X xValue)
{
	int right = aSet.size() - 1;
	ASSERT(right>=1);
	if((xValue > aSet[right]) || xValue < aSet[0])
		return -1;
	int left = 0;

	for(int i=(left+right)/2;right-left>1;i=(left+right)/2)
	{
		if(xValue > aSet[i])
			left=i;
		else
			right=i;
	}
	ASSERT((right - left) == 1);
	return left;
}

template<typename X> int TableFunctionX<X>::FindLeftX(X xValue) const
{
	return FindLeftX(m_x, xValue);
}

template<typename X> X TableFunctionX<X>::f(X _x) const
{
	if(_x<fXmin)
		return m_leftVal;
	if(_x>fXmax)
		return m_rightVal;
	X scaledX = fXscale ? fXscale->ToScale(_x) : _x;
	X scaledY = f_scaled(scaledX);
	return fYscale ? fYscale->FromScale(scaledY) : scaledY;
}

template<typename X> void TableFunctionX<X>::SetAutoLimits()
{
	m_leftVal = fYscale ? fYscale->FromScale(m_y[0]) : m_y[0];
	m_rightVal = fYscale ? fYscale->FromScale(m_y[m_y.size()-1]) : m_y[m_y.size()-1];
}

template<typename X> X TableFunctionX<X>::Xmin() const
{
	return m_leftVal == 0. ? fXmin : -DBL_MAX;
}

template<typename X> X TableFunctionX<X>::Xmax() const
{
	return m_rightVal == 0. ? fXmax : DBL_MAX;
}

template<typename X> bool TableFunctionX<X>::InTableRange(X aArg) const
{
	return aArg>=fXmin && aArg<=fXmax;
}

template<typename X> void TableFunctionX<X>::Print(std::ostream& aOut, X aXcoef, X aYcoef)
{
	int iMax = m_x.size();
	aOut << "# " << "leftVal=" << aYcoef*m_leftVal << " ; rightVal=" << aYcoef*m_rightVal << " ; Xmin=" << aXcoef*Xmin() << " ; Xmax=" << aXcoef*Xmax() << "\n";
	for(int i=0; i<iMax; i++)
	{
		X xScaled = m_x[i];
		X yScaled = m_y[i];
		X x = fXscale ? fXscale->FromScale(xScaled) : xScaled;
		X y = fYscale ? fYscale->FromScale(yScaled) : yScaled;
		aOut << aXcoef*x << "\t" << aYcoef*y << "\t" << xScaled << "\t" << yScaled << "\n";
	}
}

template class TableFunctionX<double>;
template class TableFunctionX<long double>;

template<typename X> X LinearFuncX<X>::f_scaled(X _x) const
{
	int i = TableFunctionX<X>::FindLeftX(_x);

	if(i<0)
		return (_x < this->m_x[0])?this->m_leftVal:this->m_rightVal;

	X x1 = this->m_x[i];
	X y1 = this->m_y[i];
	X y2 = this->m_y[i+1];
	X x2 = this->m_x[i+1];

	return y1+(y2-y1)/(x2-x1)*(_x - x1);
}

template class LinearFuncX<double>;
template class LinearFuncX<long double>;

GSLTableFunc::~GSLTableFunc()
{
	if(fSpline)
		gsl_spline_free (fSpline);
	if(fAcc)
		gsl_interp_accel_free (fAcc);
}

double GSLTableFunc::f_scaled(double _x) const
{
	return gsl_spline_eval (fSpline, _x, fAcc);
}

void GSLTableFunc::Init(const gsl_interp_type * aInterpType)
{
    fAcc = gsl_interp_accel_alloc ();
    fSpline = gsl_spline_alloc (aInterpType, m_x.size());
    try{
    gsl_spline_init (fSpline, &m_x[0], &m_y[0], m_x.size());
    }
    catch(Exception* ex)
    {
    	std::cerr << "GSLTableFunc::Init() error\n input data:\n";
    	Print(std::cerr);
    	throw ex;
    }
}



#ifdef USE_MPFR

template<> mpfr::mpreal LogScaleX<mpfr::mpreal>::FromScale(mpfr::mpreal aXprime)
{
	return mpfr::exp(aXprime);
}

template<> mpfr::mpreal LogScaleX<mpfr::mpreal>::ToScale(mpfr::mpreal aX)
{
	ASSERT(aX>=0);
	return aX>0?mpfr::log(aX):-1000;
}

template class TableFunctionX<mpfr::mpreal>;
template class LinearFuncX<mpfr::mpreal>;

#endif

} /* namespace Utils */
