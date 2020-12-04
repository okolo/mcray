/*
 * TableFunction.h
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


#ifndef TABLEFUNCTION_H_
#define TABLEFUNCTION_H_

#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "MathUtils.h"
#include "TableReader.h"


namespace Utils {

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
	static int FindLeftX(const std::vector<X>&, X xValue);
	TableFunctionX(const std::vector<X>& _x, const std::vector<X>& _y, X _leftVal=0., X _rightVal=0.);

	//takes ownership of aXscale and aYscale
	TableFunctionX(const std::vector<X>& _x, const std::vector<X>& _y, IScale<X>* aXscale, IScale<X>* aYscale, X _leftVal=0., X _rightVal=0.);

	//stores smart ref to aReader
	TableFunctionX(TableReaderX<X>* aReader, X _leftVal=0., X _rightVal=0.);

	bool InTableRange(X aArg) const;

	virtual X f(X _x) const;
	virtual X f_scaled(X _x) const = 0;
	virtual X Xmin() const;
	virtual X Xmax() const;

	void SetAutoLimits();
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
	GSLTableFunc(TableReader *aReader, const gsl_interp_type * aInterpType=gsl_interp_linear, double _leftVal=0., double _rightVal=0.):
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

template<typename X = double >
class MatrixFunctionX : public Function2X<X>
{
public:
	MatrixFunctionX(const std::string& aFile)
	{
		std::ifstream dataFile(aFile.c_str());
		if(!dataFile)
			Exception::Throw("Failed to open " + aFile);
		std::string header;
		std::getline(dataFile, header);
		int logscX,logscY;
		const char* headerFormat = "# minX=%lg stepX=%lg logscaleX=%d minY=%lg stepY=%lg logscaleY=%d";
		if(sscanf(header.c_str(),headerFormat,
				  &xMin,&xStep,&logscX,&yMin,&yStep,&logscY)!=6)
			Exception::Throw("Invalid header format in " + aFile + "\nExpected format: " + headerFormat);
		logScaleX = (bool)logscX;
		logScaleY = (bool)logscY;
		if((logScaleX && (xMin<=0 || xStep<=1.))||xStep<=0.)
			Exception::Throw("Invalid X scale in " + aFile);
		if((logScaleY && (yMin<=0 || yStep<=1.))||yStep<=0.)
			Exception::Throw("Invalid Y scale in " + aFile);
		fData = new TableReaderX<X>(dataFile, TableReader::Auto);
		nCols = fData->numberOfColumns();
		if(logScaleX) {
			xMax = xMin*pow(xStep,nCols-1);
			xStep = log(xStep);//convert multiplier to log step
		}
		else{
			xMax = xMin + xStep*(nCols-1);
		}
		nRows = fData->getColumn(0).size();
		if(logScaleY){
			yMax = yMin*pow(yStep,nRows-1);
			yStep = log(yStep);//convert multiplier to log step
		}
		else{
			yMax = yMin + yStep*(nRows-1);
		}
	}

	X f(double x, double y) const
	{
		double binX, binY;

		if(logScaleX)
		{//log scale
			binX = log(x/xMin)/xStep;
		}
		else
		{//linear scale
			binX = (x-xMin)/xStep;
		}
		if(binX<0. || binX>nCols-1)
			return 0.;
		if(logScaleY)
		{//log scale
			binY = log(y/yMin)/yStep;
		}
		else
		{//linear scale
			binY = (y-yMin)/yStep;
		}
		if(binY<0. || binY>nRows-1)
			return 0.;
		int iX = floor(binX);
		int iY = floor(binY);
		double weightX = 1.-binX+iX;
		double weightY = 1.-binY+iY;

		//linear f interpolation //TODO implement logscale f interpolation
		double f_11 = fData->getColumn(iX)[iY]*weightX*weightY;
		double f_21 = (weightX<1.)? (fData->getColumn(iX+1)[iY]*(1.-weightX)*weightY) : 0.;
		double f_12 = (weightY<1.)? (fData->getColumn(iX)[iY+1]*weightX*(1.-weightY)) : 0.;
		double f_22 = (weightX<1.&&weightY<1.)?(fData->getColumn(iX+1)[binY+1]*(1.-weightX)*(1.-weightY)) : 0.;
		double result = f_11 + f_21 + f_12 + f_22;
		return result;
	}

	virtual X MinArg(int aArgNo) const {return aArgNo==1?xMin:yMin;};
	virtual X MaxArg(int aArgNo) const {return aArgNo==1?xMax:yMax;}

private:
	SafePtr<TableReaderX<X> >  fData;
	X xMin;
	X yMin;
	X xMax;
	X yMax;
	X xStep;
	X yStep;
	int nCols;
	int nRows;
	bool logScaleX;
	bool logScaleY;
};

typedef MatrixFunctionX<double> MatrixFunction;

} /* namespace Utils */

#endif /* TABLEFUNCTION_H_ */
