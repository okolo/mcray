/*
 * PrecisionTests.cpp
 *
 *  Created on: Dec 25, 2011
 *      Author: ok
 */

#include "PrecisionTests.h"
#include "mpfrc++/mpreal.h"
#include "Utils.h"

using namespace mpfr;

namespace Test{

void PrecisionTests::SetPrecision(int aPrecision)
{
	mpreal::set_default_prec(aPrecision);
}

double PrecisionTests::PartialSigmaAccurate(double rMin, double s)
{
/*
Formula for partial cross section was taken from http://lanl.arxiv.org/abs/1106.5508v1 eq. (9)
*/
	mpreal yMin=1./s;//s here is in units of m^2
	mpreal yMax = 1-rMin;
	if(yMin>=yMax)
	{
		ASSERT((yMax-yMin)/(yMax+yMin)<1e-10);
		return 0.;
	}

	mpreal y_1 = 1.-yMin;
	mpreal y_1_2 = y_1*y_1;
	mpreal a=yMax-yMin;
	mpreal result = 0;

	result = yMin*(yMax-yMin)/y_1*(log(yMax/yMin)/(yMax-yMin)*(1.-4.*yMin*(1.+yMin)/y_1_2)+4.*(yMin/yMax+yMin)/y_1_2+0.5*(yMax+yMin));
	ASSERT_VALID_NO(result.toDouble());

	return result.toDouble();
}

}//namespace Test
