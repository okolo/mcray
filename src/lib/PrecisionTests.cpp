/*
 * PrecisionTests.cpp
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

/*
 * PrecisionTests.cpp
 *
 *  Created on: Dec 25, 2011
 *      Author: ok
 */

#include "PrecisionTests.h"
#include "../external/mpfrc++/mpreal.h"
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
