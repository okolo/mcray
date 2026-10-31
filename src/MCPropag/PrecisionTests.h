/*
 * PrecisionTests.h
 *
 *  Created on: Dec 25, 2011
 *      Author: ok
 */

#ifndef PRECISIONTESTS_H_
#define PRECISIONTESTS_H_

#ifdef USE_MPFR

#include "mpfrc++/mpreal.h"

namespace Test
{

class PrecisionTests {
public:
	static void SetPrecision(int aPrecision);
	static double PartialSigmaAccurate(double rMin, double s);
};

}
#endif //#ifdef USE_MPFR


#endif /* PRECISIONTESTS_H_ */
