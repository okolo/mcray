/*
 * File:   Randomizer.cpp
 * Author: ok
 *
 * Created on August 26, 2011, 7:58 PM
 */

#include "Randomizer.h"
#include <math.h>
#include "Utils.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

namespace MCPropag
{
Randomizer::Randomizer(unsigned long int aSeed):
		fRand2(0)
{
	/*
The following table shows the relative performance of a selection the available random number generators.
The fastest simulation quality generators are taus, gfsr4 and mt19937.
The generators which offer the best mathematically-proven quality are those based on the ranlux algorithm.

     1754 k ints/sec,    870 k doubles/sec, taus
     1613 k ints/sec,    855 k doubles/sec, gfsr4
     1370 k ints/sec,    769 k doubles/sec, mt19937
      565 k ints/sec,    571 k doubles/sec, ranlxs0
      400 k ints/sec,    405 k doubles/sec, ranlxs1
      490 k ints/sec,    389 k doubles/sec, mrg
      407 k ints/sec,    297 k doubles/sec, ranlux
      243 k ints/sec,    254 k doubles/sec, ranlxd1
      251 k ints/sec,    253 k doubles/sec, ranlxs2
      238 k ints/sec,    215 k doubles/sec, cmrg
      247 k ints/sec,    198 k doubles/sec, ranlux389
      141 k ints/sec,    140 k doubles/sec, ranlxd2
	*/
	fRand = gsl_rng_alloc (gsl_rng_ranlux);//gsl_rng*
	if(aSeed>0)//use aSeed=0 and GSL_RNG_SEED environment variable to overwrite seed at runtime
		gsl_rng_set ((gsl_rng*)fRand, aSeed);
}

Randomizer::~Randomizer() {
	gsl_rng_free((gsl_rng*)fRand);
	if(fRand2)
		gsl_rng_free((gsl_rng*)fRand2);
}

double Randomizer::GetFreePath(double aMeanFreePath, double& aRand)
{
	if(aRand<0 || aRand>1)
		aRand = Rand();
    if(aRand==0)
    	return 1e300;
    return aMeanFreePath * log(1./aRand);
}

double Randomizer::Rand()
{
	return gsl_rng_uniform_pos ((gsl_rng*)fRand);
}
double Randomizer::RandZero()
{
	return gsl_rng_uniform ((gsl_rng*)fRand);
}
double Randomizer::RandGauss(double sigma)
{
	return gsl_ran_gaussian ((gsl_rng*)fRand,sigma);
}

unsigned long int Randomizer::CreateIndependent()
{
	if(!fRand2)
	{
		fRand2 = gsl_rng_alloc (gsl_rng_ranlux389);//must be different type from fRand
		unsigned long int max2 = gsl_rng_max ((gsl_rng*)fRand2);
		unsigned long int max1 = gsl_rng_max ((gsl_rng*)fRand);
		fRand2Max = max1<max2?max1:max2;
		gsl_rng_set ((gsl_rng*)fRand2, gsl_rng_uniform_int((gsl_rng*)fRand, fRand2Max));
	}
	return gsl_rng_uniform_int((gsl_rng*)fRand2, fRand2Max);
}

}
