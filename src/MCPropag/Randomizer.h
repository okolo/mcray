/*
 * File:   Randomizer.h
 * Author: ok
 *
 * Created on August 26, 2011, 7:58 PM
 */

#ifndef RANDOMIZER_H
#define	RANDOMIZER_H

namespace MCPropag
{
///not expected to be thread-safe
class Randomizer {
public:
    Randomizer(unsigned long int aSeed=0);
    virtual ~Randomizer();
    double Rand();//generate random number in the range [0,1] <- EXCLUDE ZERO AND ONE
    double RandZero();//generate random number in the range [0,1) <- INCLUDE ZERO
    double RandGauss(double sigma);

	/*!
	@param[in]  aMeanFreePath mean free path of a particle
	@param[in,out]  aRand random number used for sampling. If 0<aRand<1 than input value is used for sampling,
	otherwise randomly generated number is used and returned as output parameter
	@return sampled value of the free path
	*/
    double GetFreePath(double aMeanFreePath, double& aRand);

    //generate seed for independent Randomizer
    unsigned long int CreateIndependent();
private:
    void* fRand;
    void* fRand2;//used in CreateIndependent() only
    unsigned long int fRand2Max;
};
}
#endif	/* RANDOMIZER_H */
