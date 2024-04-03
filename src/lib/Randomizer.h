/*
 * Randomizer.h
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


#ifndef RANDOMIZER_H
#define	RANDOMIZER_H

namespace mcray
{
///not expected to be thread-safe
class Randomizer {
public:
    Randomizer(unsigned long int aSeed=0);
    virtual ~Randomizer();
    double Rand();//generate random number in the range [0,1]
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
