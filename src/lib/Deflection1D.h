/*
 * Deflection1D.h
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

#ifndef DEFLECTION1D_H_
#define DEFLECTION1D_H_

#include "Interaction.h"
#include "MathUtils.h"

namespace Interactions {
using namespace mcray;
using namespace Utils;

//Small deflection approximation (assuming spherical symmetry)
class Deflection1D : public DeflectionInteraction{
    class RandomDirectionB : public Function
    {
    public:
    	RandomDirectionB(double aAbsBvalue, unsigned long int aSeed=0) : fRandomizer(aSeed), fAbsBvalue(fabs(aAbsBvalue)){}
    	double f(double aZ) const;//sign of B is disregarded since sum of squares is used
    	virtual Function* Clone() const;
    private:
    	mutable Randomizer fRandomizer;
    	double fAbsBvalue;
    };
    class ConstComovingCorrelationLength : public Function
    {
    public:
    	ConstComovingCorrelationLength(double aCorL) : fCorL(aCorL){}
    	double f(double aZ) const;
    	virtual Function* Clone() const;
    private:
    	double fCorL;;
    };
public:
	//Constructor with dependences of magnetic field abs value and correlation length on distance to the source in terms of redshift
	Deflection1D(Function* aB, Function* aLcor);

	Deflection1D(double aBgauss, double aLcorMpc, unsigned long int aRandomSeed)
	{
		fLcor = new ConstComovingCorrelationLength(aLcorMpc*units.Mpc);
		fB = new RandomDirectionB(aBgauss, aRandomSeed);
	}
	virtual ~Deflection1D();

	virtual bool Propagate(cosmo_time deltaT, Particle &aParticle, Randomizer &aRandomizer) const;
    virtual double Rate(const Particle& aParticle) const;
    virtual DeflectionInteraction* Clone() const
    {
    	return new Deflection1D(fB->Clone(), fLcor->Clone());
    }
//protected:
	//Used for test purposes only (constant B projection value is not physical)
    //This method can be used for comparison with kinetic equation based code
    //with deflection taken into account in the form of effective electron decay
	Deflection1D(double aBgauss, double aLcorMpc)
	{
		fLcor = new ConstComovingCorrelationLength(aLcorMpc*units.Mpc);
		fB = new ConstFunction(aBgauss);
	}
private:
    // fB->f(z) returns strength of (transverse) extragalactic B-field [Gauss] at redshift z at the distance from the observer
    // corresponding to light travel time from redshift z to 0
    // (spherical symmetry is assumed)
    SafePtr<Function> fB;//sign of B is disregarded since sum of squares is used
    // fLcor->f(z) returns B correlation length at redshift z at the distance from the
    // observer corresponding to light travel time from redshift z to 0
    // (spherical symmetry is assumed)
    SafePtr<Function> fLcor;//sign of B is disregarded since sum of squares is used
};

} /* namespace mcray */
#endif /* DEFLECTION1D_H_ */
