/*
 * Deflection1D.h
 *
 *  Created on: Mar 22, 2012
 *      Author: ok
 */

#ifndef DEFLECTION1D_H_
#define DEFLECTION1D_H_

#include "Interaction.h"
#include "MathUtils.h"

namespace Interactions {
using namespace MCPropag;
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

} /* namespace MCPropag */
#endif /* DEFLECTION1D_H_ */
