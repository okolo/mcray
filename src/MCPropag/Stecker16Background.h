//
// Created by ok on 01/09/2016.
//

#ifndef PROPAGATION_STECKER16BACKGROUND_H
#define PROPAGATION_STECKER16BACKGROUND_H

#include "Background.h"
//#include "Vector.h"
//#include "TableFunc.h"
#include <string>
#include <vector>
//#include "DataReader.h"
using namespace std;

#include "Background.h"

namespace Backgrounds {

    using namespace MCPropag;

    class Stecker16Background : public IBackground {
    public:
        Stecker16Background(string aTableFile, std::string aName) :
                fMatrixFunction(aTableFile),
                fName(aName)
        {}

        static int UnitTest();
/*!
	 Calculate density E*dn(E)/dE at given redshift in internal units. Physical density should be returned, not a comoving one
	 @param[in]  aE  background photon energy (in internal units)
	 @param[in]  aZ  red shift
	 @return background photon density E*dn(E)/dE in internal units
*/
        virtual double n(double aE, double aZ) const;
        virtual double MinE() const { return fMatrixFunction.MinArg(2)*units.eV; }
        virtual double MaxE() const { return fMatrixFunction.MaxArg(2)*units.eV; }
        virtual double MinZ() const { return fMatrixFunction.MinArg(1); }
        virtual double MaxZ() const { return fMatrixFunction.MaxArg(1); }
        virtual std::string Name() const { return fName; };

    private:
        MatrixFunction fMatrixFunction;
        std::string fName;
    };

    class Stecker16LowerBackground : public Stecker16Background {
    public:
        Stecker16LowerBackground();
    };

    class Stecker16UpperBackground : public Stecker16Background {
    public:
        Stecker16UpperBackground();
    };
}
#endif //PROPAGATION_STECKER16BACKGROUND_H
