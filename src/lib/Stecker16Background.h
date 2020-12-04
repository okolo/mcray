/*
 * Stecker16Background.h
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

    using namespace mcray;

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
