/*
 * Stecker16Background.cpp
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


#include "Stecker16Background.h"
#include "Units.h"

namespace Backgrounds {

    Stecker16LowerBackground::Stecker16LowerBackground() :
            Stecker16Background(TABLES_DIR  "stecker_ebl16/comoving_enerdens_lo.csv", "Stecker_16_lower") {
}

Stecker16UpperBackground::Stecker16UpperBackground() :
        Stecker16Background(TABLES_DIR "stecker_ebl16/comoving_enerdens_up.csv", "Stecker_16_upper")
{
}


double Stecker16Background::n(double E, double z) const {
    // matrix function returns the spectrum E*dn/dE in sm^-3 [E]=eV in comoving volume
    double result = fMatrixFunction.f(z, E/units.eV) / units.Hz_photon * units.erg;
    double dl = (1.+z);
    //convert to physical density in internal units
    return dl*dl*dl/units.cm3*result;
}

int Stecker16Background::UnitTest()
{
    Stecker16UpperBackground backgrUpper;
    BackgroundUtils::UnitTest(backgrUpper);
    Stecker16LowerBackground backgrLower;
    BackgroundUtils::UnitTest(backgrLower);
    return 0; // Unwarning stub. Why was there no return statement?
}

}
