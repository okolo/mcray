/*
 * Inoue12IROSpectrum.cpp
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

#include "Inoue12IROSpectrum.h"

namespace Backgrounds
{
using namespace mcray;

Inoue12IROSpectrum::Inoue12IROSpectrum(std::string aName, std::string aDataFile):
MatrixBackground(aName, aDataFile, false, true)
{
}

double Inoue12IROSpectrum::n(double aE, double z) const
{
	double E_eV = aE/units.eV;
	return MatrixBackground::n(aE,z)/E_eV;
}

Inoue12BaselineIROSpectrum::Inoue12BaselineIROSpectrum():
Inoue12IROSpectrum("Inoue12Baseline", "EBL_Inoue/EBL_proper_baseline.dat")
{
}

Inoue12LowPop3IROSpectrum::Inoue12LowPop3IROSpectrum():
Inoue12IROSpectrum("Inoue12low", "EBL_Inoue/EBL_proper_low_pop3.dat")
{
}

Inoue12UpperPop3IROSpectrum::Inoue12UpperPop3IROSpectrum():
Inoue12IROSpectrum("Inoue12hi", "EBL_Inoue/EBL_proper_up_pop3.dat")
{
}

}
