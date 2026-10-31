/*
 * Inoue12IROSpectrum.cpp
 *
 *  Created on: Feb 24, 2013
 *      Author: ok
 */

#include "Inoue12IROSpectrum.h"

namespace Backgrounds
{
using namespace MCPropag;

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
