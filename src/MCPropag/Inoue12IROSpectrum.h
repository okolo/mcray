/*
 * Inoue12IROSpectrum.h
 *
 *  Created on: Feb 24, 2013
 *      Author: ok
 */

#ifndef INOUE12IROSPECTRUM_H_
#define INOUE12IROSPECTRUM_H_

#include "Background.h"

namespace Backgrounds {

using namespace MCPropag;

class Inoue12IROSpectrum  : public MatrixBackground{
public:
	Inoue12IROSpectrum(std::string aName, std::string aDataFile);
	virtual double n(double E, double z) const;
};

class Inoue12BaselineIROSpectrum  : public Inoue12IROSpectrum{
public:
	Inoue12BaselineIROSpectrum();
};

class Inoue12LowPop3IROSpectrum  : public Inoue12IROSpectrum{
public:
	Inoue12LowPop3IROSpectrum();
};

class Inoue12UpperPop3IROSpectrum  : public Inoue12IROSpectrum{
public:
	Inoue12UpperPop3IROSpectrum();
};

}

#endif /* INOUE12IROSPECTRUM_H_ */
