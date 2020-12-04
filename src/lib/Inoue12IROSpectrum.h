/*
 * Inoue12IROSpectrum.h
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

#ifndef INOUE12IROSPECTRUM_H_
#define INOUE12IROSPECTRUM_H_

#include "Background.h"

namespace Backgrounds {

using namespace mcray;

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
