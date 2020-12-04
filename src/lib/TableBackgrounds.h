/*
 * TableBackgrounds.h
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


#ifndef TABLEBACKGROUNDS_H_
#define TABLEBACKGROUNDS_H_

#include "Background.h"

namespace Backgrounds {

using namespace mcray;

class Kneiske1001EBL : public TableBackground{
public:
	Kneiske1001EBL();
	virtual ~Kneiske1001EBL();
protected:
	//convert energy in internal units to table x scale
	double EtoRecordX(double aE, double aZ) const;
	//retrieve energy in internal units from the data record
	double recordToE(double aX, double aZ) const;
	////retrieve concentration E*dn/dE in internal units from the data record
	double recordToN(double aX, double aY, double aZ) const;
};

class Kneiske0309EBL : public TableBackground{
public:
	Kneiske0309EBL();
	virtual ~Kneiske0309EBL();
protected:
	//convert energy in internal units to table x scale
	double EtoRecordX(double aE, double aZ) const;
	//retrieve energy in internal units from the data record
	double recordToE(double aX, double aZ) const;
	////retrieve concentration E*dn/dE in internal units from the data record
	double recordToN(double aX, double aY, double aZ) const;
};

/// EBL from Table1,2 of http://arxiv.org/pdf/0805.1841v2.pdf
class Franceschini08EBL : public TableBackground{
public:
	Franceschini08EBL();
	virtual ~Franceschini08EBL();
protected:
	//convert energy in internal units to table x scale
	double EtoRecordX(double aE, double aZ) const;
	//retrieve energy in internal units from the data record
	double recordToE(double aX, double aZ) const;
	////retrieve concentration E*dn/dE in internal units from the data record
	double recordToN(double aX, double aY, double aZ) const;
};

class ElmagKneiskeBestFit : public MatrixBackground
{
public:
	/// Kneiske BestFit Background used in Elmag
	/// There is an error in Elmag 1.02 and earlier, the tables are treated as physical concentrations
	/// But in fact they are built for comoving frame (see fig 1 of http://lanl.arxiv.org/abs/astro-ph/0309141v1)
	/// The parameter here is left just for comparison to Elmag 1.02
	/// In Elmag 2.0 this bug was fixed
	ElmagKneiskeBestFit(bool aIsComoving = true);
};

class ElmagKneiskeMinimal : public MatrixBackground
{
public:
	/// Kneiske Lower limit Background used in Elmag
	/// There is an error in Elmag 1.02 and earlier, the tables are treated as physical concentrations
	/// But in fact they are built for comoving frame
	/// The parameter here is left just for comparison to Elmag 1.02
	/// In Elmag 2.0 this bug was fixed
	ElmagKneiskeMinimal(bool aIsComoving = true);
};

class Franceschini17EBL : public TableBackground{
public:
    Franceschini17EBL();
    virtual ~Franceschini17EBL();
protected:
    //convert energy in internal units to table x scale
    double EtoRecordX(double aE, double aZ) const;
    //retrieve energy in internal units from the data record
    double recordToE(double aX, double aZ) const;
    ////retrieve concentration E*dn/dE in internal units from the data record
    double recordToN(double aX, double aY, double aZ) const;
};

} /* namespace Backgrounds */

#endif /* TABLEBACKGROUNDS_H_ */
