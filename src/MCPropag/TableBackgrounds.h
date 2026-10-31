/*
 * TableBackgrounds.h
 *
 *  Created on: Oct 20, 2011
 *      Author: ok
 */

#ifndef TABLEBACKGROUNDS_H_
#define TABLEBACKGROUNDS_H_

#include "Background.h"

namespace Backgrounds {

using namespace MCPropag;

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

// Gilmore 2012 fixed
class Gilmore12fixedEBL : public MatrixBackground
{
public:
	Gilmore12fixedEBL(bool aIsComoving = false);
};

// Gilmore 2012 fiducial
class Gilmore12fiducialEBL : public MatrixBackground
{
public:
	Gilmore12fiducialEBL(bool aIsComoving = false);
};

// Gilmore 2012 fixed with CIBER minimal bump
class Gilmore12fixedCiberEBL : public MatrixBackground
{
public:
	Gilmore12fixedCiberEBL(bool aIsComoving = false);
};

// Gilmore 2012 fiducial with CIBER minimal bump
class Gilmore12fiducialCiberEBL : public MatrixBackground
{
public:
	Gilmore12fiducialCiberEBL(bool aIsComoving = false);
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
