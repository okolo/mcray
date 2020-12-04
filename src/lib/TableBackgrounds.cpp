/*
 * TableBackgrounds.cpp
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

#include "TableBackgrounds.h"
#include "Units.h"
#include <math.h>

namespace Backgrounds {


Kneiske1001EBL::Kneiske1001EBL():
		TableBackground("Kneiske1001 lower-limit")
{
	const char* files[] = {"0","0.1","0.3","0.8","2",0};
	Init("Kneiske1001", files, true);
}

Kneiske1001EBL::~Kneiske1001EBL() {
}

//convert energy in internal units to table x scale
double Kneiske1001EBL::EtoRecordX(double aE, double aZ) const
{
	double lambda = 2.*M_PI/aE;//wavelength in internal units
	lambda *= units.Lunit*1e4;//wavelength in microns
	return log10(lambda);
}

//retrieve energy in internal units from the data record
double Kneiske1001EBL::recordToE(double aX, double aZ) const
{
	double lambda = pow(10.,aX);//wavelength in microns
	lambda /= (units.Lunit*1e4);//wavelength in internal units
	return 2.*M_PI/lambda;
}

////retrieve concentration E*dn/dE in internal units from the data record
double Kneiske1001EBL::recordToN(double aX, double aY, double aZ) const
{//[y] = nW m^-2 sr^-1 in comoving frame
	double comovingFactor = (1.+aZ);
	comovingFactor *= (comovingFactor*comovingFactor);
	double aE = recordToE(aX, aZ);
	double y = pow(10.,aY);
	y *= (4.*M_PI*1e-9/(units.MeVinESU*1e-7/units.Eunit)*1e-4*units.Lunit*units.Lunit*units.Tunit);
	double result = y/aE*comovingFactor;
	return result;
}

Kneiske0309EBL::Kneiske0309EBL():
		TableBackground("Kneiske0309 best fit")
{
	const char* files[] = {"0","0.2","0.4","0.6","1","2","3","4",0};
	Init("Kneiske0309", files, true);
}

Kneiske0309EBL::~Kneiske0309EBL() {
}

//convert energy in internal units to table x scale
double Kneiske0309EBL::EtoRecordX(double aE, double aZ) const
{
	double lambda = 2.*M_PI/aE;//wavelength in internal units
	lambda *= units.Lunit*1e4;//wavelength in microns
	return log10(lambda);
}

//retrieve energy in internal units from the data record
double Kneiske0309EBL::recordToE(double aX, double aZ) const
{
	double lambda = pow(10.,aX);//wavelength in microns
	lambda /= (units.Lunit*1e4);//wavelength in internal units
	return 2.*M_PI/lambda;
}

////retrieve concentration E*dn/dE in internal units from the data record
double Kneiske0309EBL::recordToN(double aX, double aY, double aZ) const
{//[y] = nW m^-2 sr^-1 in comoving frame
	double comovingFactor = (1.+aZ);
	comovingFactor *= (comovingFactor*comovingFactor);
	double aE = recordToE(aX, aZ);
	double y = pow(10.,aY);
	y *= (4.*M_PI*1e-9/(units.MeVinESU*1e-7/units.Eunit)*1e-4*units.Lunit*units.Lunit*units.Tunit);
	double result = y/aE*comovingFactor;
	return result;
}

//////////

Franceschini08EBL::Franceschini08EBL():
		TableBackground("Franceschini 08")
{
	const char* files[] = {"0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6","1.8","2.0", 0};
	Init("Franceschini08EBL", files, true);
}

Franceschini08EBL::~Franceschini08EBL() {
}

//convert energy in internal units to table x scale
double Franceschini08EBL::EtoRecordX(double aE, double aZ) const
{
	return log10(aE/units.eV);
}

//retrieve energy in internal units from the data record
double Franceschini08EBL::recordToE(double aX, double aZ) const
{
	return pow(10.,aX)*units.eV;
}

////retrieve concentration E*dn/dE in internal units from the data record
double Franceschini08EBL::recordToN(double aX, double aY, double aZ) const
{//[y] = log(E*dn/dE) in cm^-3 in physical frame
	double y = pow(10.,aY);
	return y/units.cm3;
}

////////////

ElmagKneiskeBestFit::ElmagKneiskeBestFit(bool aIsComoving):
		MatrixBackground("Elmag best fit", "kneiskeElmagBestFit.dat", aIsComoving, true)
{
}

ElmagKneiskeMinimal::ElmagKneiskeMinimal(bool aIsComoving):
		MatrixBackground("Elmag lower-limit", "kneiskeElmagLowerLimit.dat", aIsComoving, true)
{
}

///////////////
Franceschini17EBL::Franceschini17EBL():
        TableBackground("Franceschini17")
{
    const char* files[] = {"0", "0.2", "0.4", "0.6", "0.8", "1", "1.2", "1.4", "1.8", "2", 0};
    Init("Franceschini17", files, true);
}

Franceschini17EBL::~Franceschini17EBL() {
}

//convert energy in internal units to table x scale
double Franceschini17EBL::EtoRecordX(double aE, double aZ) const
{
    return log10(aE/units.eV);
}

//retrieve energy in internal units from the data record
double Franceschini17EBL::recordToE(double aX, double aZ) const
{
    return pow(10.,aX)*units.eV;
}

////retrieve concentration E*dn/dE in internal units from the data record
double Franceschini17EBL::recordToN(double aX, double aY, double aZ) const
{
    return pow(10., aY)/units.cm3;
}

} /* namespace Backgrounds */

