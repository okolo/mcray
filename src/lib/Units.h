/*
 * Units.h
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


#ifndef UNITS_H_
#define UNITS_H_

namespace mcray
{

class Units
{
public:
	Units();
	static const double Mpc_in_cm;     /*  (cm)  */
	static const double lightspeed; /* (cm/s) */
	static const double e; /* electric charge of electron in system plank=c=1*/

	double Eunit; /* units of energy used (in MeV) */
	double Eunit3; /* Eunit^3 */
	double Lunit;   /* length unit in cm. */
	double Tunit;   /* time unit in sec.  */
	double Vunit;   /* Vunit=Lunit^3   */
	double sigmaUnit;  /* sigmaUnit=Lunit^2   */
	double SpecUnit;  /* SpecUnit=Eunit*Tunit*Vunit */
	double outEunit; /* Output data energy unit in MeV*/
	double Bunit;  /* Magnetic field internal unit in Gauss (esu) */
	static const double MeVinESU;  /* 1MeV in Ergs */
	static const double YearInSec;//number of seconds in 1 year
	static const double phTemperature_mult;

//////   conventional units expressed in internal units   ////////////////

	double Mpc;
	double kpc;
	double barn;
	double MeV;
	double GeV;
	double eV;
	double TeV;
	double PeV;
	double second;
	double year;
	double cm;
	double cm3;//cubic cm
	double erg;
	double K;//Kelvin
	double Gauss;

	double Hz_photon;// energy of 1 Hz EM wave photon (in internal units)
	double J; //Joule
	double W; //Watt
private:
	static const double plank; /* Plank constant (MeV*s)   */
	static const double plank_ESU;/* Plank constant (erg*s)   */
};
	 extern const Units units;
} /* namespace mcray */
#endif /* UNITS_H_ */
