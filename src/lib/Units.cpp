/*
 * Units.cpp
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


#include "Units.h"
#include <math.h>

namespace mcray
{
const double Units::e = 0.0854245431348626; /* sqrt(alpha) */
const double Units::plank=6.582122020e-22; /* Plank constant (MeV*s)   */
const double Units::plank_ESU=1.05457266e-27;/* Plank constant (erg*s)   */
const double Units::phTemperature_mult=1.16e10;// K/MeV
const double Units::Mpc_in_cm=3.0856775807e24;     /*  (cm)  */
#ifdef ELMAG_TEST
const double Units::lightspeed=3e10; /* (cm/s) */
#else
const double Units::lightspeed=2.99792458e10; /* (cm/s) */
#endif
const double Units::MeVinESU=1.60217733e-6;  /* 1MeV in Ergs */
const double Units::YearInSec = 31557600;//number of seconds in 1 year

const Units units;//global units object

Units::Units()
{
	Eunit=1.; /* units of energy used (in MeV) */ //TODO: BUG!!! setting Eunit=1e-6 lead to different results: try for ex --test FERMI Absorption 1 2 8 1e15 0.06
	Eunit3=0; /* Eunit^3 */
	Lunit=0;   /* length unit in sm. */
	Tunit=0;   /* time unit in sec.  */
	Vunit=0;   /* Vunit=Lunit^3   */
	sigmaUnit=0;  /* sigmaUnit=Lunit^2   */
	SpecUnit=0;  /* SpecUnit=Eunit*Tunit*Vunit */
	outEunit=1e-6; /* Output data energy unit in MeV*/
	Bunit=0;  /* Magnetic field internal unit in Gauss (esu) */
	barn = 0;/*	1 barn (cross-section unit) in internal units */

	MeV = 1./Eunit;
	K = MeV/phTemperature_mult;
	eV = 1e-6*MeV;
	GeV = 1e3*MeV;
	TeV = 1e6*MeV;
	PeV = 1e9*MeV;
	Tunit=plank/Eunit;      /* time unit in sec.  */
	Lunit=lightspeed*Tunit; /* length unit in sm. */
	Eunit3 = Eunit*Eunit*Eunit;

	sigmaUnit=Lunit*Lunit;
	Vunit=Lunit*Lunit*Lunit;
	SpecUnit=Eunit*Tunit*Vunit;

	Bunit=MeVinESU*MeVinESU*pow(plank_ESU*lightspeed,-1.5);	// 1MeV^2 in Gauss
	Bunit*=(Eunit*Eunit);	//B unit in Gauss (esu)
	Gauss=1./Bunit;
	barn = 1e-24/sigmaUnit;
	Mpc = Mpc_in_cm/Lunit;
	kpc = 1e-3*Mpc;
	second = 1./Tunit;
	year = YearInSec * second;
	cm = 1./Lunit;
	cm3 = cm*cm*cm;
	Hz_photon = 2.*M_PI/second;
	erg = MeV/MeVinESU;
	J = (1e7*erg);
	W = (J/second);
}

} /* namespace mcray */
