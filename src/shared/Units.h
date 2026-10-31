/*
 * Units.h
 *
 *  Created on: Sep 29, 2011
 *      Author: ok
 */

#ifndef UNITS_H_
#define UNITS_H_

namespace Utils
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
	double K;//Kelvin
	double Gauss;
private:
	static const double plank; /* Plank constant (MeV*s)   */
	static const double plank_ESU;/* Plank constant (erg*s)   */
};
	 extern const Units units;
} /* namespace MCPropag */
#endif /* UNITS_H_ */
