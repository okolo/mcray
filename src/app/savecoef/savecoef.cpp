/*
 * CRbeam.cpp
 *
 * Authors:
 *       Oleg Kalashev
 *       Alexandr Korochkin
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

#include <cstdlib>
#include "ParticleStack.h"
#include "Utils.h"
#include "Cosmology.h"
#include "PropagationEngine.h"
#include <iostream>
#include <GZK.h>
#include <Inoue12IROSpectrum.h>
#include "PPP.h"
#include "ProtonPP.h"
#include "TableBackgrounds.h"
#include "ICS.h"
#include "Test.h"
#include "NeutronDecay.h"
#include "Deflection3D.h"
#include "SteckerEBL.h"
#include "CmdLine.h"
#include "Stecker16Background.h"
#include <stdlib.h>
#include "MathUtils.h"
#include "savecoef.h"

#ifndef _DEBUG
#include <sys/time.h>
#endif

using namespace std;
using namespace mcray;
using namespace Utils;
using namespace Backgrounds;
using namespace Interactions;


using namespace cors::cmdline;


/*
 * This is an example of user_main function performing user defined tasks
 * the function should be provided by end user
 */
int user_main(int argc, char** argv) {
    std::cout << "Number of available threads: " << omp_thread_count() << std::endl;
	CRbeam prog(argc, argv);
	return prog.run();
}

enum ClineParams
{
	PHelp = 1,PLog,PLogFilterThread,PLogLevel,PTrajectoryLog,
    PNparticles,PBatchSize,
    PRedshift,PFilamentZ,PEnergy,PMinEnergy,PPowerLaw,PMinSourceEnergy,PPrimary,PNoEMcascade,PSourceEvolM,
	PBackground,PBackgroundMult,PExtraBackground,PExtraBackgroundPhysical,PEBLCut,PMonoCmb,PFixedCmbT,PExtDeltaZconst,PExtPowerLow,PExtDeltaZexp,
	POutput,POverwriteOutput,PThinning,PRescalePPP,PEGMF,PEGMFLmin,PEGMFLmax,PEGMFNmodes,PRandomEGMF,PTurbulentEGMF,PDeflectionAccuracy,PTauPrint,POutputSuffix
};
#define xstr(s) str(s)
#define str(s) #s
#define M_POINT_SOURCE 1000

static CmdInfo commands[] = {
		{ PHelp,    CmdInfo::FLAG_NULL,     "-h",   "--help",       NULL,   "Show help information" },
		{ PLog,    CmdInfo::FLAG_NULL,     "-log",   "--log",       NULL,   "Enable logging" },
		{ PLogFilterThread,    CmdInfo::FLAG_ARGUMENT, "-logT",   "--log-thread",   "-1",    "Log thread filter (set to -1 to disable)" },
		{ PLogLevel,    CmdInfo::FLAG_ARGUMENT, "-logL",   "--log-level",   "1",    "Log level:\n"
																					"\tError = 0,\n"
																					"\tWarning = 1,\n"
																					"\tMessage = 2,\n"
																					"\tVerbose = 3"
		},
        { PTrajectoryLog,     CmdInfo::FLAG_NULL,     "-tlog",   "--tlog",       NULL,   "Log particle trajectories" },
        { PNparticles,    CmdInfo::FLAG_ARGUMENT, "-N",   "--nparticles",   "10000",     "Number of particles" },
        { PBatchSize,    CmdInfo::FLAG_ARGUMENT, "-bN",   "--batch",   "0",     "Number of particles in minibatch (0=auto)" },

		{ PRedshift,    CmdInfo::FLAG_ARGUMENT, "-z",   "--z",       NULL,   "Source z (or maximal z if --m_z par is given)" },
        { PFilamentZ,    CmdInfo::FLAG_ARGUMENT, "-zf",   "--z-filament",       "0.",   "Filament z (kill protons with z less than this value)" },
		{ PEnergy,    CmdInfo::FLAG_ARGUMENT, "-E",   "--emax",   "1e18",    "Initial (or maximal for power law) particle energy in eV" },
		{ PMinEnergy,    CmdInfo::FLAG_ARGUMENT, "-e",   "--emin",   "1e11",     "Minimal energy in eV" },
		{ PPowerLaw,    CmdInfo::FLAG_ARGUMENT, "-pl",   "--power",   "0",     "if>0 use E^(-power) injection " },
		{ PMinSourceEnergy,    CmdInfo::FLAG_ARGUMENT, "-es",   "--emin_source",   "0", "Minimal injection energy in eV (for power law injection only)" },

		{ PPrimary,    CmdInfo::FLAG_ARGUMENT, "-p",   "--primary",   "10",    "Primary particle: \n"
																					  "\tElectron = 0,\n"
																					  "\tPositron = 1,\n"
																					  "\tPhoton = 2,\n"
																					  "\tNeutron 9,\n"
																					  "\tProton 10" },
		{ PNoEMcascade,    CmdInfo::FLAG_NULL,     "-nEM",   "--noEMcascade",   NULL,   "don't calculate spectra of electrons and photons" },
		{ PSourceEvolM, CmdInfo::FLAG_ARGUMENT, "-m",   "--m_z",   xstr(M_POINT_SOURCE),    "if parameter is set (!=" xstr(M_POINT_SOURCE) "), use continuous source distribution ~(1+z)^m" },
		{ PBackground,    CmdInfo::FLAG_ARGUMENT, "-b",   "--background",   "3",    "EBL used:\n"
																							"0 - zero,\n"
																							"1 - Kneiske best fit (ELMAG),\n"
																							"2 - Kneiske minimal (ELMAG),\n"
																							"3 - Inoue 2012 Baseline,\n"
																							"4 - Inoue 2012 lower limit,\n"
																							"5 - Inoue 2012 upper limit,\n"
																							"6 - Franceschini 2008\n"
																							"7 - Kneiske best fit (digitized)\n"
																							"8 - Kneiske minimal (digitized)\n"
																							"9 - Stecker 2005\n"
																							"10 - Stecker 2016 lower limit,\n"
																							"11 - Stecker 2016 upper limit,\n"
                                                                                            "12 - Franceschini 2017,\n"
																							  },
        { PBackgroundMult,    CmdInfo::FLAG_ARGUMENT, "-bm",   "--backgr-mult",   "1",    "EBL multiplier" },
        { PExtraBackground, CmdInfo::FLAG_ARGUMENT, "-badd", "--add-backgr", "", "Additional matrix background file path (relative to ./tables)"},
        { PExtraBackgroundPhysical, CmdInfo::FLAG_NULL, "-badd_p",   "--add-backgr-phys",   NULL,   "Treat additional background data as physical (not comoving) density" },
		{ PEBLCut,    CmdInfo::FLAG_ARGUMENT,     "-bcut",   "--backgr-cut",   "100",   "Cut EBL above this energy [eV]" },
		{ PMonoCmb,    CmdInfo::FLAG_NULL,     "-mcmb",   "--monocmb",   NULL,   "Use monochromatic 6.3e-4 eV background instead of CMB" },
		{ PFixedCmbT,    CmdInfo::FLAG_NULL,     "-fcmb",   "--fixedcmb",   NULL,   "use CMB with fixed temperature (1+Zmax)2.73K" },

		{ PExtDeltaZconst,    CmdInfo::FLAG_ARGUMENT, "-bc",   "--eblDeltaZconst",   "-1",    "EBL extension DeltaZconst param (set to negative to disable EBL extension)"  },
		{ PExtPowerLow,    CmdInfo::FLAG_ARGUMENT, "-bp",   "--eblPowerLow",   "0",    "EBL extension PowerLow param"  },
		{ PExtDeltaZexp,    CmdInfo::FLAG_ARGUMENT, "-be",   "--eblDeltaZexp",   "1",    "EBL extension DeltaZexp param"  },

		{ POutput,    CmdInfo::FLAG_ARGUMENT, "-o",   "--output",   NULL,    "Output directory path (must not exist, will be created)" },
		{ POverwriteOutput,    CmdInfo::FLAG_NULL,     "-oo",   "--overwrite",   NULL,   "overwrite output dir if exists" },
		{ PThinning,    CmdInfo::FLAG_ARGUMENT, "-t",   "--thinning",   "0.9",    "Alpha thinning" },
		{ PRescalePPP,    CmdInfo::FLAG_ARGUMENT,   "-pt",   "--ppp-thinning",   "10",   "PPP Rescale Coefficient (double > 0, actual number of secondaries per interaction)" },

		{ PEGMF,    CmdInfo::FLAG_ARGUMENT, "-mf",   "--EGMF",   "1e-15",    "Extragalactic magnetic field in Gauss" },
		{ PEGMFLmin,    CmdInfo::FLAG_ARGUMENT, "-mflmin",   "--lminEGMF",   "0.05",    "Extragalactic magnetic field minimum scale in Mpc at z=0" },
		{ PEGMFLmax,    CmdInfo::FLAG_ARGUMENT, "-mflmax",   "--lmaxEGMF",   "5.0",    "Extragalactic magnetic field maximum scale in Mpc at z=0" },
		{ PEGMFNmodes,    CmdInfo::FLAG_ARGUMENT, "-mfnm",   "--nmodesEGMF",   "300",    "Number of modes in magnetic field spectrum" },
		{ PRandomEGMF,    CmdInfo::FLAG_NULL,     "-mfr",   "--randomEGMF",   NULL,   "randomize EGMF (use different EGMF configurations for different initial particles)" },
		{ PTurbulentEGMF,    CmdInfo::FLAG_NULL,     "-mft",   "--turbulentEGMF",   NULL,   "use turbulent EGMF with Kolmogorov spectrum" },
        { PDeflectionAccuracy,    CmdInfo::FLAG_ARGUMENT, "-mfac",   "--accEGMF",   "10",    "Deflection simulation accuracy factor" },
		{ PTauPrint,    CmdInfo::FLAG_NULL,     "-ptau",   "--print-tau",   NULL,   "print tau" },
        { POutputSuffix, CmdInfo::FLAG_ARGUMENT, "-os", "--output-suffix", "", "Output dir name suffix (appended to autogenerated names)"},
		{ 0 },
};

CRbeam::CRbeam(int argc, char** argv):
		fEmin_eV(0),
		fEmax_eV(0),
		fPowerLaw(0),
		fEminSource_eV(0),
		fZmax(0),
        fZfilament(0),
		fAlpha(0),
		fPPPrescaleCoef(10),
		fNoParticles(0),
        fBatchSize(0),
		fPrimary(Proton),
		fNoEM(false),
		fMz(M_POINT_SOURCE),
		fBackgroundModel(0),
        fEBLmult(1),
        fCustomBackground(""),
        fCustomBackgroundPhysical(false),
		fMonoCMB(false),
		fOutputDir(0),
		fOverwriteOutput(false),
		fFixedCmb(false),
		fEGMF(0.),
		fLminEGMF(0.),
		fLmaxEGMF(0.),
		fNmodesEGMF(0),
		fRandomizeEGMF(false),
		fTurbulentEGMF(false),
        fDeflectionAccuracy(10),
		fLogging(false),
        fTrajectoryLogging(false),
		fLogThread(-1),
		fLogLevel(WarningLL),
		fEBLMaxE(100.),
        fTauPrint(false),
		cmd(argc,argv,commands)
{
    std::cout << "sizeof(coord_type) = " << sizeof(coord_type) << std::endl;
	if( cmd.has_param("-h") )
	{
		cmd.printHelp(std::cerr);
		exit(255);//enable binary_check
	}
	fLogging = cmd(PLog);
    fTrajectoryLogging = cmd(PTrajectoryLog);
	fLogThread = cmd(PLogFilterThread);
	fLogLevel = (LogLevel)((int)cmd(PLogLevel));
	fEmin_eV = cmd(PMinEnergy);
	fEmax_eV = cmd(PEnergy);
	fPowerLaw = cmd(PPowerLaw);
	fEminSource_eV = cmd(PMinSourceEnergy);
	if(fEminSource_eV<fEmin_eV)
		fEminSource_eV=fEmin_eV;
	fZmax = cmd(PRedshift);
	fAlpha = cmd(PThinning);
	fPPPrescaleCoef = cmd(PRescalePPP);

	fEGMF = cmd(PEGMF);
	fLminEGMF = cmd(PEGMFLmin);
	fLmaxEGMF = cmd(PEGMFLmax);
	fNmodesEGMF = cmd(PEGMFNmodes);
	fRandomizeEGMF = cmd(PRandomEGMF);
	fTurbulentEGMF = cmd(PTurbulentEGMF);
	fDeflectionAccuracy = (int)cmd(PDeflectionAccuracy);
	//fMaxDeflection = cmd(PMaxDeflection);
	//fMaxDeflection *= (M_PI/180.);

	fNoParticles = cmd(PNparticles);
    fBatchSize = cmd(PBatchSize);
	int primary = cmd(PPrimary);
	fNoEM = cmd(PNoEMcascade);
	fMz = cmd(PSourceEvolM);
	fBackgroundModel = cmd(PBackground);
    fEBLmult = cmd(PBackgroundMult);
    fCustomBackground = cmd(PExtraBackground);
    fCustomBackgroundPhysical = cmd(PExtraBackgroundPhysical);
	fMonoCMB = cmd(PMonoCmb);
	fFixedCmb = cmd(PFixedCmbT);
	fOutputDir = cmd(POutput);
	fOverwriteOutput = cmd(POverwriteOutput);
	fExtDeltaZconst = cmd(PExtDeltaZconst);
	fExtPowerLow = cmd(PExtPowerLow);
	fExtDeltaZexp = cmd(PExtDeltaZexp);
	fEBLMaxE = cmd(PEBLCut);
    fTauPrint = cmd(PTauPrint);
    fZfilament = cmd(PFilamentZ);

	if(fOutputDir[0]=='\0')
		fOutputDir = 0;

	if(primary<Electron || primary>Proton)
		Exception::Throw("invalid primary particle type");
	fPrimary = (ParticleType)primary;
}

int CRbeam::run()
{
	int kAc = 1;
	fEBLMaxE *= units.eV;
	double Emax = fEmax_eV*units.eV;//eV
	if(!cosmology.IsInitialized())
		cosmology.Init(fZmax + 10);
	double zMin = 0.;
	double logStepK = pow(10,0.05/kAc);
	double Emin = fEmin_eV*units.eV;
	double EminSource = fEminSource_eV*units.eV;
	double alphaThinning = fAlpha;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	std::string outputDir;

	if(fOutputDir)
		outputDir = fOutputDir;
	else
	{
		outputDir = ".";
	}

	CompoundBackground backgr;

	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	if(fFixedCmb)
		cmbTemp *= (1.+fZmax);
	IBackground* b1 = new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp, 0., fZmax + 1., fFixedCmb);
	IBackground* b2 = 0;
	switch(fBackgroundModel)
	{
		case 0:
			b2 = 0;
			break;
		case 1:
			b2 = new ElmagKneiskeBestFit();
			break;
		case 2:
			b2 = new ElmagKneiskeMinimal();
			break;
		case 3:
			b2 = new Inoue12BaselineIROSpectrum();
			break;
		case 4:
			b2 = new Inoue12LowPop3IROSpectrum();
			break;
		case 5:
			b2 = new Inoue12UpperPop3IROSpectrum();
			break;
		case 6:
			b2 = new Franceschini08EBL();
			break;
		case 7:
			b2 = new Kneiske0309EBL();
			break;
		case 8:
			b2 = new Kneiske1001EBL();
			break;
		case 9:
			b2 = new Stecker2005EBL();
			break;
		case 10:
			b2 = new Stecker16LowerBackground();
			break;
		case 11:
			b2 = new Stecker16UpperBackground();
			break;
        case 12:
            b2 = new Franceschini17EBL();
            break;
		default:
			Exception::Throw("Invalid background model specified");
	}

	if(b2 && fEBLMaxE < b2->MaxE())
		b2 = new CuttedBackground(b2, 0, fEBLMaxE);

	if(b2 && fZmax>b2->MaxZ() && fExtDeltaZconst>=0)
	{
		b2 = new HighRedshiftBackgrExtension(b2, fExtDeltaZconst, fExtPowerLow, fExtDeltaZexp);
	}

	backgr.AddComponent(b1);
	if(b2)
	{
		backgr.AddComponent(b2, fEBLmult);
		double dens = BackgroundUtils::CalcIntegralDensity(*b2)*units.cm3;
		std::cerr << "EBL integral density [cm^{-3}]:" << dens << std::endl;
		std::ofstream backgrOut;
		backgrOut.open((outputDir + "/backgr").c_str(),std::ios::out);
		BackgroundUtils::Print(b2, 0., backgrOut, 100);
	}
    if(strlen(fCustomBackground)){
        MatrixBackground* bc = new MatrixBackground(fCustomBackground, fCustomBackground, !fCustomBackgroundPhysical, false);
        backgr.AddComponent(bc);
        double dens = BackgroundUtils::CalcIntegralDensity(*bc)*units.cm3;
        std::cerr << "Custom background integral density [cm^{-3}]: " << dens << std::endl;
        std::ofstream backgrOut((outputDir + "/cbackgr").c_str());
        BackgroundUtils::Print(bc, 0., backgrOut, 100);
    }

	double stepZ = fZmax<0.05 ? fZmax/2 : 0.025;
	double epsRel = 1e-3;
	double centralE1 = 6.3e-10/units.Eunit;//6.3e-4eV
	double n1 = 413*units.Vunit;//413cm^-3
	ConstFunction k1(centralE1);
	ConstFunction c1(n1);
	PBackgroundIntegral backgrI;
	if(fMonoCMB)
	{
		backgrI = new MonochromaticBackgroundIntegral(c1, k1, logStepK, epsRel);
	}
	else
	{
		backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, fZmax, epsRel);
	}

    Interactions::GammaPP gammaPP(backgrI);
    BinaryDataStorage ds(outputDir.c_str());
    double minFracE = 1e-5;
    double logStep = 0.01;
    gammaPP.SaveCoefficients(ds, Emin, Emax, minFracE, logStep, fZmax);
    std::cout << "output dir: " << outputDir << std::endl;

    // test of BinaryDataStorage
//    double x[] = {1.0,2,3,4,5};
//    std::vector<double> test(std::begin(x), std::end(x));
//    ds.save_array(test, "test");
//    std::vector<double> test2 = ds.load_array("test");
//    vector<double>::iterator iter = test2.begin();
//    for(iter; iter < test2.end(); iter++)
//    {
//        cout << *iter << " ";
//    }
//    cout << endl;
	return 0;
}

CRbeam::~CRbeam() {

}


