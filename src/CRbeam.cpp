
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
#include "CRbeam.h"
#include "Test.h"
#include "NeutronDecay.h"
#include "Deflection3D.h"
#include "SteckerEBL.h"
#include "CmdLine.h"
#include "Stecker16Background.h"
#include <omp.h>
#include <stdlib.h>
#include "MathUtils.h"
#include "Randomizer.h"
#include "Jet.h"
#include <math.h>

#ifndef _DEBUG
#include <sys/time.h>
#endif

//extern "C" {
//#include <fftw3.h>
//}

using namespace std;
using namespace MCPropag;
using namespace Utils;
using namespace Backgrounds;
using namespace Interactions;


using namespace cors::cmdline;


/*
 * This is an example of user_main function performing user defined tasks
 * the function should be provided by end user
 */
int user_main(int argc, char** argv) {
	CRbeam prog(argc, argv);
	return prog.run();
}

enum ClineParams
{
	PHelp = 1,PLog,PLogFilterThread,PLogLevel,PTrajectoryLog,
    PNparticles,PBatchSize,PNumThreads,
    PRedshift,PZmin,PEnergy,PMinEnergy,PPowerLaw,PMinSourceEnergy,PPrimary,PNoEMcascade,PSourceEvolM,
	PBackground,PBackgroundMult,PExtraBackground,PExtraBackgroundPhysical,PEBLCut,PMonoCmb,PFixedCmbT,PNoCmb,PExtDeltaZconst,PExtPowerLow,PExtDeltaZexp,
	PMFFile,
	POutput,POverwriteOutput,PSaveSourceDir,PSaveMagnetic,PThinning,PEThinning,PRescalePPP,PEGMF,PEGMFL,PRandomEGMF,PTurbulentEGMF,PTurbulentFFT, PTurbulentEGMFstep,
	PTurbulentEGMFMinVarScale, PTauPrint,
	PJetAngle,PJetTheta,PJetPhi,PJetType,PJetParam,
	PObserver, PObserverX,PObserverY,PObserverZ,PObserverRadius,
	PRandSeed,PTraceParticle,PAccuracy,PMFCut,
	PSourceX,PSourceY,PSourceZ,PInteractions
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
        { PTrajectoryLog,     CmdInfo::FLAG_ARGUMENT,     "-tlog",   "--tlog",       NULL,   "Log particle trajectories" },

        { PNparticles,    CmdInfo::FLAG_ARGUMENT, "-N",   "--nparticles",   "10000",     "Number of particles" },
        { PBatchSize,    CmdInfo::FLAG_ARGUMENT, "-bN",   "--batch",   "0",     "Number of particles in minibatch (0=auto)" },
		{ PNumThreads,     CmdInfo::FLAG_ARGUMENT,     "-nth",   "--nthreads",       "0",   "Number of threads; 0=Max available" },

		{ PRedshift,    CmdInfo::FLAG_ARGUMENT, "-z",   "--z",       NULL,   "Source z (or maximal z if --m_z par is given)" },
		{ PZmin,    CmdInfo::FLAG_ARGUMENT, "-zmin",   "--zmin",       NULL,   "Observer z" },
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
																							"7 - Kneiske best fit 2004 (digitized)\n"
																							"8 - Kneiske minimal 2010 (digitized)\n"
																							"9 - Stecker 2005\n"
																							"10 - Stecker 2016 lower limit,\n"
																							"11 - Stecker 2016 upper limit,\n"
                                                                                            "12 - Franceschini 2017,\n"
                                                                                            "13 - Gilmore 2012 fixed,\n"
                                                                                            "14 - Gilmore 2012 fiducial,\n"
                                                                                            "15 - Gilmore 2012 fixed + CIBER min bump,\n"
                                                                                            "16 - Gilmore 2012 fiducial + CIBER min bump,\n"
																							  },
        { PBackgroundMult,    CmdInfo::FLAG_ARGUMENT, "-bm",   "--backgr-mult",   "1",    "EBL multiplier" },
        { PExtraBackground, CmdInfo::FLAG_ARGUMENT, "-badd", "--add-backgr", "", "Additional matrix background file path (relative to ./tables)"},
        { PExtraBackgroundPhysical, CmdInfo::FLAG_NULL, "-badd_p",   "--add-backgr-phys",   NULL,   "Treat additional background data as physical (not comoving) density" },
		{ PEBLCut,    CmdInfo::FLAG_ARGUMENT,     "-bcut",   "--backgr-cut",   "100",   "Cut EBL above this energy [eV]" },
		{ PMonoCmb,    CmdInfo::FLAG_NULL,     "-mcmb",   "--monocmb",   NULL,   "Use monochromatic 6.3e-4 eV background instead of CMB" },
		{ PFixedCmbT,    CmdInfo::FLAG_NULL,     "-fcmb",   "--fixedcmb",   NULL,   "use CMB with fixed temperature (1+Zmax)2.73K" },
		{ PNoCmb,    CmdInfo::FLAG_NULL,     "-ncmb",   "--nocmb",   NULL,   "zero CMB" },

		{ PExtDeltaZconst,    CmdInfo::FLAG_ARGUMENT, "-bc",   "--eblDeltaZconst",   "-1",    "EBL extension DeltaZconst param (set to negative to disable EBL extension)"  },
		{ PExtPowerLow,    CmdInfo::FLAG_ARGUMENT, "-bp",   "--eblPowerLow",   "0",    "EBL extension PowerLow param"  },
		{ PExtDeltaZexp,    CmdInfo::FLAG_ARGUMENT, "-be",   "--eblDeltaZexp",   "1",    "EBL extension DeltaZexp param"  },

		{ PMFFile,    CmdInfo::FLAG_ARGUMENT, "-mffile",   "--mf-file",   NULL,    "Input magnetic field file. Cancel turbulent MF configuration." },
		{ POutput,    CmdInfo::FLAG_ARGUMENT, "-o",   "--output",   NULL,    "Output directory path (must not exist, will be created)" },
		{ POverwriteOutput,    CmdInfo::FLAG_NULL,     "-oo",   "--overwrite",   NULL,   "overwrite output dir if exists" },
		{ PSaveSourceDir,    CmdInfo::FLAG_NULL,     "-os",   "--outsource",   NULL,   "write source direction in output file" },
		{ PSaveMagnetic,    CmdInfo::FLAG_NULL,     "-omf",   "--outmagnetic",   NULL,   "write magnetic field to file" },
		{ PThinning,    CmdInfo::FLAG_ARGUMENT, "-t",   "--thinning",   "0",    "Alpha thinning" },
		{ PEThinning,    CmdInfo::FLAG_ARGUMENT, "-te",   "--thinningEnergy",   "0",    "Thinning occurs below this energy eV" },
		{ PRescalePPP,    CmdInfo::FLAG_ARGUMENT,   "-pt",   "--ppp-thinning",   "10",   "PPP Rescale Coefficient (double > 0, actual number of secondaries per interaction)" },

		{ PEGMF,    CmdInfo::FLAG_ARGUMENT, "-mf",   "--EGMF",   "1e-15",    "Extragalactic magnetic field in Gauss" },
		{ PEGMFL,    CmdInfo::FLAG_ARGUMENT, "-mfl",   "--lEGMF",   "1",    "Extragalactic magnetic field correlation length in Mpc at z=0" },
		{ PRandomEGMF,    CmdInfo::FLAG_NULL,     "-mfr",   "--randomEGMF",   NULL,   "randomize EGMF (use different EGMF configurations for different initial particles)" },
		{ PTurbulentEGMF,    CmdInfo::FLAG_NULL,     "-mft",   "--turbulentEGMF",   NULL,   "use turbulent EGMF with Kolmogorov spectrum" },
		{ PTurbulentFFT,    CmdInfo::FLAG_NULL,     "-mfft",   "--turbulentFFT",   NULL,   "use FFTw library for EGMF generation" },
		{ PTurbulentEGMFstep, CmdInfo::FLAG_ARGUMENT, "-mfstep", "--EGMFstep", "1.01", "Ratio of two neighbourhood modes in EGMF spectrum"},
		{ PTurbulentEGMFMinVarScale, CmdInfo::FLAG_ARGUMENT, "-mfvarsc", "--EGMFMinVarScale", "0.1", "Minimal variability scale of the MF as fraction of Lc"},
        { PTauPrint,    CmdInfo::FLAG_NULL,     "-ptau",   "--print-tau",   NULL,   "print tau" },
        { PJetAngle,    CmdInfo::FLAG_ARGUMENT,     "-jet",   "--jet-angle",   "0",   "jet scan angle of cone from 0 (|| Oz) to Pi (isotropic emission). Note: this param uses as cut off for non-cone jet types" },
        { PJetTheta,    CmdInfo::FLAG_ARGUMENT,     "-jett",   "--jet-theta",   "0",   "jet angle between Oz and mean jet direction" },
        { PJetPhi,    CmdInfo::FLAG_ARGUMENT,     "-jetp",   "--jet-phi",   "0",   "jet angle between Ox and mean jet direction" },
        { PJetType,    CmdInfo::FLAG_ARGUMENT,     "-jettp",   "--jet-type",   "0",   "jet type:\n"
        																						"0 - cone type\n"
        																						"1 - Gauss type\n"
        																						"2 - von Miser-Fisher (vMF) type\n"
        																						},
        { PJetParam,    CmdInfo::FLAG_ARGUMENT,     "-jetprm",   "--jet-param",   "0",   "jet distribution param: 1-sigma for Gauss or kappa for vMF" },
        { PObserver,    CmdInfo::FLAG_ARGUMENT,     "-obs",   "--observer-type",   "-1",   "observer type: -1 - default built-in observer, 1 - sphere observer with radius in (0,0,0) default coords, 2 - inversed sphere observer for destroy particles after movement not travel time" },
        { PObserverX,    CmdInfo::FLAG_ARGUMENT,     "-obsx",   "--observer-x",   "0",   "observer X coord, kpc" },
        { PObserverY,    CmdInfo::FLAG_ARGUMENT,     "-obsy",   "--observer-y",   "0",   "observer Y coord, kpc" },
        { PObserverZ,    CmdInfo::FLAG_ARGUMENT,     "-obsz",   "--observer-z",   "0",   "observer Z coord, kpc" },
        { PObserverRadius,    CmdInfo::FLAG_ARGUMENT,     "-obsr",   "--observer-radius",   "100",   "observer radius, kpc" },
        { PRandSeed,    CmdInfo::FLAG_ARGUMENT,     "-rseed",   "--random-seed",   "-1",   "initial seed for randomizer, generates by current time if rseed<0 (default)" },
		{ PTraceParticle,    CmdInfo::FLAG_NULL,     "-otrace",   "--out-trace-particles",   NULL,   "write trajectory of all particles, check free disk space before using" },
		{ PAccuracy,    CmdInfo::FLAG_ARGUMENT,     "-acc",   "--accuracy",   "10",   "propagation accuracy" },
		{ PMFCut,    CmdInfo::FLAG_ARGUMENT,     "-mfcut",   "--mf-cut",   "-1",   "reset magnetic field to 0 after *this param* corr lengths, debug feature, not need for physics" },
        { PSourceX,    CmdInfo::FLAG_ARGUMENT,     "-srcx",   "--source-x",   "0",   "source X coord, check redshift, kpc" },
        { PSourceY,    CmdInfo::FLAG_ARGUMENT,     "-srcy",   "--source-y",   "0",   "source Y coord, check redshift, kpc" },
        { PSourceZ,    CmdInfo::FLAG_ARGUMENT,     "-srcz",   "--source-z",   "0",   "source Z coord, check redshift, kpc"},
		{ PInteractions,    CmdInfo::FLAG_NULL,     "-inter",   "--enable-interactions",   NULL,   "add interactions GZK, Neutron decay and ProtonPPcel" },
 		{ 0 },
};

/*void test_fftw(){
	const int N = 8;
	fftw_complex *in, *out;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

	fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	for (int i = 0; i < N; ++i) {
		in[i][0] = sin(2.0 * M_PI * i / N); // Вещественная часть
		in[i][1] = 0.0;                     // Мнимая часть
	}

	fftw_execute(plan);

	std::cout << "FFTW TESTING:\nАмплитуды:" << std::endl;
	for (int i = 0; i < N/2 + 1; ++i) {
		double amplitude = sqrt(out[i][0]*out[i][0] + out[i][1]*out[i][1]) / N * 2;
		std::cout << "Частота " << i << ": " << amplitude << std::endl;
	}

	fftw_destroy_plan(plan);
	fftw_free(in);
	fftw_free(out);

	fftw_cleanup();
}//*/

class CRbeamLogger : public Logger{
    int fB_slot;
    MagneticField* fB;
public:
    CRbeamLogger(std::ostream& aOut, int aB_slot=-1, MagneticField* aB=0):
            Logger(aOut),
            fB_slot(aB_slot),
            fB(aB){
    }
    virtual void print_data(const Particle &p, std::ostream& aOut){
        Logger::print_data(p, aOut);
        MagneticField* pB = fB;
        if(fB_slot>=0){
            if(!pB)
                pB = (MagneticField*)(p.interactionData[fB_slot]);
        }
        double Bperp = 0.;
        if(pB){
            std::vector<double> B0(3);
            pB->GetValueGauss(p.X, p.Time, B0);
            Bperp = sqrt(B0[0]*B0[0]+B0[1]*B0[1]);
        }
        aOut << "\t" << Bperp;
    }

    virtual void print_header(std::ostream& aOut){
        Logger::print_header(aOut);
        aOut << "\tB_perp/G";
    }
};

CRbeam::CRbeam(int argc, char** argv):
		fNumThreads(0),
		fEmin_eV(0),
		fEmax_eV(0),
		fPowerLaw(0),
		fEminSource_eV(0),
		fZmax(0),
		fZmin(0),
		fAlpha(0),
		fEThinning(0),
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
		fPMFFile(0),
		fOverwriteOutput(false),
		fFixedCmb(false),
		fNoCmb(false),
		fEGMF(0.),
		fLcorEGMF(1.),
		fRandomizeEGMF(false),
		fTurbulentEGMF(false),
		fTurbulentFFT(false),
		//fMaxDeflection(0.5/180.*M_PI),
		fLogging(false),
        fTrajectoryLogging(false),
		fLogThread(-1),
		fLogLevel(WarningLL),
		fEBLMaxE(100.),
        fTauPrint(false),
        fJetAngle(0),
        fJetTheta(0),
        fJetPhi(0),
        fJetType(0),
        fJetParam(0),
		fObserverType(-1),
		fObserverX(0),
		fObserverY(0),
		fObserverZ(0),
		fObserverRadius(100),
		fSaveSourceDir(false),
		fSaveMagnetic(false),
		fRandSeed(-1),
		fTraceParticle(false),
		fAccuracy(10),
		fPMFCut(-1),
		fSourceX(0),
		fSourceY(0),
		fSourceZ(0),
		fInteractions(false),
//		fObserverZ(10);
		cmd(argc,argv,commands)
{
	if( cmd.has_param("-h") )
	{
		cmd.printHelp(std::cerr);
		std::cerr << std::endl << "Build on " << __DATE__ << " " << __TIME__ << std::endl;
		exit(255);//enable binary_check
	}
	fNumThreads = cmd(PNumThreads);
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
	fZmin = cmd(PZmin);
	fAlpha = cmd(PThinning);
	fEThinning = cmd(PEThinning);
	fPPPrescaleCoef = cmd(PRescalePPP);

	fEGMF = cmd(PEGMF);
	fLcorEGMF = cmd(PEGMFL);
	fRandomizeEGMF = cmd(PRandomEGMF);
	fTurbulentEGMF = cmd(PTurbulentEGMF);
	fTurbulentFFT = cmd(PTurbulentFFT);
	fTurbulentEGMFstep = cmd(PTurbulentEGMFstep);
	fTurbulentEGMFMinVarScale = cmd(PTurbulentEGMFMinVarScale);
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
	fNoCmb = cmd(PNoCmb);
	fPMFFile = cmd(PMFFile);
	fOutputDir = cmd(POutput);
	fOverwriteOutput = cmd(POverwriteOutput);
	fExtDeltaZconst = cmd(PExtDeltaZconst);
	fExtPowerLow = cmd(PExtPowerLow);
	fExtDeltaZexp = cmd(PExtDeltaZexp);
	fEBLMaxE = cmd(PEBLCut);
    fTauPrint = cmd(PTauPrint);
    fJetAngle = cmd(PJetAngle);
    fJetTheta = cmd(PJetTheta);
    fJetPhi = cmd(PJetPhi);
    fJetType = cmd(PJetType);
    fJetParam = cmd(PJetParam);
	fObserverType = cmd(PObserver);
	fObserverX = cmd(PObserverX);
	fObserverY = cmd(PObserverY);
	fObserverZ = cmd(PObserverZ);
	fObserverRadius = cmd(PObserverRadius);
	fSaveSourceDir = cmd(PSaveSourceDir);
	fSaveMagnetic = cmd(PSaveMagnetic);
	fRandSeed = cmd(PRandSeed);
	fTraceParticle = cmd(PTraceParticle);
	fAccuracy = cmd(PAccuracy);
	fPMFCut = cmd(PMFCut);
	fSourceX = cmd(PSourceX);
	fSourceY = cmd(PSourceY);
	fSourceZ = cmd(PSourceZ);
	fInteractions = cmd(PInteractions);

	if(fOutputDir[0]=='\0')
		fOutputDir = 0;

	if(primary<Electron || primary>Proton)
		Exception::Throw("invalid primary particle type");
	fPrimary = (ParticleType)primary;
}

int CRbeam::run()
{
	double aOmegaVac=0.7;
	double aHubbleKm_sec_Mpc=70;
	int kAc = 1;
	fEBLMaxE *= units.eV;
	double Emax = fEmax_eV*units.eV;//eV
	if(!cosmology.IsInitialized())
		cosmology.Init(fZmax + 10, aOmegaVac, aHubbleKm_sec_Mpc);
	double zMin = fZmin;
	double logStepK = pow(10,0.05/kAc);
	double Emin = fEmin_eV*units.eV;
	double EminSource = fEminSource_eV*units.eV;
	double alphaThinning = fAlpha;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	double EThinning = fEThinning*units.eV;
	std::string outputDir;

//	test_fftw(); // TEST FFTW

	if(fOutputDir)
		outputDir = fOutputDir;
	else
	{
		outputDir = Particle::Name(fPrimary);
		outputDir += ("_E" + ToString(fEmax_eV) + "_z" + ToString(fZmax) + "_N" + ToString(fNoParticles)) + "_EBL" + ToString(fBackgroundModel);
        if(strlen(fCustomBackground)){
            outputDir += 'c';
        }
		if(fEGMF!=0)
		{
			outputDir += "_B";
			if(fRandomizeEGMF)
				outputDir += "rand";
			if(fTurbulentEGMF)
				outputDir += "Turb";
			outputDir += ToString(fEGMF);
            outputDir += "_Lc";
            outputDir += ToString(fLcorEGMF);
		}
		else
		{
			outputDir += "_B0";
		}
		if(fNoEM)
			outputDir += "_noEM";
		if(fMonoCMB)
			outputDir += "_mono";
		if(fFixedCmb)
			outputDir += "_fixedCmb";
		if(fPowerLaw>0)
			outputDir += ("_p" + ToString(fPowerLaw));
		if(fMz!=M_POINT_SOURCE)
			outputDir += ("_m" + ToString(fMz));
	}
    std::cout << "output dir: " << outputDir << std::endl;
    SmartPtr<RawOutput3D> pOutput = new RawOutput3D(outputDir, fOverwriteOutput, fMz!=M_POINT_SOURCE, fPowerLaw>0, fTrajectoryLogging, fObserverType!=-1,fSaveSourceDir>0);//this creates folder outputDir, later we will set output to outputDir/z0
	if(fLogging){
		debug.SetOutputFile(outputDir + "/log.txt");
		debug.EnableTimestamp();
		debug.SetThreadFilter(fLogThread);
		debug.SetLevel(fLogLevel);
	}
	CompoundBackground backgr;

	double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	// if(fFixedCmb)
	// 	cmbTemp *= (1.+fZmax);

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
        case 13:
        	b2 = new Gilmore12fixedEBL();
        	break;
        case 14:
        	b2 = new Gilmore12fiducialEBL();
        	break;
        case 15:
        	b2 = new Gilmore12fixedCiberEBL();
        	break;
        case 16:
        	b2 = new Gilmore12fiducialCiberEBL();
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

	if(!fNoCmb)
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

	unsigned long int seed = 0;
#ifdef _DEBUG
	seed=2015;
#else
	{
	    struct timeval tv;
	    gettimeofday(&tv, NULL);
	    seed = tv.tv_sec * 1000 + tv.tv_usec / 1000;
	}
#endif
	if(fRandSeed>=0){
		seed=fRandSeed;
		cerr<<"Setup seed: "<<seed<<endl;
	}
	Randomizer rand(seed);
	std::vector<int> FilterParticleType;
	std::vector<double> FilterParticleEnergy;

	// FilterParticleType.push_back(0);
	// FilterParticleEnergy.push_back((1-1e-10)*Emax);

	// FilterParticleType.push_back(1);
	// FilterParticleEnergy.push_back((1-1e-10)*Emax);

	// FilterParticleType.push_back(2);
	// FilterParticleEnergy.push_back(1*units.eV);

	// SmartPtr<IFilter> resultFilter = new TypeSpecificEnergyBasedFilter(FilterParticleType, FilterParticleEnergy, M_PI);
	SmartPtr<IFilter> resultFilter = new EnergyBasedFilter(Emin, M_PI);
	bool useResultFilterInRuntime = true;
	Result result(resultFilter, useResultFilterInRuntime);
	result.SetAutoFlushInterval(10);//save every 10 particle cascades
	result.EnableFlushOnCtrlC();//make sure results are not lost if program is interrupted by Ctrl-C
	result.AddOutput(pOutput);
	ParticleStack particles;
	PropagationEngine pe(particles, result, rand.CreateIndependent());

	// EnergyBasedThinning thinning(alphaThinning, BlackList);
	// thinning.EnableFor(Electron);
	// thinning.EnableFor(Positron);

	PhotonThinning thinning(alphaThinning, EThinning);
	pe.SetThinning(&thinning);

	// pe.AddInteraction(new RedShift());

	if(fNoEM) {
		ParticleType discardedParticles[] = {Electron, Positron, Photon, ParticleTypeEOF};
		pe.SetPropagationFilter(new ParticleTypeFilter(discardedParticles));
	}
	else
	{
		pe.AddInteraction(new Interactions::ICS(backgrI, Emin));
		// pe.AddInteraction(new Interactions::IcsCEL(backgrI, Emin));
		// pe.AddInteraction(new Interactions::GammaPP(backgrI));
		// pe.AddInteraction(new PPP(backgrI, fPPPrescaleCoef));//secondaries only
	}
	if(fInteractions){
		pe.AddInteraction(new GZK(backgrI));
		pe.AddInteraction(new NeutronDecay());
		pe.AddInteraction(new ProtonPPcel(backgrI));
	}


	double Beta = 0.;
	double Alpha = 0;
	double Theta = 0;
	double Phi = 0;
	int Bslot = -1;
    SmartPtr<MagneticField> mf,mf2,mf3;
    SmartPtr<ComboMF> cmf;
	if(fEGMF>0. || (strlen(fPMFFile)>0 && !fTurbulentEGMF )) {
		if (fRandomizeEGMF) {
			Deflection3D* defl = new Deflection3D(fAccuracy);
			Bslot = defl->MFSlot();
			pe.AddInteraction(defl);
		}
		else {
			if (fTurbulentEGMF){
//				cout<<"cosmology.z2t0(0.00235)/units.Mpc = "<<cosmology.z2t0(0.00235)/units.Mpc<<endl;
				if(fTurbulentFFT){
					double cube_size = fabs(cosmology.z2t0(fZmax)/units.Mpc)*1.01;
					cout<<"Generate FFT cube for "<<cube_size<<" Mpc max distance"<<endl;
					mf2 = ((MagneticField*) new FFT_TurbulentMF(rand, fLcorEGMF, fEGMF, cube_size, fTurbulentEGMFMinVarScale, fNumThreads,
							fSourceX/1000,fSourceY/1000,fSourceZ/1000));
					cmf = new ComboMF(mf2); // EXPERIMENTAL FEATURE
					//mf3 = ((MagneticField*)new FFT_TurbulentMF(rand, 100*1e-6, 3e-6, 10*1e-3, fTurbulentEGMFMinVarScale, fNumThreads, fSourceX/1000,fSourceY/1000,fSourceZ/1000));
					//mf3 = ((MagneticField*)new TurbulentMF(rand, 100*1e-6, 3e-6, 1.02, fTurbulentEGMFMinVarScale,100));
					//cmf->AddField(mf3);
					mf = (MagneticField*) cmf;
				}else
					mf = ((MagneticField*) new TurbulentMF(rand, fLcorEGMF, fEGMF, fTurbulentEGMFstep, fTurbulentEGMFMinVarScale,fPMFCut));
			}else if( strlen(fPMFFile)>0 ){ // Load MF from file: x y z Bx By Bz columns
			    std::cout<<"Usage input magnetic field from file: "<<fPMFFile<<std::endl;
				mf = ((MagneticField*) new SavedMF(fPMFFile,SavedMF::FILLING_MIRROR));
			}else
				// mf = ((MagneticField*) new CellMF(rand, fLcorEGMF, fEGMF));
				mf = ((MagneticField*) new MonochromaticMF(rand, fLcorEGMF, fEGMF));
				// mf = ((MagneticField*) new MonochromaticMF(fLcorEGMF, fEGMF, Beta, Alpha, Theta, Phi));
			pe.AddInteraction(new Deflection3D(mf,fAccuracy));
			mf->get_div(0,0,0);
			mf->get_div(1,0,0);
			mf->get_div(0,1,0);
			mf->get_div(0,0,1);
			mf->get_div(-1,0,0);
			mf->get_div(0,-1,0);
			mf->get_div(0,0,-1);
			double div1 = fabs(mf->get_integral_div(0,0,0,1,1,1,0.1));
			double div2 = fabs(mf->get_integral_div(0,0,0,1,1,1,0.05));
			double div3 = fabs(mf->get_integral_div(0,0,0,1,1,1,0.009));
			if(div1<div3 || div3>1){
				cerr<<"Div B selftest failed! Stop execution"<<endl;
				exit(250);
			}
		}
	}

	// generating initial state

	CosmoTime tEnd;
	tEnd.setZ(zMin);
	pOutput->SetOutputDir(outputDir + "/z" + ToString(zMin));
	fstream paramsOut((outputDir + "/params").c_str(), std::ios_base::out);
	paramsOut<<"# Build at "<<__DATE__ << " " << __TIME__ << std::endl;
	cmd.printParamValues(paramsOut);
	if(fRandSeed<0){
		paramsOut<<"--random-seed\t"<<seed<<std::endl;
		paramsOut.flush();
	}

	result.SetEndTime(tEnd);
	CosmoTime tStart;
	tStart.setZ(fZmax);
	double dt = tEnd.t()-tStart.t();
	double t0 = tStart.t();

 // Uncomment to check MF realization ////////////////////////////
	if(fSaveMagnetic && !fRandomizeEGMF){
	 	cout<<"Save magnetic field to file..."<<endl;
		std::ofstream mfout((outputDir + "/magneticField.out").c_str());
		double stepScale = mf->MinVariabilityScale(tEnd)/fAccuracy/units.Mpc * 4;
		std::cout << "Old scale: " << stepScale << std::endl;
		stepScale = 0.05;
		std::cout << "New scale: " << stepScale << std::endl;
		mf->Print(stepScale,dt/units.Mpc, mfout);
		mfout.close();
	}else if(fSaveMagnetic){
		cout<<"Saving magnetic field to file is not possible for --randomEGMF key!"<<endl;
	}
	if(!fRandomizeEGMF){
		std::ofstream mfh((outputDir + "/magneticField.hash").c_str());
		mf->PrintUniqHash(mfh);
		mfh.close();
	}
 // ///////////////////////////////// ////////////////////////////

	if(fObserverType!=-1){
		/*if(fObserverZ==0){ // For backward compatibility
			fObserverZ=dt/units.kpc;
			cerr<<"Note: observer Z coords set as not zero, observer will be on shell:"<<dt<<endl;
		}//*/
		pe.fResult.fObserver = Observer(fObserverType,  fObserverX*units.kpc,
														fObserverY*units.kpc,
														fObserverZ*units.kpc,
														fObserverRadius*units.kpc);// for update dt
		result.fObserver = Observer(fObserverType,  fObserverX*units.kpc,
													fObserverY*units.kpc,
													fObserverZ*units.kpc,
													fObserverRadius*units.kpc);// for check hit
	}

	double normW = 1.;
	if(fPowerLaw>0 && fPowerLaw!=1.){//make sure total weight is roughly equal to number of particles
		normW = log(Emax/EminSource)*(1.-fPowerLaw)/(pow(Emax,1.-fPowerLaw)-pow(EminSource,1.-fPowerLaw));
	}
    SafePtr<CRbeamLogger> tr_logger;
    std::ofstream trLogOut;
    if(fTrajectoryLogging){
        trLogOut.open((outputDir + "/tr.log").c_str());
        if(!trLogOut.is_open()){
            Exception::Throw("Failed to create " + outputDir + "/tr.log");
        }

        tr_logger = new CRbeamLogger(trLogOut, Bslot, mf);
        pe.SetLogger(tr_logger);
		fBatchSize = fNoParticles; // this is needed to make sure particle Id's are unique
    }
	if(fBatchSize==0){
		fBatchSize = omp_get_max_threads()*10;
		std::cout << "max threads = " << omp_get_max_threads() << std::endl;
		std::cerr << "using auto batch size " << fBatchSize << std::endl;
	}

    std::vector<SafePtr<MagneticField> > fields(fBatchSize);
    if(fTauPrint){
        class Kern : public Function{
            PropagationEngine& fPe;
            mutable Particle   fParticle;
        public:
            Kern(PropagationEngine& pe, const Particle aParticle):
                    fPe(pe),
                    fParticle(aParticle){};
            virtual double f(double t) const{
                fParticle.Time = t;
                return fPe.GetInteractionRate(fParticle);
            }
        };
        CosmoTime tSource(fZmax);
        //CosmoTime tEnd(0);
        if(fPowerLaw<=0.)
            EminSource = Emax;
        double step = pow(10.,0.1);
        MathUtils mu;
        ofstream out_tau((outputDir + "/tau").c_str());
        for (double E=EminSource; E<=Emax; E*=step){
            Particle particle(fPrimary, tSource.z());
            particle.Energy = E;
            Kern k(pe, particle);
            double tau = mu.Integration_qag(k, tSource.t(),tEnd.t(), 0, 1e-3, 1000);
            out_tau << E/units.eV << "\t" << tau << std::endl;
        }
        out_tau.close();
    }
	if(fSourceX*fSourceX+fSourceY*fSourceY+fSourceZ*fSourceZ > 0.1){
		std::cerr<<"Warning: you are using not tested feature. Check redshift and Universe time before analyse"<<std::endl;
	}
	for(int i=0; i<fNoParticles;)
	{
		CosmoTime tSource(fZmax);
		double weight = 1.;
		if(fMz!=M_POINT_SOURCE)
		{
			double t = t0 + rand.Rand()*dt;
			tSource.setT(t);
			double z = tSource.z();
			weight = pow(1.+z,fMz);
		}
		Particle particle(fPrimary, tSource.z());
		if(fPowerLaw<=0.){
			particle.Energy = Emax;
			particle.Weight = weight;
		}
		else{
			particle.Energy = EminSource*pow(Emax/EminSource, rand.Rand());//choose E randomly in log scale
			particle.Weight = weight*normW*pow(particle.Energy, 1.-fPowerLaw);
		}
		/*********************/
		if(fTraceParticle){
			//std::ofstream pfout((outputDir + "/particle"+ToString(i)+".dat").c_str());
			particle.Tracer=new std::ofstream((outputDir + "/particle"+ToString(i)+".dat").c_str());
		}
		particle.X[0]=fSourceX*units.kpc;
		particle.X[1]=fSourceY*units.kpc;
		particle.X[2]=fSourceZ*units.kpc;
		if(fJetAngle>1.e-14){// Customize source direction
			double cos_theta_max = cos(fJetAngle);
			double rndval = rand.RandZero();
			double cos_theta = 1.-rndval*(1.-cos_theta_max); // Random from cos=1...maxcos
			if(fJetType==1){ // TODO: Bad codestyle, rewrite it
				cos_theta = 1.-fabs(rand.RandGauss(fJetParam));
				if(cos_theta<-1.){
					cos_theta+=2.; // Hot fix
					cerr<<"Warning: bad random fix for Gauss distribution"<<endl;
				}
			}
			double sin_theta = sqrt(1.-cos_theta*cos_theta);
			double phi = rand.Rand()*2.*M_PI; // Random from 0 to 2pi
			double Px = cos(phi)*sin_theta;
			double Py = sin(phi)*sin_theta;
			double Pz = cos_theta;
			// Now rotate it:
			// By theta:
			particle.Pdir[0] = Px;
			particle.Pdir[1] = cos(fJetTheta)*Py-sin(fJetTheta)*Pz;
			particle.Pdir[2] = sin(fJetTheta)*Py+cos(fJetTheta)*Pz;
			// By phi:
			fJetPhi+=M_PI/2.;
			particle.Pdir[0] = cos(fJetPhi)*Px-sin(fJetPhi)*particle.Pdir[1];
			particle.Pdir[1] = sin(fJetPhi)*Px+cos(fJetPhi)*particle.Pdir[1];
			fJetPhi-=M_PI/2.;
			// Copy to start direction, this place don't need high speed
			particle.PdirStart[0] = particle.Pdir[0];
			particle.PdirStart[1] = particle.Pdir[1];
			particle.PdirStart[2] = particle.Pdir[2];
			switch(fJetType){
				case 0:
					particle.Weight *= coneZ_weight(cos_theta,cos_theta_max); // 1.
				break;
				case 1:
					particle.Weight *= 1./(1.e-10+GaussConeZ_weight(fJetTheta, fJetPhi,particle.Pdir,fJetParam));
					if(isnan(particle.Weight)){
						particle.Weight=0;
//						cerr<<"WTF?!"<<endl;
					}
				break;
				case 2:
					particle.Weight *= vMF_weight(fJetTheta, fJetPhi,particle.Pdir,fJetParam);
				break;
				default:
					Exception::Throw("Invalid jet type. Maybe this is a future feature?");
			}
		}
		/********************/
		if(fRandomizeEGMF)
		{
            MagneticField* mf = fTurbulentEGMF ? ((MagneticField*) new TurbulentMF(rand, fLcorEGMF, fEGMF, fTurbulentEGMFstep, fTurbulentEGMFMinVarScale,fPMFCut) /*new TurbulentMF(rand, fLcorEGMF, fEGMF)*/) :
						((MagneticField*) new CellMF(rand, fLcorEGMF, fEGMF));
			if(i==0){
				std::ofstream mfh((outputDir + "/magneticField.hash").c_str());
				mfh<<"# Only first realization, there is --randomEGMF key!"<<std::endl;
				mf->PrintUniqHash(mfh);
				mfh.close();
			}
            fields[i%fBatchSize] = mf; //store till the end of current batch
			particle.interactionData[Bslot] = mf;
		}
		particles.AddPrimary(particle);
        i++;
        if(i%fBatchSize==0 || i==fNoParticles){
            std::cerr << "batch " << (i-1)/fBatchSize + 1 << " of " << ceil(fNoParticles/(double)fBatchSize) << std::endl;
            pe.RunMultithread(0, fNumThreads);
            result.Flush();
            if(result.AbortRequested())
                break; //handling Ctrl-C
        }
	}

	return 0;
}

CRbeam::~CRbeam() {

}
