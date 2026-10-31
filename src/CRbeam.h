/*
 * PointUHECR.h
 *
 *  Created on: Aug 10, 2015
 *      Author: ok
 */

#ifndef CRBEAM_H_
#define CRBEAM_H_

#include <cstdlib>
#include <iostream>
#include "MCPropag/ParticleStack.h"
#include "MCPropag/Utils.h"
#include "MCPropag/Cosmology.h"
#include "MCPropag/PropagationEngine.h"
#include "MCPropag/Output.h"
#include "MCPropag/TableBackgrounds.h"
#include "MCPropag/ICS.h"
#include "MCPropag/GammaPP.h"
#include "CmdLine.h"

using namespace std;
using namespace MCPropag;
using namespace Utils;
using namespace Backgrounds;
using namespace Interactions;

class CRbeam {
public:
	CRbeam(int argc, char** argv);
	virtual ~CRbeam();
	void Test();
	int run();
private:
	int 			fNumThreads;
	double			fEmin_eV;
	double			fEmax_eV;
	double			fPowerLaw;
	double	 		fEminSource_eV;
	double			fZmax;
	double			fZmin;
	double			fAlpha;
	double			fEThinning;
	int 			fPPPrescaleCoef;
	int				fNoParticles;
	int				fBatchSize;
	ParticleType 	fPrimary;
	bool 			fNoEM;
	double			fMz;
	int				fBackgroundModel;
    double          fEBLmult;
	const char*		fCustomBackground;
	bool 			fCustomBackgroundPhysical;
	bool			fMonoCMB;
	double			fExtDeltaZconst;
	double			fExtPowerLow;
	double			fExtDeltaZexp;
	const char*		fOutputDir;
	const char*		fPMFFile;
	bool 			fOverwriteOutput;
	bool 			fSaveSourceDir;
	bool            fSaveMagnetic;
	bool			fFixedCmb;
	bool			fNoCmb;
	double			fEGMF;
	double			fLcorEGMF;
	bool 			fRandomizeEGMF;
	bool 			fTurbulentEGMF;
	bool			fTurbulentFFT;
	double 			fTurbulentEGMFstep;
	double 			fTurbulentEGMFMinVarScale;
	bool            fLogging;
	bool 			fTrajectoryLogging;
	int             fLogThread;
	LogLevel        fLogLevel;
	double			fEBLMaxE;
	double			fJetAngle;
	double			fJetTheta;
	double			fJetPhi;
	double			fJetParam;
	int				fJetType;
	int				fObserverType;
	double			fObserverX;
	double			fObserverY;
	double			fObserverZ;
	double			fObserverRadius;
//	double			fObserverZ;
    bool            fTauPrint;
    int             fRandSeed;
    bool            fTraceParticle;
    int				fAccuracy;
    double			fPMFCut;
	double			fSourceX;
	double			fSourceY;
	double			fSourceZ;
    bool            fInteractions;
	//double			fMaxDeflection;
	cors::cmdline::CmdLine 		cmd;
};

#endif /* CRBEAM_H_ */
