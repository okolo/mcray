/*
 * CRbeam.h
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

#ifndef CRBEAM_H_
#define CRBEAM_H_

#include <cstdlib>
#include <iostream>
#include "ParticleStack.h"
#include "Utils.h"
#include "Cosmology.h"
#include "PropagationEngine.h"
#include "Output.h"
#include "TableBackgrounds.h"
#include "ICS.h"
#include "GammaPP.h"
#include "CmdLine.h"
#include "Configuration.h"

using namespace std;
using namespace mcray;
using namespace Utils;
using namespace Backgrounds;
using namespace Interactions;

class CRbeam {
	class FilamentFilter : public Result{
	public:
		FilamentFilter(cosmo_time aRedshift, IFilter*	aOutputFilter, bool aIsPropagationFilter):
				Result(aOutputFilter, aIsPropagationFilter),
				fFilamentLocation(aRedshift)
				{}
		virtual cosmo_time GetMinTravelTimeLeft(const Particle& aParticles) const;
	private:
		CosmoTime fFilamentLocation;
	};

public:
	CRbeam(int argc, char** argv);
	virtual ~CRbeam();
	void Test();
	int run();
    bool IsNuclearPrimary() const { return fNuclearPrimary; }
    unsigned int PrimaryA() const { return fPrimaryA; }
    unsigned int PrimaryZ() const { return fPrimaryZ; }
    bool HasExplicitSeed() const { return fExplicitSeed; }
    unsigned long ConfiguredSeed() const { return fSeed; }
    const YAML::Node& EffectiveConfiguration() const { return fConfiguration.Effective(); }
    const std::string& PDMissingDataPolicy() const { return fPDMissingData; }
    const std::string& PDPhotonMode() const { return fPDPhotons; }
    const std::string& DecayPhotonMode() const { return fDecayPhotons; }
    const std::string& PDPhotonYieldPolicy() const { return fPDPhotonYields; }
    const std::string& NuclearGZKMode() const { return fNuclearGZK; }
    const std::string& GZKResidualPolicy() const { return fGzkResiduals; }
    const std::string& NuclearPPPMode() const { return fNuclearPPP; }
    const std::string& PhotodisintegrationMode() const { return fPhotodisintegration; }
private:
	double			fEmin_eV;
	double			fEmax_eV;
	double			fPowerLaw;
	double	 		fEminSource_eV;
	double			fZmax;
    double          fZfilament;
	double			fAlpha;
	double 			fPPPrescaleCoef;
	int				fNoParticles;
	int				fBatchSize;
	ParticleType 	fPrimary;
    bool fNuclearPrimary;
    unsigned int fPrimaryA, fPrimaryZ;
    bool fExplicitSeed;
    unsigned long fSeed;
    std::string fPhotodisintegration, fTalysDataPath, fNuclearDecay, fNuclearPPP, fNuclearGZK, fGzkResiduals, fPDPhotons, fDecayPhotons, fPDPhotonYields, fPDMissingData;
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
	bool 			fOverwriteOutput;
	bool			fFixedCmb;
	double			fEGMF;
	double			fLminEGMF;
	double			fLmaxEGMF;
	double			fNmodesEGMF;
	bool 			fRandomizeEGMF;
	bool 			fTurbulentEGMF;
	uint            fDeflectionAccuracy;
	bool            fLogging;
	bool 			fTrajectoryLogging;
	int             fLogThread;
	LogLevel        fLogLevel;
	double			fEBLMaxE;
    bool            fTauPrint;
	//double			fMaxDeflection;
    std::string fResolvedCustomBackground;
    Configuration fConfiguration;
	cors::cmdline::CmdLine 		cmd;
};

#endif /* CRBEAM_H_ */
