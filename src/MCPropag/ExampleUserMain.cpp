/*
 * ExampleUserMain.cpp
 *
 *  Created on: Oct 13, 2012
 *      Author: ok
 */

#include <cstdlib>
#include "ParticleStack.h"
#include "Utils.h"
#include "Cosmology.h"
#include "PropagationEngine.h"
#include <iostream>
#include "Output.h"
#include "TableBackgrounds.h"
#include "ICS.h"
#include "GammaPP.h"
#include "Test.h"
#include "PrecisionTests.h"

using namespace std;
using namespace MCPropag;
using namespace Utils;
using namespace Backgrounds;
using namespace Interactions;

void BuildInitialState(ParticleStack& particles, double aZmax)
{
	int nParticles = 10000;
	for(int i=0; i<nParticles; i++)
	{
		Particle gamma(Photon, aZmax);
		gamma.Energy = 1e7/units.Eunit;//MeV
		gamma.Weight = 1;
		particles.AddPrimary(gamma);
	}
}

/*
 * This is an example of user_main function performing user defined tasks
 * the function should be provided by end user
 */
int user_main(int argc, char** argv) {
	double Zmax = 0.01;
	double Emin = 1./units.Eunit;
	double alphaThinning = 1;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
	ParticleStack particles;
	Result result(Emin);
	result.AddOutput(new SpectrumOutput("out", Emin, pow(10, 0.05)));

	cosmology.Init();

	//CompoundBackground backgr;
	//double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;
	//backgr.AddComponent(new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp));
	//backgr.AddComponent(new Backgrounds::Kneiske1001EBL());

	//backgr.AddComponent(new Backgrounds::GaussianBackground(1e-6/units.Eunit, 1e-7/units.Eunit, 1., 0.1*units.Vunit));
	//double n = BackgroundUtils::CalcIntegralDensity(backgr)/units.Vunit;
	//backgr.AddComponent(new Backgrounds::GaussianBackground(6.3e-10/units.Eunit, 1e-11/units.Eunit, 1., 413*units.Vunit));
	//n = BackgroundUtils::CalcIntegralDensity(backgr)/units.Vunit;

	GaussianBackground b1(1e-6/units.Eunit, 1e-7/units.Eunit, 1., 0.1*units.Vunit);
	GaussianBackground b2(6.3e-10/units.Eunit, 1e-11/units.Eunit, 1., 413*units.Vunit);


	double stepZ = Zmax<0.2 ? Zmax/2 : 0.1;
	double logStepK = pow(10,0.05);
	double epsRel = 1e-3;

	//BackgroundIntegral backgrI(backgr, stepZ, logStepK, Zmax, epsRel);
	PBackgroundIntegral backgrI1 = new ContinuousBackgroundIntegral(b1, stepZ, logStepK, Zmax, epsRel);
	PBackgroundIntegral backgrI2 = new ContinuousBackgroundIntegral(b1, stepZ, logStepK, Zmax, epsRel);


	PropagationEngine pe(particles, result, 2011);
	EnergyBasedThinning thinning(alphaThinning);
	pe.SetThinning(&thinning);
	pe.AddInteraction(new RedShift());
	pe.AddInteraction(new Interactions::ICS(backgrI1,Emin));
	pe.AddInteraction(new Interactions::IcsCEL(backgrI1,Emin));
	pe.AddInteraction(new Interactions::GammaPP(backgrI1));
	pe.AddInteraction(new Interactions::ICS(backgrI2,Emin));
	pe.AddInteraction(new Interactions::IcsCEL(backgrI2,Emin));
	pe.AddInteraction(new Interactions::GammaPP(backgrI2));

	/// build initial state
	BuildInitialState(particles, Zmax);

	pe.Run();

    return 0;
}
