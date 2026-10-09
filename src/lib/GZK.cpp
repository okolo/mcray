/*
 * GZK.cpp
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

#include "GZK.h"
#include "TableFunction.h"
#include "Sophia.h"

//used for unit test only
#include "PropagationEngine.h"
#include "Inoue12IROSpectrum.h"
#include "TableBackgrounds.h"
#include "NeutronDecay.h"
#include "Test.h"

using namespace mcray;
namespace Interactions{

GZK::GZK(BackgroundIntegral* aBackground, int aRandSeed,
         std::shared_ptr<const NuclearDecayData> nuclearData)
    :fBackground(aBackground),fNuclearData(std::move(nuclearData))
{
    fMpSophia = SOPHIA::Mass(Proton);
    fMnSophia = SOPHIA::Mass(Neutron);
    fSigmaN = InitSigma(Neutron);
    fSigmaP = InitSigma(Proton);
    if(aRandSeed)
    	SOPHIA::SetRandomSeed(aRandSeed);
}

Function* GZK::InitSigma(ParticleType aPrim)
{
    std::string tablesDir = tables_dir + "sophia2/";
    Utils::TableReader reader(tablesDir+(aPrim==Proton?"p":"n"), 2);
    std::vector<double>& s = reader.getColumn(0);
    std::vector<double>& sigma = reader.getColumn(1);
    double M = aPrim==Proton?fMpSophia:fMnSophia;
    for(int i=0; i<s.size(); i++)
    {
    	double val = M*(M+2.*s[i]*units.GeV);
        s[i] = M*(M+2.*s[i]*units.GeV);
        sigma[i] *= (units.barn*1e-6);
    }
    double rightVal = sigma[sigma.size()-1];
    return new LinearFunc(s, sigma, new LogScale(), new LogScale(), 0, rightVal);
}

double GZK::Rate(const mcray::Particle &aParticle) const {
    if (aParticle.Type==Nucleus && fNuclearData) {
        const auto rates=NuclearTargetRates(aParticle);
        const double total=rates[0]+rates[1];
        if (!std::isfinite(total)) Utils::Exception::Throw("Invalid nuclear GZK total rate");
        return total;
    }
    if(aParticle.Type!=Proton && aParticle.Type!=Neutron)
        return 0.;
    const Function* sigma = (aParticle.Type==Proton) ? fSigmaP : fSigmaN;
    return fBackground->GetRateS(*sigma, aParticle);
}

RandomInteraction *GZK::Clone() const {
    return new GZK(fBackground->Clone(),0,fNuclearData);
}

double GZK::NuclearModifier(unsigned int a,unsigned int count) {
    if (!a || count>a) Utils::Exception::Throw("Invalid bound-nucleon count");
    if (!count) return 0.;
    if (a==1) return 1.;
    return 0.85*(a<=8 ? std::pow(static_cast<double>(count),2./3.) : count);
}
std::array<double,2> GZK::NuclearTargetRates(const Particle& p) const {
    if (!fNuclearData || p.Type!=Nucleus) return {{0.,0.}};
    p.ValidateIdentity();
    if (p.MassNumber()>56 || p.MassModel()!=Particle::NuclearMassModel::TALYS ||
        p.Mass()!=fNuclearData->Nuclei().Get(p.MassNumber(),p.AtomicNumber()).mass ||
        !std::isfinite(p.Energy) || p.Energy<p.Mass())
        Utils::Exception::Throw("Nuclear GZK requires A<=56 and the selected nuclear-data mass");
    std::array<double,2> rates{};
    for (size_t i=0;i<2;++i) {
        const unsigned int count=i ? p.MassNumber()-p.AtomicNumber() : p.AtomicNumber();
        if (!count) continue;
        Particle target=p; target.SetType(i ? Neutron : Proton);
        target.Energy=p.Energy*(target.Mass()/p.Mass()); // Rate at equal gamma.
        rates[i]=NuclearModifier(p.MassNumber(),count)*Rate(target);
        if (!std::isfinite(rates[i]) || rates[i]<0.) Utils::Exception::Throw("Invalid nuclear GZK target rate");
    }
    return rates;
}
bool GZK::GetSecondaries(Particle& p,std::vector<Particle>& products,Randomizer& random) const {
    if (p.Type!=Nucleus) return RandomInteractionS::GetSecondaries(p,products,random);
    const auto rates=NuclearTargetRates(p);
    const double total=rates[0]+rates[1];
    if (!std::isfinite(total)) Utils::Exception::Throw("Invalid nuclear GZK total rate");
    if (!(total>0.)) return false;
    const bool onProton=random.Rand()*total<rates[0];
    Particle target=p; target.SetType(onProton ? Proton : Neutron);
    target.Energy=p.Energy/p.MassNumber(); // CRPropa's bound-nucleon kinematics.
    double s=0.;
    if (!SampleS(target,s,random)) return false;
    std::vector<Particle> generated;
    SampleSecondaries(target,generated,s,random);
    if (generated.empty()) Utils::Exception::Throw("SOPHIA returned no products for nuclear GZK");
    fNuclearData->AppendResidual(p,p.MassNumber()-1,p.AtomicNumber()-(onProton?1:0),
                                p.Energy-target.Energy,generated);
    products.insert(products.end(),generated.begin(),generated.end());
    return true;
}

bool GZK::SampleS(const mcray::Particle &aParticle, double &aS,
                                mcray::Randomizer &aRandomizer) const {
    if(aParticle.Type!=Proton && aParticle.Type!=Neutron)
        return false;
    const Function* sigma = aParticle.Type==Proton?fSigmaP:fSigmaN;
    if(fBackground->GetRateAndSampleS(*sigma, aParticle, aRandomizer, aS)!=0)
    {
    	const double sophiaThreshold=1.1646*units.GeV*units.GeV;
    	if(aS<sophiaThreshold)
    	{
    		ASSERT(aS>=sophiaThreshold*0.999);
    		aS = 1.001*sophiaThreshold;
    	}
    	return true;
    }
    else
    	return false;
}

void GZK::SampleSecondaries(mcray::Particle &aParticle,
                                          std::vector<mcray::Particle> &aSecondaries, double aS,
                                          mcray::Randomizer &aRandomizer) const {
    double M = aParticle.Type==Proton?fMpSophia:fMnSophia;
    double epsPrimeGeV = 0.5*(aS/M-M)/units.GeV;
    int noSecondaries = 0;
    std::vector<double> secEfrac(SOPHIA::MaxNumberOfProducts);
    std::vector<int> secTypes(SOPHIA::MaxNumberOfProducts);
#pragma omp critical (SOPHIA)
    {//SOPHIA class is not thread-safe
        //todo: separately collect protons and neutrons in two queues and try to handle the queues one after another with critical section disabled
        SOPHIA::SamplePhotopionRel(aParticle.Type, epsPrimeGeV, noSecondaries, secEfrac.data(), secTypes.data());
    }
    if (noSecondaries<0 || noSecondaries>SOPHIA::MaxNumberOfProducts)
        Utils::Exception::Throw("Invalid SOPHIA product count");
    for (int i=0;i<noSecondaries;++i)
        if (secTypes[i]<Electron || secTypes[i]>=Nucleus ||
            !std::isfinite(secEfrac[i]) || secEfrac[i]<0.)
            Utils::Exception::Throw("Unsupported SOPHIA product or invalid energy fraction");
    while(--noSecondaries>=0)
    {
        Particle sec = aParticle;
        sec.SetType((ParticleType)secTypes[noSecondaries]);
        sec.Energy *= secEfrac[noSecondaries];
        sec.fCascadeProductionTime = aParticle.Time;
        aSecondaries.push_back(sec);
    }
}

    void GZK::UnitTest()
    {
        double Zmax = 1;
        double Emin = 1e9*units.eV;
        double Emax = 1e21*units.eV;
        int Nparticles = 10000;
        std::string outputDir = "testGZK";
        unsigned int kAc = 1;
        double epsRel = 1e-3;
        double stepZ = Zmax<0.05 ? Zmax/2 : 0.025;
        double logStepK = pow(10,0.05/kAc);
        if(!cosmology.IsInitialized())
            cosmology.Init(Zmax + 10);
        double cmbTemp = 2.73*units.K;
        double alphaThinning = 0.5;
        CompoundBackground backgr;
        backgr.AddComponent(new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp, 0., Zmax + 1.));
        //IR/O component
        backgr.AddComponent(new Backgrounds::Inoue12BaselineIROSpectrum());  //new GaussianBackground(0.1*units.eV, 0.01*units.eV, 5, 1./units.cm3, 0, Zmax + 1.));
        PBackgroundIntegral backgrI(new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel));
        int seed = 2015;
        Result result(new EnergyBasedFilter(Emin, M_PI), true);

        SmartPtr<RawOutput> pOutput = new RawOutput(outputDir, false);
        result.AddOutput(pOutput);
        ParticleStack particles;
        PropagationEngine pe(particles, result, seed);
        EnergyBasedThinning thinning(alphaThinning);
        pe.SetThinning(&thinning);

        Particle proton(Proton, Zmax);
        int nSteps = log(Emax/Emin)/log(logStepK)+1.;
        double mult = pow(Emax/Emin, 1./nSteps);
        {
            std::ofstream rateOut;
            rateOut.open((outputDir + "/GZK").c_str(),std::ios::out);
            PRandomInteraction i = new GZK(backgrI);
            pe.AddInteraction(i);
            int step=0;
            for(proton.Energy=Emax/100; step<=nSteps; proton.Energy*=mult, step++)
                rateOut << proton.Energy/units.eV << "\t" << i->Rate(proton)*units.Mpc << "\n";
            rateOut.close();
        }
        Randomizer rand;
        CosmoTime tEnd;
        proton.Energy = Emax;
        pOutput->SetOutputDir(outputDir + "/z0");
        result.SetEndTime(tEnd);

        for(int i=0; i<Nparticles; i++)
        {
            particles.AddPrimary(proton);
        }
        pe.RunMultithreadReleaseOnly();
    }

    void GZK::UnitTestEconserv(double aE, double aZ){
        int nParticles = 100;
        Test::EConservationTest test(aZ,1e6,2015);
        test.alphaThinning = 0.0;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
        test.Engine().AddInteraction(new GZK(test.Backgr(),(int)(test.GetRandomizer().CreateIndependent())));
        //test.Engine().AddInteraction(new NeutronDecay());
        Particle p(Proton, aZ);
        p.Energy = aE*units.eV;
        test.SetPrimary(p, nParticles);
        test.Run();
    }
}
