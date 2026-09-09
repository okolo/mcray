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
#include "Randomizer.h"

//used for unit test only
#include "PropagationEngine.h"
#include "Inoue12IROSpectrum.h"
#include "TableBackgrounds.h"
#include "NeutronDecay.h"
#include "Test.h"

using namespace mcray;
namespace Interactions{

GZK::GZK(BackgroundIntegral* aBackground, int aRandSeed):fBackground(aBackground)
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
    std::string tablesDir = "tables/sophia2/";
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
    if(aParticle.Type==Proton || aParticle.Type==Neutron)
    {
        const Function* sigma = (aParticle.Type==Proton) ? fSigmaP : fSigmaN;
        return fBackground->GetRateS(*sigma, aParticle);
    }
    else if (aParticle.Type >= H2 && aParticle.Type < ParticleTypeEOF)
    {
        Particle innerProton = aParticle;
        innerProton.Type = Proton;
        innerProton.Energy *= (innerProton.Mass()/aParticle.Mass());
        double ProtonRate = fBackground->GetRateS(*fSigmaP, innerProton);

        Particle innerNeutron = aParticle;
        innerNeutron.Type = Neutron;
        innerNeutron.Energy *= (innerNeutron.Mass()/aParticle.Mass());
        double NeutronRate = fBackground->GetRateS(*fSigmaN, innerNeutron);

        int aPrimZ = Particle::ElectricCharge(aParticle.Type);
        int aPrimN = Particle::AtomicMass(aParticle.Type)-aPrimZ;

        //std::cout<<"GZK:\t"<<aParticle.Type<<"\t"<<ProtonRate<<"\t"<<NeutronRate<<"\t"<<aPrimZ*ProtonRate+aPrimN*NeutronRate<<std::endl;
        return (aPrimZ*ProtonRate + aPrimN*NeutronRate);
    }
    else
    {
        return 0;
    }
}

RandomInteraction *GZK::Clone() const {
    return new GZK(fBackground->Clone());
}

bool GZK::SampleS(const mcray::Particle &aParticle, double &aS,
                                mcray::Randomizer &aRandomizer) const {
    ParticleType aType;
    if(aParticle.Type==Proton || aParticle.Type==Neutron)
    {
        aType = aParticle.Type;
        const Function* sigma = aType==Proton?fSigmaP:fSigmaN;
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
    else if (aParticle.Type >= H2 && aParticle.Type < ParticleTypeEOF)
    {
        int aPrimZ = Particle::ElectricCharge(aParticle.Type);
        double Proportion = aPrimZ/Particle::AtomicMass(aParticle.Type)

        double aRand = aRandomizer.Rand();
        aType = aRand<Proportion?Proton:Neutron;

        Particle innerNucleon1 = aParticle;
        innerNucleon1.Type = aType;
        innerNucleon1.Energy *= (innerNucleon1.Mass()/aParticle.Mass());
        const Function* sigma1 = aType==Proton?fSigmaP:fSigmaN;

        Particle innerNucleon2 = aParticle;
        innerNucleon2.Type = aType;
        innerNucleon2.Energy *= (innerNucleon2.Mass()/aParticle.Mass());
        const Function* sigma2 = aType==Proton?fSigmaP:fSigmaN;


        if(fBackground->GetRateAndSampleS(*sigma1, innerNucleon1, aRandomizer, aS)!=0)
        {
            const double sophiaThreshold=1.1646*units.GeV*units.GeV;
            if(aS<sophiaThreshold)
            {
                ASSERT(aS>=sophiaThreshold*0.999);
                aS = 1.001*sophiaThreshold;
            }
            aS *= aType==Proton?1:-1;
            return true;
        }
        else if(fBackground->GetRateAndSampleS(*sigma2, innerNucleon2, aRandomizer, aS)!=0)
        {
            const double sophiaThreshold=1.1646*units.GeV*units.GeV;
            if(aS<sophiaThreshold)
            {
                ASSERT(aS>=sophiaThreshold*0.999);
                aS = 1.001*sophiaThreshold;
            }
            aS *= aType==Proton?1:-1;
            return true;
        }
        else
            return false;
    }
    else
    {
        return false;
    }
}

void GZK::SampleSecondaries(mcray::Particle &aParticle,
                                          std::vector<mcray::Particle> &aSecondaries, double aS,
                                          mcray::Randomizer &aRandomizer) const {
    if(aParticle.Type==Proton || aParticle.Type==Neutron)
    {
        double M = aParticle.Type==Proton?fMpSophia:fMnSophia;
        double epsPrimeGeV = 0.5*(aS/M-M)/units.GeV;
        int noSecondaries = 0;
        double secEfrac[250];
        int secTypes[250];
    #pragma omp critical (SOPHIA)
        {//SOPHIA class is not thread-safe
            //todo: separately collect protons and neutrons in two queues and try to handle the queues one after another with critical section disabled
            SOPHIA::SamplePhotopionRel(aParticle.Type, epsPrimeGeV, noSecondaries, secEfrac, secTypes);
        }
        while(--noSecondaries>=0)
        {
            Particle sec = aParticle;
            sec.Type = (ParticleType)secTypes[noSecondaries];
            sec.Energy *= secEfrac[noSecondaries];
            sec.fCascadeProductionTime = aParticle.Time;
            aSecondaries.push_back(sec);
        }
    }
    else if (aParticle.Type >= H2 && aParticle.Type < ParticleTypeEOF)
    {   
        ParticleType aType;
        double M;
        int aPrimZ = aParticle.ElectricCharge();
        int aPrimN = aParticle.AtomicMass()-aPrimZ;
        int dZ = 0;
        int dN = 0;

        if(aS<0)
        {
            aS *= -1;
            M = fMnSophia;
            aType = Neutron;
            dN = 1;
        }
        else
        {
            M = fMpSophia;
            aType = Proton;
            dZ = 1;
        }

        double NucleonEnergy = aParticle.Energy * M/aParticle.Mass();
        double epsPrimeGeV = 0.5*(aS/M-M)/units.GeV;
        int noSecondaries = 0;
        double secEfrac[250];
        int secTypes[250];

    #pragma omp critical (SOPHIA)
        {//SOPHIA class is not thread-safe
            //todo: separately collect protons and neutrons in two queues and try to handle the queues one after another with critical section disabled
            SOPHIA::SamplePhotopionRel(aType, epsPrimeGeV, noSecondaries, secEfrac, secTypes);
        }
        
        for (int i=0; i<noSecondaries; i++)
        {
            int Type = secTypes[i];
            if (Type==9)
                dN -= 1;
            else if (Type ==10)
                dZ -= 1;
        }

        //std::cout<<dZ<<"\t"<<dN<<std::endl;
        int mainSecZ = aPrimZ-dZ;
        int mainSecN = aPrimN-dN;

        ParticleType mainSecType;
        bool hasMainSec = (mainSecZ>0 && mainSecN>0);
        bool invalidBranch = false;


        if (hasMainSec)
        {
            if (dZ==0 && dN==0)
                mainSecType = aParticle.Type;
            else if (mainSecZ==1 && mainSecN == 0)
                mainSecType = Proton;
            else if (mainSecZ==0 && mainSecN == 1)
                mainSecType = Neutron;
            else
            {
                for (int i=66; ParticleType(i)<EndRealNuclei; i++)
                {
                    int Z = Particle::ElectricCharge(ParticleType(i));
                    int N = Particle::AtomicMass(ParticleType(i))-Z;
                    if (mainSecZ==Z && mainSecN == N)
                    {
                        mainSecType = ParticleType(i);
                        break;
                    }
                    if (i==EndRealNuclei-1)
                    {
                        hasMainSec = false;
                    }

                }
            }
        }
        ASSERT(mainSecType.Type<StartRealNuclei || mainSecType.Type>EndRealNuclei);

        if (!hasMainSec)
        {
            if(mainSecZ>0 || mainSecN>0)
            {
                if(aType==9)
                    std::cout<<std::endl<<"GZK:\tinvalid Neutron channel for "<<aParticle.ToString()<<std::endl;
                else
                    std::cout<<std::endl<<"GZK:\tinvalid Proton channel for "<<aParticle.ToString()<<std::endl;
                std::cout<<"\tMay be secondary particle is an unstable short-living nuclei?"<<std::endl;
                std::cout<<"\tGZK process cancelled"<<std::endl<<std::endl;
                Particle survivedPrim = aParticle;
                survivedPrim.Weight *= 0;                               //allows to exclude all future secondaries, but not previous
                aParticle.Ninteractions -= 1;
                aSecondaries.push_back(survivedPrim);
                invalidBranch = true;
            }
        }
        else
        {
            double EnergyLoss = 0.;
            while(--noSecondaries>=0)
            {
                if (secTypes[noSecondaries]!=9 && secTypes[noSecondaries]!=10)
                {
                    Particle sec = aParticle;
                    sec.Type = (ParticleType)secTypes[noSecondaries];
                    sec.Energy = NucleonEnergy * secEfrac[noSecondaries];
                    sec.fCascadeProductionTime = aParticle.Time;
                    aSecondaries.push_back(sec);
                    EnergyLoss += sec.Energy;
                }
            }
            Particle mainSec = aParticle;
            mainSec.Type = mainSecType;
            mainSec.fCascadeProductionTime = aParticle.Time;
            mainSec.Energy -= EnergyLoss;  
            aSecondaries.push_back(mainSec);
        }
        
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
