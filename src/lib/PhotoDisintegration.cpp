/*
 * PhotoDisintegration.cpp
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


#include "PhotoDisintegration.h"
#include "ParticleStack.h"
#include "TableBackgrounds.h"
#include "PropagationEngine.h"

Interactions::PhotoDisintegration::PhotoDisintegration(mcray::BackgroundIntegral *aBackground, bool aSampleK):
        fBackground(aBackground),
        fSampleK(aSampleK)
{
}

mcray::RandomInteraction *Interactions::PhotoDisintegration::Clone() const {
    return new PhotoDisintegration(fBackground);
}

Interactions::PhotoDisintegration::~PhotoDisintegration() {

}

double Interactions::PhotoDisintegration::Rate(const mcray::Particle &aParticle) const {
    if(aParticle.Type<A02 || aParticle.Type>=EndNuclei)
        return 0.;
    const CPDMLine* canals = fMap.getOutcome(CNucleus::getA(aParticle.Type));
    int noOfCanals = canals?canals->size():0;
    if(!noOfCanals)
        return 0.;
    if(fSampleK){
        SigmaK sigma(canals);
        return fBackground->GetRateK(sigma, aParticle);
    }else{
        SigmaS sigma(canals);
        return fBackground->GetRateS(sigma, aParticle);
    }
}

int Interactions::PhotoDisintegration::UnitTest()
{
    UnitTest(A04, true, true, false, 5e18*units.eV, 1e20*units.eV, 1, 10000);
    return 0;
}

void Interactions::PhotoDisintegration::UnitTest(ParticleType aPrimary, bool aPrintRates, bool aPropagate, bool aSplit, double Emin, double Emax, double Zmax, int nParticles)
{
    //double alphaThinning = 0;//alpha = 1 conserves number of particles on the average; alpha = 0 disables thinning
    ParticleStack particles;
    Result result(Emin/CNucleus::getA(aPrimary));
    if(aPropagate)
    {
        std::string out = Particle::Name(aPrimary);
        out += aSplit?"_split" : "";
        out += "_Emin";
        out += ToString(Emin);
        out += "_Emax";
        out += ToString(Emax);
        out += "_Zmax";
        out += ToString(Zmax);
        //result.AddOutput(new SpectrumOutput("out_" + out, Emin, pow(10, 0.05)));
        result.AddOutput(new RawOutput(out, false, true));
        debug.SetOutputFile("debug_" + out);
    }
    if(!cosmology.IsInitialized())
        cosmology.Init();

    CompoundBackground backgr;
    double cmbTemp = 2.73/Units::phTemperature_mult/units.Eunit;

    backgr.AddComponent(new PlankBackground(cmbTemp, 1e-3*cmbTemp, 1e3*cmbTemp));
    //backgr.AddComponent(new PlankBackground(cmbTemp, 6.2e-10/units.Eunit, 6.4e-10/units.Eunit));//only central value

    double stepZ = 0.05;
    double logStepK = pow(10,0.05);
    double epsRel = 1e-3;

    Backgrounds::Kneiske0309EBL* kn0309 = new Backgrounds::Kneiske0309EBL();
    PBackgroundIntegral backgrIebl;
    PRandomInteraction interactionEbl;
    if(!aSplit)
    {
        backgr.AddComponent(kn0309);
    }
    else
    {
        backgrIebl = new ContinuousBackgroundIntegral(*kn0309, stepZ, logStepK, Zmax, epsRel);
        interactionEbl = new Interactions::PhotoDisintegration(backgrIebl);
    }

    PBackgroundIntegral backgrI = new ContinuousBackgroundIntegral(backgr, stepZ, logStepK, Zmax, epsRel);
    PRandomInteraction interaction = new Interactions::PhotoDisintegration(backgrI);

    PCELInteraction rs = new RedShift();
    double E = Emin;
    if(aPrintRates)
    {
        std::ofstream ratesOut;
        std::string ratesFile = ToString(Particle::Name(aPrimary)) + "_rates";
        ratesFile += aSplit ? "_split" : "";
        ratesOut.open(ratesFile.c_str(),std::ios::out);

        ratesOut << "#E\tRedshiftRate\tRate\t\t[E]=1eV [Rate]=1/Mpc\n";

        for(double E=Emin; E<Emax; E*=logStepK)
        {
            Particle particle(aPrimary, 0.);
            particle.Energy = E / units.Eunit;//MeV
            double rate = interaction->Rate(particle);
            double redshiftRate = rs->Rate(particle);
            if(aSplit)
            {
                rate += interactionEbl->Rate(particle);
            }
            double mult = 1./units.Lunit*Units::Mpc_in_cm;
            ratesOut << E*1e6 << "\t" << redshiftRate*mult << "\t" << rate * mult << std::endl;
        }
    }
    if(!aPropagate)
        return;

    PropagationEngine pe(particles, result);

    pe.AddInteraction(interaction);

    if(aSplit)
    {
        pe.AddInteraction(interactionEbl);
    }

    E=Emax;
    Particle primary(aPrimary, Zmax);
    primary.Energy = E;//MeV
    //primary.Weight = exp(-Emax/E);
    //if(primary.Weight>0)
    for(int i=0; i < nParticles; i++)
        particles.AddPrimary(primary);
    //}

    pe.RunMultithreadReleaseOnly();
    if(aSplit)
        delete kn0309;
}


bool Interactions::PhotoDisintegration::SampleS(const mcray::Particle &aParticle, double &aS,
                                                mcray::Randomizer &aRandomizer) const {
    ASSERT(aParticle.Type>=A02 && aParticle.Type<EndNuclei);
    const CPDMLine* canals = fMap.getOutcome(CNucleus::getA(aParticle.Type));
    int noOfCanals = canals?canals->size():0;
    ASSERT(noOfCanals>0);
    SigmaS sigma(canals);
    return (fBackground->GetRateAndSampleS(sigma, aParticle, aRandomizer, aS)!=0);
}

bool Interactions::PhotoDisintegration::SampleK(const mcray::Particle &aParticle, double &aK,
                                                mcray::Randomizer &aRandomizer) const {
    ASSERT(aParticle.Type>=A02 && aParticle.Type<EndNuclei);
    const CPDMLine* canals = fMap.getOutcome(CNucleus::getA(aParticle.Type));
    int noOfCanals = canals?canals->size():0;
    ASSERT(noOfCanals>0);
    SigmaK sigma(canals);
    return (fBackground->GetRateAndSampleK(sigma, aParticle, aRandomizer, aK)!=0);
}

bool Interactions::PhotoDisintegration::GetSecondaries(Particle& aParticle, std::vector<Particle>& aSecondaries, Randomizer& aRandomizer) const
{
    if(fSampleK){
        double k = 0.;
        if(!SampleK(aParticle, k, aRandomizer))
            return false;
        SampleSecondariesWithK(aParticle, aSecondaries, k, aRandomizer);
    }
    else{
        double s = 0.;
        if(!SampleS(aParticle, s, aRandomizer))
            return false;
        SampleSecondariesWithS(aParticle, aSecondaries, s, aRandomizer);
    }
    return true;
}

void Interactions::PhotoDisintegration::SampleSecondariesWithS(mcray::Particle &aParticle,
                                                               std::vector<mcray::Particle> &aSecondaries, double aS,
                                                               mcray::Randomizer &aRandomizer) const {
    ASSERT(aParticle.Type>=A02 && aParticle.Type<EndNuclei);
    double m=aParticle.Mass();
    double k = 0.5*(aS/m-m);
    SampleSecondariesWithK(aParticle, aSecondaries, k, aRandomizer);
}

void Interactions::PhotoDisintegration::SampleSecondariesWithK(mcray::Particle &aParticle,
                                                               std::vector<mcray::Particle> &aSecondaries, double k,
                                                               mcray::Randomizer &aRandomizer) const {
    ASSERT(aParticle.Type>=A02 && aParticle.Type<EndNuclei);
    const CPDMLine* canals = fMap.getOutcome(CNucleus::getA(aParticle.Type));
    int noOfCanals = canals?canals->size():0;
    ASSERT(noOfCanals>0);

    //select canal
    double totSigma = 0.;
    std::vector<double> sigma(canals->size(),0.);
    size_t iCanal=0;
    for (std::vector<CPhotoDisintegrationCanal *>::const_iterator it = canals->begin();
         it != canals->end(); it++)
    {
        Function* curSigma = (*it)->Sigma();
        double kMinCur = curSigma->Xmin();
        if (kMinCur < k)
            totSigma += curSigma->f(k);
        sigma[iCanal]=totSigma;
        iCanal++;
    }
    double randSigma = aRandomizer.Rand()*totSigma;
    CPhotoDisintegrationCanal* selectedCanal=0;
    for (iCanal=0; iCanal<canals->size(); iCanal++)
    {
        if(sigma[iCanal]>=randSigma) {
            selectedCanal = canals->at(iCanal);
            break;
        }
    }

    double secE=0;// energy conservation test (debug only)

    const TIsotope &isotope = selectedCanal->isotope();
    int maxDeltaA = selectedCanal->getMaxDeltaA();
    int A = isotope.A;
    ASSERT(CNucleus::getA(aParticle.Type)==A);
    ASSERT(maxDeltaA < A);
    double autoN = 0.;
    double autoP = 0.;
    for (int deltaA = 1; deltaA <= maxDeltaA; deltaA++) {
        double rate = selectedCanal->getRate(deltaA);
        ParticleType prodType=(ParticleType)(aParticle.Type-deltaA);//in case A=2,3 secondary canal for proton (H1) may be duplicated, but it is ok, since summary rate is correct
        int deltaZ = isotope.Z - CNucleus::getZ(prodType);
        if (deltaZ < 0) {//processes with neutron emission followed by beta decay
            deltaZ = 0;
        }
        int deltaN = deltaA - deltaZ;
        autoN += deltaN * rate;
        autoP += deltaZ * rate;
        if (rate > 0.) {
            Particle product = aParticle;
            product.Type = prodType;
            product.Weight *= rate;
            product.Energy *= (product.Mass()/aParticle.Mass());//same gamma factor
            aSecondaries.push_back(product);
            secE += (product.Energy * product.Weight);// energy conservation test (debug only)
        }
    }
    double rateP = selectedCanal->getRateP();
    double rateN = selectedCanal->getRateN();
    if (rateP == CPhotoDisintegrationCanal::AutoRate) {
        rateP = autoP;
    };
    if (rateN == CPhotoDisintegrationCanal::AutoRate) {
        rateN = autoN;
    };
    if (rateN > 0.) {
        Particle product = aParticle;
        product.Type = Neutron;
        product.Weight *= rateN;
        product.Energy *= (product.Mass()/aParticle.Mass());//same gamma factor
        aSecondaries.push_back(product);
        secE += (product.Energy * product.Weight);// energy conservation test (debug only)
    }
    if (rateP > 0.) {
        Particle product = aParticle;
        product.Type = Proton;
        product.Weight *= rateP;
        product.Energy *= (product.Mass()/aParticle.Mass());//same gamma factor
        aSecondaries.push_back(product);
        secE += (product.Energy * product.Weight);
    }
    double primE = aParticle.Energy * aParticle.Weight;// energy conservation test (debug only)
    ASSERT(MathUtils::RelDifference(secE,primE)<1e-3);// verify energy conservation
}

Interactions::PhotoDisintegration::SigmaS::SigmaS(const CPDMLine* canals) :
fCanals(canals),
fXmin(1e300),
fXmax(0),
fM(-1.)
{
    if(fCanals->size()) {
        double kMin = 1e300;
        double kMax = 0;
        fM = Particle::Mass(fCanals->at(0)->getParticle());
        for (std::vector<CPhotoDisintegrationCanal *>::const_iterator it = fCanals->begin();
             it != fCanals->end(); it++)
        {
            double kMinCur = (*it)->Sigma()->Xmin();
            if (kMinCur < kMin)
                kMin = kMinCur;
            double kMaxCur = (*it)->Sigma()->Xmax();
            if (kMaxCur > kMax)
                kMax = kMaxCur;
        }
        fXmin = fM * (fM + 2. * kMin);
        fXmax = fM * (fM + 2. * kMax);
    }
}

double Interactions::PhotoDisintegration::SigmaS::f(double s) const {
    if(s<=fXmin)
        return 0.;
    double k = 0.5*(s/fM-fM);
    double sigma = 0.;
    for (std::vector<CPhotoDisintegrationCanal *>::const_iterator it = fCanals->begin();
         it != fCanals->end(); it++)
    {
        Function* curSigma = (*it)->Sigma();
        if (curSigma->Xmin() <= k && curSigma->Xmax() >= k) {
            sigma += curSigma->f(k);
        }
    }
    return sigma;
}


Interactions::PhotoDisintegration::SigmaK::SigmaK(const CPDMLine* canals) :
        fCanals(canals),
        fXmin(1e300),
        fXmax(0)
{
    if(fCanals->size()) {
        double kMin = 1e300;
        double kMax = 0;
        for (std::vector<CPhotoDisintegrationCanal *>::const_iterator it = fCanals->begin();
             it != fCanals->end(); it++)
        {
            double kMinCur = (*it)->Sigma()->Xmin();
            if (kMinCur < kMin)
                kMin = kMinCur;
            double kMaxCur = (*it)->Sigma()->Xmax();
            if (kMaxCur > kMax)
                kMax = kMaxCur;
        }
        fXmin = kMin;
        fXmax = kMax;
        ASSERT(kMin>0);
    }
}

double Interactions::PhotoDisintegration::SigmaK::f(double k) const {
    if(k<=fXmin || k>fXmax)
        return 0.;

    double sigma = 0.;
    for (std::vector<CPhotoDisintegrationCanal *>::const_iterator it = fCanals->begin();
         it != fCanals->end(); it++)
    {
        Function* curSigma = (*it)->Sigma();
        if (curSigma->Xmin() <= k && curSigma->Xmax() >= k) {
            sigma += curSigma->f(k);
        }
    }
    return sigma;
}