/*
 * NuclearDecay.cpp
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


#include "NuclearDecay.h"


namespace Interactions {
    using namespace mcray;

    bool ND_init = 1;

    bool UnstableDecay(int aPrimZ, int aPrimN, Particle &aParticle, std::vector<Particle> &aSecondaries,
                                      Randomizer &aRandomizer) {
    
    double Me = Particle::Mass(Electron);
    double delta_m = Particle::Mass(Neutron)-Particle::Mass(Proton);
    double decay_k = 0.5*(delta_m-Me*Me/delta_m);
    double decay_E = sqrt(decay_k*decay_k+Me*Me);

    if (ND_init)
    {
        // ToDo: Don't read file every time
    }
    double aMass = (aPrimZ*Particle::Mass(Proton) + aPrimN*Particle::Mass(Neutron));
    double aRand = aRandomizer.Rand();
    double gamma = aParticle.Energy/aMass;
    double beta = aParticle.beta();


    int aBranch = 0;
    int afirstPhoton = 0;
    int aNPhotons = 0;
    double aPhotonEnergy = 0;

    double decayM = 0; 
    std::vector<double> branchDecayM(MaxNBranches, 0.0);
    std::vector<int> branch(MaxNBranches, 0);
    std::vector<int> NPossiblePhotons(MaxNBranches, 0);

    std::vector<double> PhotonEnergy;
    std::vector<double> PhotonProbability;

    std::ifstream file(TableLifetime.c_str());
    if (file.eof())
            Exception::Throw("Failed to open file " + TableLifetime);

    std::string line;
    int counter = 0;

    while (std::getline(file, line)) {
        if (line[0] == '#')
            continue;
        std::istringstream iss(line);

        int Z, N;
        iss >> Z;
        iss >> N;

        if(Z != aPrimZ || N != aPrimN)
            continue;
        
        double aLifetime;

        iss >> branch[counter];
        iss >> aLifetime;

        double E;
        double P;

        while(iss >> E >> P)
        {
            PhotonEnergy.push_back(E);
            PhotonProbability.push_back(P);
        }

        NPossiblePhotons[counter] = PhotonProbability.size();


        branchDecayM[counter] = aMass/aLifetime/units.second;
        decayM += branchDecayM[counter];
        counter ++;
    }

    for (int i=0; i <= MaxNBranches; i++)
    {
        if(aRand>0 && aRand<branchDecayM[i]/decayM)
        {
            aBranch = branch[i];
            aNPhotons = NPossiblePhotons[i];
            break;
        }
        afirstPhoton += NPossiblePhotons[i];
        aRand -= branchDecayM[i]/decayM;
    }

    int nBetaMinus = aBranch /10000 %10;
    int nBetaPlus = aBranch /1000 %10;
    int nAlpha = aBranch /100 %10;
    int nProton = aBranch /10 %10;
    int nNeutron = aBranch %10;

    int dZ = - nBetaMinus + nBetaPlus + 2*nAlpha + nProton;
    int dN = + nBetaMinus - nBetaPlus + 2*nAlpha + nNeutron;

    int mainSecZ = aPrimZ-dZ;
    int mainSecN = aPrimN-dN;
    double aSecMass = (mainSecZ*Particle::Mass(Proton) + mainSecN*Particle::Mass(Neutron));

    ParticleType mainSecType = Particle::FindNucleus(mainSecZ, mainSecN);
    Particle mainSec = aParticle;
    

    aRand = aRandomizer.Rand();
    for (int i=afirstPhoton; i < afirstPhoton+aNPhotons; i++)
    {
        if(aRand>0 && aRand<PhotonProbability[i])
        {
            aPhotonEnergy = PhotonEnergy[i]*1e3*units.eV;

            double cosTheta = 2 * aRandomizer.Rand() - 1;
            mainSec.Energy -= aPhotonEnergy* gamma * (1. - cosTheta);   
            //std::cout<<std::endl<<"Photon energy: "<<aPhotonEnergy* gamma * (1. - cosTheta);
            break;
        }
        aRand -= PhotonProbability[i];
    }

    if (mainSecType == ParticleTypeEOF)
    {
        UnstableDecay(mainSecZ, mainSecN, mainSec, aSecondaries, aRandomizer);
    }
    else
    {
        mainSec.Type = mainSecType;
        mainSec.fCascadeProductionTime = aParticle.Time;
    }

    double aEnergyLoss = 0;

    if (nBetaMinus)
    {
        //std::cout<<std::endl<<"beta-";
        
        double angle_e = aRandomizer.Rand()*M_PI;
        double cos_e = cos(angle_e);
        
        Particle electron = aParticle;
        electron.Type = Electron;
        electron.Weight *= nBetaMinus;
        electron.fCascadeProductionTime = aParticle.Time;
        electron.Energy = gamma*(decay_E-beta*cos_e*decay_k);
        aEnergyLoss += electron.Energy*nBetaMinus;

        //neutrino
        Particle nu = aParticle;
        nu.Type = NeutrinoAE;
        nu.Weight *= nBetaMinus;
        nu.fCascadeProductionTime = aParticle.Time;
        nu.Energy = gamma*decay_k*(1.0+beta*cos_e);//cos_nu=-cos_e
        aEnergyLoss += nu.Energy*nBetaMinus;

        aSecondaries.push_back(electron);
        aSecondaries.push_back(nu);
    }

    if (nBetaPlus)
    {
        //std::cout<<std::endl<<"beta+";

        double angle_e = aRandomizer.Rand()*M_PI;
        double cos_e = cos(angle_e);
        
        Particle positron = aParticle;
        positron.Type = Positron;
        positron.Weight *= nBetaPlus;
        positron.fCascadeProductionTime = aParticle.Time;
        positron.Energy = gamma*(decay_E-beta*cos_e*decay_k);
        aEnergyLoss += positron.Energy*nBetaPlus;

        //neutrino
        Particle nu = aParticle;
        nu.Type = NeutrinoE;
        nu.Weight *= nBetaPlus;
        nu.fCascadeProductionTime = aParticle.Time;
        nu.Energy = gamma*decay_k*(1.0+beta*cos_e);//cos_nu=-cos_e
        aEnergyLoss += nu.Energy*nBetaPlus;

        aSecondaries.push_back(positron);
        aSecondaries.push_back(nu);
    }

    if (nNeutron)
    {
        Particle secNeutron = aParticle;
        secNeutron.Type = Neutron;
        secNeutron.Weight *= nNeutron;
        secNeutron.Energy *= (secNeutron.Mass()/aMass);//same gamma factor
        secNeutron.fCascadeProductionTime = aParticle.Time;
        aSecondaries.push_back(secNeutron);
        aEnergyLoss += secNeutron.Energy*nNeutron;
    }
    if (nProton)
    {
        Particle secProton = aParticle;
        secProton.Type = Proton;
        secProton.Weight *= nProton;
        secProton.Energy *= (secProton.Mass()/aMass);//same gamma factor
        secProton.fCascadeProductionTime = aParticle.Time;
        aSecondaries.push_back(secProton);
        aEnergyLoss += secProton.Energy*nProton;
    }
    if (nAlpha)
    {
        Particle secHe4 = aParticle;
        secHe4.Type = He4;
        secHe4.Weight *= nAlpha;
        secHe4.Energy *= (secHe4.Mass()/aMass);//same gamma factor
        secHe4.fCascadeProductionTime = aParticle.Time;
        aSecondaries.push_back(secHe4);
        aEnergyLoss += secHe4.Energy*nAlpha;
    }

    mainSec.Energy -= aEnergyLoss;
    aSecondaries.push_back(mainSec);
    
    //std::cout<<std::endl<<"Decayed: "<<aParticle.ToString()<<" to "<<mainSec.ToString()<<std::endl;

    return true;
    }

    NuclearDecay::NuclearDecay()
    {
        double Me = Particle::Mass(Electron);
        double delta_m = Particle::Mass(Neutron)-Particle::Mass(Proton);
    }

    double NuclearDecay::Rate(const Particle &aParticle) const {
        if(aParticle.Type<StartRealNuclei || aParticle.Type>EndRealNuclei)
            return 0.;

        double decayM = 0; 
        std::vector<double> branchDecayM(MaxNBranches, 0.0);
        std::vector<int> branch(MaxNBranches, 0);

        int aPrimZ = Particle::ElectricCharge(aParticle.Type);
        int aPrimN = Particle::AtomicMass(aParticle.Type)-aPrimZ;
        
        std::ifstream file(TableLifetime.c_str());
        if (file.eof())
                Exception::Throw("Failed to open file " + TableLifetime);

        std::string line;
        int counter = 0;

        while (std::getline(file, line)) {
            if (line[0] == '#')
                continue;
            std::istringstream iss(line);

            int Z, N;
            iss >> Z;
            iss >> N;

            if(Z != aPrimZ || N != aPrimN)
                continue;
            
            double aLifetime;

            iss >> branch[counter];
            iss >> aLifetime;

            branchDecayM[counter] = aParticle.Mass()/aLifetime/units.second;
            decayM += branchDecayM[counter];
            counter ++;
        }
        return decayM/aParticle.Energy;
    }


    bool NuclearDecay::GetSecondaries(Particle &aParticle, std::vector<Particle> &aSecondaries,
                                      Randomizer &aRandomizer) const {
        if(aParticle.Type<StartRealNuclei || aParticle.Type>EndRealNuclei)
            return 0.;
        
        int aPrimZ = Particle::ElectricCharge(aParticle.Type);
        int aPrimN = Particle::AtomicMass(aParticle.Type)-aPrimZ;

        UnstableDecay(aPrimZ, aPrimN, aParticle, aSecondaries, aRandomizer);

        return true;
    }


    RandomInteraction *NuclearDecay::Clone() const {
        return new NuclearDecay();
    }
}
