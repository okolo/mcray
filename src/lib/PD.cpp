/*
 * PhotoDisintegration.h
 */
#include "PD.h"
#include "TableFunction.h"
#include "Randomizer.h"

#include <sstream>
#include <fstream>
#include <iostream>
//#include <typeinfo>

using namespace mcray;
using namespace Utils; 
namespace Interactions {

PD::PD(mcray::BackgroundIntegral *aBackground):
        fBackground(aBackground)
{
    TablesDir = "tables/PD/";
    TableSumXS = TablesDir+"xs_pd_sum.txt";
    TableThinXS = TablesDir+"xs_pd_thin.txt";
    TableEps = TablesDir+"eps.txt";

    NuclSigma.resize(EndRealNuclei-StartRealNuclei);

    for (int i=StartRealNuclei; ParticleType(i)<EndRealNuclei; i++)
    {
        int Z = Particle::ElectricCharge(ParticleType(i));
        int N = Particle::AtomicMass(ParticleType(i))-Z;
        double M = Particle::Mass(ParticleType(i));

        NuclSigma[i-StartRealNuclei] = InitSigma(Z, N, M);
    }
}

Function* PD::InitSigma(int aPrimZ, int aPrimN, double aPrimM)                          //TODO: initialise all sigmas at once
{
    Utils::TableReader reader(TableEps, 1);
    std::vector<double>& s = reader.getColumn(0);
    std::vector<double> sigma(s.size(), 0.0);
    int len = s.size();
    int n_i = 0;
    int n_f = 0;

    std::ifstream file(TableSumXS.c_str());
    if (file.eof())
			Exception::Throw("Failed to open file " + TableSumXS);

    std::string line;
    while (std::getline(file, line)) {
		if (line[0] == '#')
			continue;
		std::istringstream iss(line);

        int Z, N;
		iss >> Z;
		iss >> N;
        if(Z != aPrimZ || N != aPrimN)
            continue;
        
        for(int i=0; i<s.size(); i++)
        {
            iss >> sigma[i];
            sigma[i] *= (units.barn*1e-3);
            s[i] = aPrimM*(aPrimM+2.*s[i]*units.MeV);
            
            if (sigma[i] != 0 && n_i == 0)
                n_i = i;
            if (sigma[i] == 0 && n_i != 0 && n_f == 0)
                n_f = i-1;
            if (sigma[i] != 0)
                n_f = 0;
            //sigma[i] *= 1e40;
            //std::cout<<s[i]<<"\t"<<sigma[i]<<std::endl;
        }
	}
    //return new LinearFunc(s, sigma, new LogScale(), new LogScale(), 0., 0.);
    //Needed to avoid:
    //>>GSL error 1 occurred in interp.c(150)
    //>>reason:interpolation error
    if (n_f == 0)
        n_f = s.size()-1;
    len = n_f-n_i+1;
    std::vector<double> sCut(len, 0.0);
    std::vector<double> sigmaCut(len, 0.0);
        
    for(int i=0; i<s.size(); i++)
    {
        if (i>=n_i && i<=n_f)
        {
            sCut[i-n_i] = s[i];
            sigmaCut[i-n_i] = sigma[i];
        }
    }
	return new LinearFunc(sCut, sigmaCut, new LogScale(), new LogScale(), 0., 0.);
}

int PD::Branching(int aPrimZ, int aPrimN, double aPrimM, double aS, Randomizer &aRandomizer) const
{
    Utils::TableReader reader(TableEps, 1);

    std::vector<double>& s = reader.getColumn(0);
    for(int i=0; i<s.size(); i++)
    {
        s[i] = aPrimM*(aPrimM+2.*s[i]*units.MeV);
    }
    
    std::vector<double> branchSigma(MaxNBranches, 0.0);
    std::vector<int> branch(MaxNBranches, 0);

    double TotalSigma = 0.;
    int counter = 0;

    std::ifstream file(TableThinXS.c_str());
    if (file.eof())
			Exception::Throw("Failed to open file " + TableThinXS);

    std::string line;

    while (std::getline(file, line)) {
	    if (line[0] == '#')
			continue;
		std::istringstream iss(line);

        int Z, N;
		iss >> Z;
		iss >> N;

        if(Z != aPrimZ || N != aPrimN)
            continue;
        
        iss >> branch[counter];
        std::vector<double> sigma(s.size(), 0.0);
        for(int i=0; i<s.size(); i++)
        {
            iss >> sigma[i];
            sigma[i] *= (units.barn*1e-3);
        }
        const Function* bsigma = new LinearFunc(s, sigma, new LogScale(), new LogScale(), 0., 0.);
        branchSigma[counter] = bsigma->f(aS);
        TotalSigma += branchSigma[counter];

        counter++;
	}


    double aRand = aRandomizer.Rand();
    int aBranch = 0;
    
    for (int i=0; i <= 17; i++)
    {
        branchSigma[i] /= TotalSigma;
        if(aRand>0 && aRand<branchSigma[i])
        {
            aBranch = branch[i];
            break;
        }
        aRand -= branchSigma[i];
    }
    file.close();
    return aBranch;
}


double PD::Rate(const mcray::Particle &aParticle) const {
    if(aParticle.Type<StartRealNuclei || aParticle.Type>EndRealNuclei)
        return 0.;
    
    if(aParticle.Type<Li6)
        return 0.;                                                                  //UNTILL proper crossections are added


    const Function* sigma = NuclSigma[aParticle.Type-StartRealNuclei];
    //std::cout<<fBackground->GetRateS(*sigma, aParticle)<<std::endl;
    return fBackground->GetRateS(*sigma, aParticle);
}

RandomInteraction *PD::Clone() const {
    return new PD(fBackground->Clone());
}

bool PD::SampleS(const mcray::Particle &aParticle, double &aS,
                                mcray::Randomizer &aRandomizer) const {
    if(aParticle.Type<StartRealNuclei || aParticle.Type>EndRealNuclei)
        return false;
    const Function* sigma = NuclSigma[aParticle.Type-StartRealNuclei];;
    double Rate = (fBackground->GetRateAndSampleS(*sigma, aParticle, aRandomizer, aS)!=0);
    return Rate;
}

void PD::SampleSecondaries(mcray::Particle &aParticle,
                                          std::vector<mcray::Particle> &aSecondaries, double aS,
                                          mcray::Randomizer &aRandomizer) const {
    ASSERT(aParticle.Type<StartRealNuclei || aParticle.Type>EndRealNuclei);
    ASSERT(aParticle.Type<Li6 || aParticle.Type>EndRealNuclei);                         //UNTILL proper crossections are added

    int aPrimZ = Particle::ElectricCharge(aParticle.Type);
    int aPrimN = Particle::AtomicMass(aParticle.Type)-aPrimZ;
    double aPrimM = Particle::Mass(aParticle.Type);
    //std::cout<<"S = "<<aS<<std::endl;
    //std::cout<<"Sigma(S) = "<<NuclSigma[aParticle.Type-StartRealNuclei]->f(aS)<<std::endl;

    int Branch = Branching(aPrimZ, aPrimN, aPrimM, aS, aRandomizer);
    //std::cout<<Branch<<std::endl;
    int nNeutron = Branch /100000 %10;
	int nProton = Branch /10000 %10;
	int nH2 = Branch /1000 %10;
	int nH3 = Branch /100 %10;
	int nHe3 = Branch /10 %10;
	int nHe4 = Branch %10;

    int dZ = nProton + nH2 + nH3 + 2 * nHe3 + 2 * nHe4;
    int dN = nNeutron + nH2 + 2 * nH3 + nHe3 + 2 * nHe4;

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

    if (hasMainSec)
    {
        Particle mainSec = aParticle;
        mainSec.Type = mainSecType;
        mainSec.fCascadeProductionTime = aParticle.Time;
        mainSec.Energy *= (mainSec.Mass()/aParticle.Mass());        //same gamma factor
        aSecondaries.push_back(mainSec);
    }
    else 
    {
        if(mainSecZ>0 || mainSecN>0)
        {
            std::cout<<std::endl<<"PD:\tinvalid channel "<<Branch<<" for "<<aParticle.ToString()<<std::endl;
            std::cout<<"\tMay be secondary particle is an unstable short-living nuclei?"<<std::endl;
            std::cout<<"\tPhotoDesintegration cancelled"<<std::endl<<std::endl;
            Particle survivedPrim = aParticle;
            survivedPrim.Ninteractions -= 1;
            survivedPrim.Weight *= 0;                               //allows to exclude all future secondaries, but not previous
            aSecondaries.push_back(survivedPrim);
            invalidBranch = true;
        }
    }
    if (!invalidBranch)
    {
        if (nNeutron)
        {
            Particle secNeutron = aParticle;
            secNeutron.Type = Neutron;
            secNeutron.Weight *= nNeutron;
            secNeutron.Energy *= (secNeutron.Mass()/aParticle.Mass());//same gamma factor
            secNeutron.fCascadeProductionTime = aParticle.Time;
            aSecondaries.push_back(secNeutron);
        }
        if (nProton)
        {
            Particle secProton = aParticle;
            secProton.Type = Proton;
            secProton.Weight *= nProton;
            secProton.Energy *= (secProton.Mass()/aParticle.Mass());//same gamma factor
            secProton.fCascadeProductionTime = aParticle.Time;
            aSecondaries.push_back(secProton);
        }
        if (nH2)
        {
            Particle secH2 = aParticle;
            secH2.Type = H2;
            secH2.Weight *= nH2;
            secH2.Energy *= (secH2.Mass()/aParticle.Mass());//same gamma factor
            secH2.fCascadeProductionTime = aParticle.Time;
            aSecondaries.push_back(secH2);
        }
        if (nH3)
        {
            Particle secH3 = aParticle;
            secH3.Type = H3;
            secH3.Weight *= nH3;
            secH3.Energy *= (secH3.Mass()/aParticle.Mass());//same gamma factor
            secH3.fCascadeProductionTime = aParticle.Time;
            aSecondaries.push_back(secH3);
        }
        if (nHe3)
        {
            Particle secHe3 = aParticle;
            secHe3.Type = He3;
            secHe3.Weight *= nHe3;
            secHe3.Energy *= (secHe3.Mass()/aParticle.Mass());//same gamma factor
            secHe3.fCascadeProductionTime = aParticle.Time;
            aSecondaries.push_back(secHe3);
        }
        if (nHe4)
        {
            Particle secHe4 = aParticle;
            secHe4.Type = He4;
            secHe4.Weight *= nHe4;
            secHe4.Energy *= (secHe4.Mass()/aParticle.Mass());//same gamma factor
            secHe4.fCascadeProductionTime = aParticle.Time;
            aSecondaries.push_back(secHe4);
        }
    }
}
}