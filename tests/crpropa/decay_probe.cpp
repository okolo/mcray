#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace crpropa;
std::string Name(int id) {
    if (isNucleus(id)) {
        const auto a=massNumber(id),z=chargeNumber(id);
        if (a==1) return z ? "proton" : "neutron";
        std::ostringstream out; out<<'A'<<std::setw(2)<<std::setfill('0')<<a<<'Z'<<std::setw(2)<<z;
        return out.str();
    }
    switch(id) {
        case 11:return "electron"; case -11:return "positron";
        case 12:return "neutrinoE"; case -12:return "neutrinoAE"; case 22:return "photon";
        default:throw std::runtime_error("Unexpected secondary identity");
    }
}
double PhotoRate(unsigned a,unsigned z,double gamma) {
    const double lg=std::log10(gamma);
    if(lg<=6. || lg>=14.) return 0.;
    std::ifstream input(getDataPath("Photodisintegration/rate_Shared.txt"));
    if(!input) throw std::runtime_error("Cannot read shared photodisintegration rate table");
    for(std::string line;std::getline(input,line);) {
        if(line.empty() || line[0]=='#') continue;
        std::istringstream row(line);unsigned zz,nn;row>>zz>>nn;
        if(zz!=z || zz+nn!=a) continue;
        std::vector<double> values(201);
        for(double& value:values) if(!(row>>value)) throw std::runtime_error("Incomplete PD rate row");
        return interpolateEquidistant(lg,6.,14.,values);
    }
    return 0.; // Approved missing-PD terminal behavior; preflight is in the mcray bundle.
}
int main(int argc,char** argv) {
    try {
        std::string mode="rates",identities,output,process="decay";
        unsigned a=12,z=5,n=1000,seed=20261011,channel=0;
        double gamma=1e9,length=1e-7,stepFraction=0.1;
        for(int i=1;i<argc;++i) {
            if(i+1==argc) throw std::runtime_error("Missing argument");
            std::string key=argv[i],value=argv[++i];
            if(key=="--mode") mode=value;
            else if(key=="--process") process=value;
            else if(key=="--identities") identities=value;
            else if(key=="--output") output=value;
            else if(key=="--a") a=std::stoul(value);
            else if(key=="--charge") z=std::stoul(value);
            else if(key=="--primaries") n=std::stoul(value);
            else if(key=="--channel") channel=std::stoul(value);
            else if(key=="--seed") seed=std::stoul(value);
            else if(key=="--gamma") gamma=std::stod(value);
            else if(key=="--distance-mpc") length=std::stod(value);
            else if(key=="--step-fraction") stepFraction=std::stod(value);
            else throw std::runtime_error("Unknown option: "+key);
        }
        if(output.empty() || !n || a<2 || !z || z>a || !std::isfinite(gamma) || gamma<1. ||
           !std::isfinite(length) || length<=0. || !std::isfinite(stepFraction) || stepFraction<=0. || stepFraction>1.)
            throw std::runtime_error("Invalid probe parameters");
        NuclearDecay decay(true,false,true);
        if(process!="decay" && process!="talys" && process!="talys-decay") throw std::runtime_error("Unknown process");
        ref_ptr<PhotoDisintegration> photo;
        if(process!="decay") photo=new PhotoDisintegration(new TabularPhotonField("Shared",false),false,0.1);

        std::ofstream out(output); if(!out) throw std::runtime_error("Cannot create output");
        out<<std::setprecision(17)<<"# c_m_s="<<c_light<<" electron_mass_MeV="<<mass_electron*c_squared/MeV<<'\n';
        if(mode=="rates") {
            std::ifstream input(identities); if(!input) throw std::runtime_error("Cannot read identity inventory");
            out<<"# A Z gamma mass_MeV rate_per_s rate_per_Mpc\n";
            for(std::string line;std::getline(input,line);) {
                if(line.empty() || line[0]=='#') continue;
                std::istringstream row(line); unsigned aa,zz; if(!(row>>aa>>zz)) throw std::runtime_error("Invalid inventory");
                for(double g:{1e6,1e8,1e10,1e12}) {
                    const double mfp=decay.meanFreePath(nucleusId(aa,zz),g);
                    const double rate=mfp==std::numeric_limits<double>::max()?0.:1./mfp;
                    out<<aa<<'\t'<<zz<<'\t'<<g<<'\t'<<nuclearMass(aa,zz)*c_squared/MeV<<'\t'<<rate*c_light<<'\t'<<rate*Mpc<<'\n';
                }
            }
        } else if(mode=="events" || mode=="transport") {
            Random::instance().seed(seed);
            const int id=nucleusId(a,z);
            const double mfp=decay.meanFreePath(id,gamma);
            const double rate=(process=="talys" || mfp==std::numeric_limits<double>::max()?0.:Mpc/mfp)
                +(photo?PhotoRate(a,z,gamma):0.);
            out<<"# primary="<<Name(id)<<" gamma="<<gamma<<" L_mpc="<<length<<" event_rate_per_Mpc="<<rate
               <<" mode="<<mode<<" process="<<process<<" conditional_channel="<<channel<<" step_fraction="<<stepFraction<<" seed="<<seed<<'\n';
            out<<"# particle_rows=1\n# batch history species yield energy_eV\n";
            for(unsigned history=0;history<n;++history) {
                ref_ptr<Candidate> initial=new Candidate(id,nuclearMass(a,z)*c_squared*gamma,Vector3d(0,0,0),Vector3d(1,0,0));
                if(mode=="events") {
                    if(!channel) throw std::runtime_error("Conditional event mode requires --channel");
                    decay.performInteraction(initial,channel);
                    for(const auto& product:initial->secondaries)
                        out<<0<<'\t'<<history<<'\t'<<Name(product->current.getId())<<'\t'<<product->getWeight()<<'\t'
                           <<product->getWeight()*product->current.getEnergy()/eV<<'\n';
                    out<<0<<'\t'<<history<<'\t'<<Name(initial->current.getId())<<'\t'<<initial->getWeight()<<'\t'
                       <<initial->getWeight()*initial->current.getEnergy()/eV<<'\n';
                    continue;
                }
                std::vector<ref_ptr<Candidate>> pending{initial};
                while(!pending.empty()) {
                    auto candidate=pending.back();pending.pop_back();
                    const double maxStep=length*Mpc*stepFraction;
                    SimplePropagation propagation(0.,maxStep);
                    candidate->setNextStep(maxStep);
                    unsigned long steps=0;
                    while(candidate->current.getPosition().x<length*Mpc) {
                        if(++steps>1000000) throw std::runtime_error("Propagation step guard exceeded");
                        const double remaining=length*Mpc-candidate->current.getPosition().x;
                        candidate->limitNextStep(remaining);
                        propagation.process(candidate);
                        // Free neutron decay is excluded from this isolated nuclear fixture.
                        if(isNucleus(candidate->current.getId()) && massNumber(candidate->current.getId())>1) {
                            if(photo) photo->process(candidate);
                            if(process!="talys") decay.process(candidate);
                        }
                        pending.insert(pending.end(),candidate->secondaries.begin(),candidate->secondaries.end());
                        candidate->secondaries.clear();
                    }
                    out<<0<<'\t'<<history<<'\t'<<Name(candidate->current.getId())<<'\t'<<candidate->getWeight()<<'\t'
                       <<candidate->getWeight()*candidate->current.getEnergy()/eV<<'\n';
                }
            }
            out<<"# batch=0 primaries="<<n<<'\n';
        } else throw std::runtime_error("Unknown mode");
        if(!out) throw std::runtime_error("Output write failed");
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
