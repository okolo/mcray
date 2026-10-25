#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"
#include "crpropa/Common.h"
#include "crpropa/Units.h"
#include "weighted_common.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
using namespace crpropa;
namespace {
using Identity=std::pair<unsigned,unsigned>;
struct Branch {unsigned code;std::vector<double> values;};
class Tables {
public:
    std::map<Identity,std::vector<double>> rates;
    std::map<Identity,std::vector<Branch>> branches;
    std::map<Identity,std::vector<std::pair<unsigned,double>>> decays;
    Tables() {
        auto read=[&](const std::string& name,auto consume) {
            std::ifstream input(getDataPath(name));if(!input) throw std::runtime_error("Missing official table: "+name);
            for(std::string line;std::getline(input,line);) {
                if(line.empty() || line[0]=='#') continue;
                std::istringstream row(line);unsigned z,n;if(!(row>>z>>n)) throw std::runtime_error("Invalid official identity");
                consume(Identity{z+n,z},row);
            }
        };
        auto grid=[](std::istringstream& row) {
            std::vector<double> result(201);for(double& v:result) if(!(row>>v) || !std::isfinite(v) || v<0.)
                throw std::runtime_error("Invalid official PD grid");return result;
        };
        read("Photodisintegration/rate_Shared.txt",[&](Identity id,auto& row){
            auto values=grid(row);for(double& v:values) v/=Mpc;rates[id]=values;
        });
        read("Photodisintegration/branching_Shared.txt",[&](Identity id,auto& row){unsigned code;row>>code;branches[id].push_back({code,grid(row)});});
        read("nuclear_decay.txt",[&](Identity id,auto& row){unsigned code;double tau;
            if(!(row>>code>>tau) || !std::isfinite(tau) || tau<=0.) throw std::runtime_error("Invalid official lifetime");
            decays[id].push_back({code,1./tau/c_light});});
    }
};
std::string name(int id) {
    if(isNucleus(id)) {
        const auto a=massNumber(id),z=chargeNumber(id);
        if(a==1) return z?"proton":"neutron";
        std::ostringstream out;out<<'A'<<std::setw(2)<<std::setfill('0')<<a<<'Z'<<std::setw(2)<<z;return out.str();
    }
    switch(id) {case 11:return "electron";case -11:return "positron";
        case 12:return "neutrinoE";case -12:return "neutrinoAE";case 22:return "photon";}
    throw std::runtime_error("Unexpected weighted secondary");
}
Identity identity(Candidate* p) {return {massNumber(p->current.getId()),chargeNumber(p->current.getId())};}
std::vector<ref_ptr<Candidate>> photoStep(ref_ptr<Candidate> initial,const Tables& tables,
    PhotoDisintegration& photo,const weighted::Guide& guide,double fraction,unsigned copies,
    std::vector<ref_ptr<Candidate>>& pending,double maxStep,unsigned long long& continuations) {
    struct Work {ref_ptr<Candidate> p;double remaining;};
    std::vector<Work> work{{initial,initial->getCurrentStep()}};
    std::vector<ref_ptr<Candidate>> completed;
    unsigned events=0;
    while(!work.empty()) {
        auto w=work.back();work.pop_back();auto p=w.p;double step=w.remaining;
        if(step<=0.) {completed.push_back(p);continue;}
        while(true) {
            if(++events>1000000) throw std::runtime_error("Weighted PD event guard exceeded");
            const auto id=identity(p);const auto found=tables.rates.find(id);
            const double lg=std::log10(p->current.getLorentzFactor());
            if(id.first<=1 || found==tables.rates.end() || lg<=6. || lg>=14.) break;
            const double rate=interpolateEquidistant(lg,6.,14.,found->second);
            if(rate<0. || !std::isfinite(rate)) throw std::runtime_error("Invalid official PD rate");
            if(rate==0.) break;
            // Same z=0 waiting law and next-step limit as pinned PhotoDisintegration::process.
            const double flight=-std::log(Random::instance().rand())/rate;
            if(step<flight) {p->limitNextStep(.1/rate);break;}
            const auto& channels=tables.branches.at(id);
            const int node=std::lround((lg-6.)/8.*200.);
            std::vector<double> raw,g;
            for(const auto& c:channels) {raw.push_back(c.values[node]);g.push_back(guide.score('P',id.first,id.second,weighted::channelCode(c.code,6)));}
            const auto physical=weighted::officialProbabilities(raw);
            const unsigned n=id.first==56?copies:1;
            continuations+=n-1;
            if(n>1) {
                // Previously emitted rewards precede this likelihood change.
                // Keep them exactly once, at their old weights, outside the clones.
                for(auto& child:p->secondaries) {child->setNextStep(maxStep);pending.push_back(child);}
                p->secondaries.clear();
            }
            for(unsigned k=0;k<n;++k) {
                auto branch=n==1?p:p->clone(false);
                const auto selected=weighted::selectDensity(physical,g,fraction,Random::instance().rand());
                branch->setWeight(weighted::reweight(p->getWeight(),selected.likelihood,n));
                photo.performInteraction(branch,channels[selected.index].code);
                if(n>1) work.push_back({branch,step-flight});
            }
            if(n>1) {p=nullptr;break;}
            step-=flight;if(step<=0.) break;
        }
        if(p) completed.push_back(p);
    }
    return completed;
}
void decayStep(Candidate* p,const Tables& tables,NuclearDecay& decay,const weighted::Guide& guide,double fraction) {
    double step=p->getCurrentStep();unsigned events=0;
    do {
        if(++events>1000000) throw std::runtime_error("Weighted decay event guard exceeded");
        const auto id=identity(p);const auto found=tables.decays.find(id);
        if(id.first<=1 || found==tables.decays.end()) return;
        double flight=std::numeric_limits<double>::max(),sum=0.;
        std::vector<double> physical,g;
        for(const auto& c:found->second) {
            const double rate=c.second/p->current.getLorentzFactor();
            // Preserve the original independent exponential race. Only the
            // winner label is resampled, independent of the minimum time.
            flight=std::min(flight,-std::log(Random::instance().rand())/rate);sum+=rate;physical.push_back(rate);
            g.push_back(guide.score('D',id.first,id.second,weighted::channelCode(c.first,5)));
        }
        if(!std::isfinite(sum) || !(sum>0.)) throw std::runtime_error("Invalid official decay rate");
        const double original=1./decay.meanFreePath(p->current.getId(),p->current.getLorentzFactor());
        if(std::abs(original-sum)>5e-13*sum) throw std::runtime_error("Weighted decay hazard differs from official module");
        if(step<flight) {p->limitNextStep(.1/sum);return;}
        for(double& v:physical) v/=sum;
        const auto selected=weighted::selectDensity(physical,g,fraction,Random::instance().rand());
        p->setWeight(weighted::reweight(p->getWeight(),selected.likelihood));
        decay.performInteraction(p,found->second[selected.index].first);step-=flight;
    } while(step>0.);
}
}
int main(int argc,char** argv) {
    try {
        std::string output,guidePath;unsigned n=100,seed=3024000000u,copies=1,a=56,z=26;bool forcePD=false;
        double fraction=0.,gamma=4279918658.3774633,length=1.,stepFraction=.03/9;
        for(int i=1;i<argc;++i) {
            if(i+1==argc) throw std::runtime_error("Missing argument");
            const std::string key=argv[i],v=argv[++i];
            if(key=="--output") output=v;else if(key=="--guide") guidePath=v;
            else if(key=="--primaries") n=std::stoul(v);else if(key=="--seed") seed=std::stoul(v);
            else if(key=="--copies") copies=std::stoul(v);else if(key=="--fraction") fraction=std::stod(v);
            else if(key=="--a") a=std::stoul(v);else if(key=="--charge") z=std::stoul(v);
            else if(key=="--gamma") gamma=std::stod(v);else if(key=="--distance-mpc") length=std::stod(v);
            else if(key=="--step-fraction") stepFraction=std::stod(v);
            else if(key=="--force-pd-adapter" && (v=="0" || v=="1")) forcePD=v=="1";
            else throw std::runtime_error("Unknown weighted reference argument: "+key);
        }
        if(!n || a<2 || a>56 || !z || z>=a || !copies || copies>64 || !std::isfinite(fraction) || fraction<0. || fraction>=1. ||
            !std::isfinite(gamma) || gamma<1. || !std::isfinite(length) || length<=0. ||
            !std::isfinite(stepFraction) || stepFraction<=0. || stepFraction>1.) throw std::runtime_error("Invalid weighted reference settings");
        const weighted::Guide guide(guidePath);const Tables tables;
        NuclearDecay decay(true,true,true);PhotoDisintegration photo(new TabularPhotonField("Shared",false),true,.1);
        Random::instance().seed(seed);
        std::ofstream out(output);if(!out) throw std::runtime_error("Cannot write reference output");
        out<<std::setprecision(17)<<"# backend=crpropa process=talys-decay primary="<<name(nucleusId(a,z))<<" seed="<<seed
           <<" primaries="<<n<<" fraction="<<fraction<<" copies="<<copies<<" L_mpc="<<length
           <<" gamma="<<gamma<<" step_fraction="<<stepFraction<<'\n'
           <<"# PD_photons=On decay_photons=On EM_transport=Off free_neutron_decay=Off nuclear_energy_debit=Off\n"
           <<"# weighted_energy=1 genealogy=primary batch history species weight weighted_energy_eV\n";
        unsigned long long continuations=0;
        for(unsigned h=0;h<n;++h) {
            std::vector<ref_ptr<Candidate>> pending{new Candidate(nucleusId(a,z),nuclearMass(a,z)*c_squared*gamma,Vector3d(0,0,0),Vector3d(1,0,0))};
            const double maxStep=length*Mpc*stepFraction;pending.front()->setNextStep(maxStep);
            unsigned long steps=0;
            while(!pending.empty()) {
                auto p=pending.back();pending.pop_back();
                if(p->current.getPosition().x>=length*Mpc) {
                    out<<0<<'\t'<<h<<'\t'<<name(p->current.getId())<<'\t'<<p->getWeight()<<'\t'<<p->getWeight()*p->current.getEnergy()/eV<<'\n';continue;
                }
                if(++steps>10000000) throw std::runtime_error("Weighted reference step guard exceeded");
                SimplePropagation propagation(0.,maxStep);
                p->limitNextStep(length*Mpc-p->current.getPosition().x);propagation.process(p);
                std::vector<ref_ptr<Candidate>> variants;
                if(fraction==0. && copies==1 && !forcePD) {
                    if(isNucleus(p->current.getId()) && massNumber(p->current.getId())>1) {photo.process(p);decay.process(p);}
                    variants.push_back(p);
                } else {
                    if(isNucleus(p->current.getId()) && massNumber(p->current.getId())>1) variants=photoStep(p,tables,photo,guide,fraction,copies,pending,maxStep,continuations);
                    else variants.push_back(p);
                    for(auto& branch:variants) if(isNucleus(branch->current.getId()) && massNumber(branch->current.getId())>1) {
                        if(fraction==0.) decay.process(branch);
                        else decayStep(branch,tables,decay,guide,fraction);
                    }
                }
                for(auto& branch:variants) {
                    for(auto& child:branch->secondaries) {child->setNextStep(maxStep);pending.push_back(child);}
                    branch->secondaries.clear();pending.push_back(branch);
                }
            }
        }
        out<<"# extra_continuations="<<continuations<<"\n# batch=0 primaries="<<n<<'\n';if(!out) throw std::runtime_error("Failed reference output write");
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
