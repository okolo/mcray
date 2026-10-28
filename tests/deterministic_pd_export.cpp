// Export production coefficients only; no interaction sampling or propagation.
#include "TalysPhotoDisintegration.h"
#include "NuclearDecay.h"
#include "shared_photon_table.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
using namespace mcray;
using namespace Interactions;
int main(int argc,char** argv) {
    try {
        std::string dataPath,photonsPath,output,photons="On",gridPath,combined="Off",decayPhotons="Off";
        double energy=0.,step=.0025,eps=1e-6;
        for(int i=1;i<argc;++i) {
            if(i+1==argc) throw std::runtime_error("Missing option value");
            const std::string key=argv[i],value=argv[++i];
            if(key=="--data") dataPath=value;
            else if(key=="--photon-table") photonsPath=value;
            else if(key=="--output") output=value;
            else if(key=="--photons") photons=value;
            else if(key=="--energy-per-nucleon") energy=std::stod(value);
            else if(key=="--background-step") step=std::stod(value);
            else if(key=="--eps-rel") eps=std::stod(value);
            else if(key=="--grid") gridPath=value;
            else if(key=="--combined") combined=value;
            else if(key=="--decay-photons") decayPhotons=value;
            else throw std::runtime_error("Unknown option: "+key);
        }
        if(dataPath.empty() || photonsPath.empty() || output.empty() ||
           !std::isfinite(energy) || energy<=0 || !std::isfinite(step) || step<=0 || step>.1 ||
           !std::isfinite(eps) || eps<1e-8 || eps>.01 || (photons!="On" && photons!="Off") ||
           (combined!="On" && combined!="Off") || (decayPhotons!="On" && decayPhotons!="Off") ||
           (decayPhotons=="On" && combined!="On"))
            throw std::runtime_error("Invalid export parameters");
        auto data=std::make_shared<TalysData>(dataPath,true,false,photons=="On");
        std::shared_ptr<NuclearDecayData> decay;
        if(combined=="On")decay=std::make_shared<NuclearDecayData>(dataPath,data,decayPhotons=="On");
        const auto states=decay ? decay->ValidateNetwork(56,26,true) : data->ValidateNetwork(56,26);
        std::map<TalysData::Identity,double> energies;
        if(!gridPath.empty()) {
            std::ifstream grid(gridPath);unsigned int a,z;double e;
            if(!grid)throw std::runtime_error("Cannot read state-energy grid");
            while(grid>>a>>z>>e)if(!states.count({a,z}) || !std::isfinite(e) || e<=0. || !energies.emplace(TalysData::Identity{a,z},e).second)
                throw std::runtime_error("Invalid state-energy grid");
            if(!grid.eof() || energies.size()!=states.size())throw std::runtime_error("Incomplete state-energy grid");
        }
        SharedPhotonTable field(photonsPath);
        TalysPhotoDisintegration process(new ContinuousBackgroundIntegral(field,.02,std::pow(10.,step),.02,eps),data,photons=="On",false);
        std::ofstream out(output);if(!out)throw std::runtime_error("Cannot open export");
        out<<std::setprecision(17)<<"# PD-only fixed-E/A coefficient export; no histories\n"
           <<"# E_per_A_eV="<<energy<<" background_step="<<step<<" eps_rel="<<eps<<" photons="<<photons
           <<" combined="<<combined<<" decay_photons="<<decayPhotons<<'\n';
        for(const auto& id:states) {
            const double stateEnergy=gridPath.empty() ? energy : energies.at(id);
            Particle p(Nucleus,0.);data->SetParticle(p,id.first,id.second);p.Energy=id.first*stateEnergy*units.eV;
            const auto rates=process.ChannelRates(p);
            out<<"S "<<id.first<<' '<<id.second<<' '<<p.Mass()/units.MeV<<' '<<rates.total*units.Mpc;
            if(!gridPath.empty())out<<' '<<stateEnergy;
            out<<'\n';
            if(id.first==1)continue;
            if(decay) {
                NuclearDecay decayProcess(decay,decayPhotons=="On");
                out<<"T "<<id.first<<' '<<id.second<<' '<<decayProcess.Rate(p)*units.Mpc<<'\n';
                for(const auto& channel:decay->Get(id.first,id.second)) {
                    out<<"D "<<id.first<<' '<<id.second<<' '<<channel.code<<' '<<channel.restRate*p.Mass()/p.Energy*units.Mpc<<'\n';
                    if(decayPhotons=="On")for(const auto& line:channel.photons)
                        out<<"L "<<id.first<<' '<<id.second<<' '<<channel.code<<' '<<line.energy/units.eV<<' '<<line.probability<<'\n';
                }
            }
            const auto& isotope=data->Get(id.first,id.second);std::set<TalysData::Identity> residuals;
            for(size_t i=0;i<isotope.channels.size();++i) {
                out<<"C "<<id.first<<' '<<id.second<<' '<<isotope.channels[i].code<<' '<<rates.exclusive[i]*units.Mpc<<'\n';
                residuals.insert(isotope.channels[i].residual);
            }
            if(photons=="On")for(const auto& residual:residuals)
                for(const auto& line:process.PhotonYields(p,residual))
                    out<<"P "<<id.first<<' '<<id.second<<' '<<residual.first<<' '<<residual.second<<' '
                       <<line.energy/units.eV<<' '<<line.probability<<'\n';
        }
        out<<"# complete states="<<states.size()<<'\n';
        if(!out)throw std::runtime_error("Coefficient export write failed");
    }catch(Utils::Exception* e){std::cerr<<e->Message()<<'\n';delete e;return 1;}
     catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
