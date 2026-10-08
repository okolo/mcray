#include "NuclearDecay.h"
#include "shared_photon_table.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace mcray;
using namespace Interactions;
int main(int argc,char** argv) {
    try {
        std::string directory,output,mode="rates",channel,photons;
        unsigned a=12,z=5,n=1000,seed=20261012;
        double gamma=1e9;
        for(int i=1;i<argc;++i) {
            if(i+1==argc) throw std::runtime_error("Missing argument");
            const std::string key=argv[i],value=argv[++i];
            if(key=="--data") directory=value;
            else if(key=="--photon-table") photons=value;
            else if(key=="--output") output=value;
            else if(key=="--mode") mode=value;
            else if(key=="--a") a=std::stoul(value);
            else if(key=="--charge") z=std::stoul(value);
            else if(key=="--gamma") gamma=std::stod(value);
            else if(key=="--primaries") n=std::stoul(value);
            else if(key=="--seed") seed=std::stoul(value);
            else if(key=="--channel") channel=value;
            else throw std::runtime_error("Unknown argument: "+key);
        }
        if(directory.empty() || output.empty() || !n || !std::isfinite(gamma) || gamma<1.)
            throw std::runtime_error("Invalid probe parameters");
        auto nuclei=std::make_shared<TalysData>(directory);
        auto data=std::make_shared<NuclearDecayData>(directory,nuclei);
        NuclearDecay decay(data);
        std::ofstream out(output); if(!out) throw std::runtime_error("Cannot create output");
        out<<std::setprecision(17)<<"# electron_mass_MeV="<<Particle::Mass(Electron)/units.MeV<<'\n';
        if(mode=="rates") {
            out<<"# A Z gamma mass_MeV rate_per_s rate_per_Mpc\n";
            for(const auto& id:data->ValidateNetwork(56,26,true)) {
                if(id.first==1) continue;
                Particle p(Nucleus,0.);nuclei->SetParticle(p,id.first,id.second);
                for(double g:{1e6,1e8,1e10,1e12}) {
                    p.Energy=g*p.Mass();const double rate=decay.Rate(p);
                    out<<id.first<<'\t'<<id.second<<'\t'<<g<<'\t'<<p.Mass()/units.MeV<<'\t'
                       <<rate*units.second<<'\t'<<rate*units.Mpc<<'\n';
                }
            }
        } else if(mode=="photo-rates") {
            if(photons.empty()) throw std::runtime_error("Photo rates require --photon-table");
            SharedPhotonTable field(photons);
            auto* integral=new ContinuousBackgroundIntegral(field,0.02,std::pow(10.,0.01),0.02,1e-4);
            TalysPhotoDisintegration photo(integral,nuclei);
            Particle p(Nucleus,0.);nuclei->SetParticle(p,a,z);
            out<<"# A Z gamma channel(full,sum,or-exclusive-code) rate_per_Mpc\n";
            for(int i=0;i<=60;++i) {
                const double g=std::pow(10.,7.+i/10.);p.Energy=g*p.Mass();
                const auto rates=photo.ChannelRates(p);
                out<<a<<'\t'<<z<<'\t'<<g<<"\tfull\t"<<rates.total*units.Mpc<<'\n';
                out<<a<<'\t'<<z<<'\t'<<g<<"\tsum\t"<<rates.exclusiveSum*units.Mpc<<'\n';
                const auto& channels=nuclei->Get(a,z).channels;
                for(size_t c=0;c<channels.size();++c)
                    out<<a<<'\t'<<z<<'\t'<<g<<'\t'<<channels[c].code<<'\t'<<rates.exclusive[c]*units.Mpc<<'\n';
            }
        } else if(mode=="events") {
            const auto& channels=data->Get(a,z);
            if(channels.size()!=1 || channels[0].code!=channel)
                throw std::runtime_error("Conditional events require a separately prepared one-channel dataset");
            Particle p(Nucleus,0.);nuclei->SetParticle(p,a,z);p.Energy=gamma*p.Mass();
            Randomizer random(seed);
            out<<"# primary="<<p.Name()<<" gamma="<<gamma<<" conditional_channel="<<channel<<" seed="<<seed<<'\n';
            out<<"# particle_rows=1\n# batch history species yield energy_eV\n";
            for(unsigned history=0;history<n;++history) {
                std::vector<Particle> products;
                if(!decay.GetSecondaries(p,products,random)) throw std::runtime_error("Selected state has no decay");
                for(const auto& product:products)
                    out<<0<<'\t'<<history<<'\t'<<product.Name()<<'\t'<<product.Weight<<'\t'<<product.Weight*product.Energy/units.eV<<'\n';
            }
            out<<"# batch=0 primaries="<<n<<'\n';
        } else throw std::runtime_error("Unknown mode");
        if(!out) throw std::runtime_error("Output write failed");
    } catch(Utils::Exception* error) { std::cerr<<error->Message()<<'\n';delete error;return 1; }
      catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
