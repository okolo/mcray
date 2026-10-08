#include "Background.h"
#include "frozen_inoue.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace mcray;
int main(int argc,char** argv) {
    if(argc<3 || argc>4) return 2;
    try {
        const std::string kind=argv[1];
        if(kind!="cmb" && kind!="cmb-inoue") return 2;
        const unsigned nodes=argc==4?std::stoul(argv[3]):4097;
        if(nodes<3 || nodes>1000000) return 2;
        const double temperature=2.7255*units.K;
        CompoundBackground photons;
        photons.AddComponent(new PlankBackground(temperature,1e-3*temperature,1e3*temperature,0.,0.02,true));
        if(kind=="cmb-inoue") photons.AddComponent(new FrozenInoueBaseline);
        std::ofstream out(argv[2]);if(!out) return 2;
        out<<std::setprecision(17)<<"# static z=0 source="<<kind<<" nominal_T_K=2.7255 K_per_MeV="<<Units::phTemperature_mult<<'\n';
        out<<"# E_eV dn_dE_cm^-3_eV^-1\n";
        for(unsigned i=0;i<nodes;++i) {
            const double energy=photons.MinE()*std::pow(photons.MaxE()/photons.MinE(),double(i)/(nodes-1));
            const double e=energy/units.eV;
            out<<e<<'\t'<<photons.n(energy,0.)*units.cm3/e<<'\n';
        }
        return out?0:1;
    } catch(Utils::Exception* error) {std::cerr<<error->Message()<<'\n';delete error;return 1;}
      catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
