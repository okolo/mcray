// Isolated native PD conditional-yield diagnostics; no network or transport claim.
#include "TalysPhotoDisintegration.h"
#include <iomanip>
#include <iostream>
#include <set>
using namespace mcray;
using namespace Interactions;
int main(int argc,char** argv) {
    try {
        if(argc!=6) {std::cerr<<"Usage: nuclear_photon_probe DATA A Z E_PER_A_EV Strict|CRPropa\n";return 2;}
        const auto a=std::stoul(argv[2]),z=std::stoul(argv[3]);const double e=std::stod(argv[4]);
        const std::string policy=argv[5];
        if(a<2 || a>56 || !z || z>a || !std::isfinite(e) || e<=0. || (policy!="Strict" && policy!="CRPropa"))
            Utils::Exception::Throw("Invalid photon probe parameters");
        auto data=std::make_shared<TalysData>(argv[1],true,false,true);
        Particle p(Nucleus,0.);data->SetParticle(p,a,z);p.Energy=a*e*units.eV;
        const double t=2.7255*units.K;CompoundBackground background;
        background.AddComponent(new PlankBackground(t,1e-3*t,1e3*t,0.,0.1));
        TalysPhotoDisintegration process(new ContinuousBackgroundIntegral(background,0.01,std::pow(10.,0.01),0.1,1e-4),data,true,policy=="Strict");
        std::cout<<std::setprecision(17)<<"# A="<<a<<" Z="<<z<<" E_per_A_eV="<<e<<" CMB_K=2.7255 z=0 policy="<<policy<<"\n";
        std::cout<<"# residual_A residual_Z line_energy_MeV raw_yield probability\n";
        std::set<TalysData::Identity> transitions;
        for(const auto& channel:data->Get(a,z).channels)transitions.insert(channel.residual);
        for(const auto& residual:transitions)for(const auto& line:process.PhotonYields(p,residual))
            std::cout<<residual.first<<'\t'<<residual.second<<'\t'<<line.energy/units.MeV<<'\t'<<line.rawYield<<'\t'<<line.probability<<'\n';
        return 0;
    }catch(Utils::Exception* error){std::cerr<<error->Message()<<'\n';delete error;return 1;}
     catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
