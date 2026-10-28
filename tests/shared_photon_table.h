#pragma once
#include "Background.h"
#include "TableFunction.h"
#include <fstream>
#include <memory>
#include <sstream>

// Static test-only physical density: E/eV, dn/dE in cm^-3 eV^-1.
// Production photon backgrounds and constants are deliberately unchanged.
class SharedPhotonTable : public mcray::IBackground {
public:
    explicit SharedPhotonTable(const std::string& file) {
        std::ifstream input(file);
        if(!input) Utils::Exception::Throw("Cannot read shared photon table: "+file);
        std::vector<double> energy,density;
        for(std::string line;std::getline(input,line);) {
            if(line.empty() || line[0]=='#') continue;
            std::istringstream row(line);double e,n;std::string extra;
            if(!(row>>e>>n) || row>>extra || !std::isfinite(e) || !std::isfinite(n) || e<=0. || n<0. ||
               (!energy.empty() && e<=energy.back())) Utils::Exception::Throw("Invalid shared photon table row");
            energy.push_back(e);density.push_back(n);
        }
        if(energy.size()<2) Utils::Exception::Throw("Shared photon table needs at least two nodes");
        fMin=energy.front()*mcray::units.eV;fMax=energy.back()*mcray::units.eV;
        fDensity=std::make_unique<Utils::LinearFunc>(energy,density,new Utils::LogScale,new Utils::LogScale,0.,0.,true);
    }
    double n(double energy,double) const override {
        const double e=energy/mcray::units.eV;
        return e*fDensity->f(e)/mcray::units.cm3;
    }
    double MinE() const override {return fMin;}
    double MaxE() const override {return fMax;}
    double MinZ() const override {return 0.;}
    double MaxZ() const override {return 0.02;}
    std::string Name() const override {return "SharedStaticPhotonTable";}
private:
    double fMin=0.,fMax=0.;
    std::unique_ptr<Utils::LinearFunc> fDensity;
};
