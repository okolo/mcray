#include "NuclearDecay.h"
#include "PropagationEngine.h"
#include "Output.h"
#include "shared_photon_table.h"
#include "crpropa/weighted_common.h"
#include <fstream>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace mcray;
using namespace Interactions;
namespace {
class WeightedPhoto : public RandomInteraction {
public:
    WeightedPhoto(PBackgroundIntegral background,std::shared_ptr<const TalysData> data,
                  std::shared_ptr<const weighted::Guide> guide,double fraction,unsigned copies,
                  std::shared_ptr<std::atomic<unsigned long long>> continuations,bool forcePD)
        : background(background),data(data),guide(guide),fraction(fraction),copies(copies),continuations(continuations),forcePD(forcePD),photo(background,data,true,false) {}
    double Rate(const Particle& p) const override {return photo.Rate(p);}
    RandomInteraction* Clone() const override {return new WeightedPhoto(background->Clone(),data,guide,fraction,copies,continuations,forcePD);}
    bool GetSecondaries(Particle& p,std::vector<Particle>& products,Randomizer& rng) const override {
        if (fraction==0. && copies==1 && !forcePD) return photo.GetSecondaries(p,products,rng);
        const auto rates=photo.ChannelRates(p);
        if (rates.total==0.) return false;
        const auto& channels=data->Get(p.MassNumber(),p.AtomicNumber()).channels;
        std::vector<double> physical,importance;
        for (size_t i=0;i<channels.size();++i) {
            physical.push_back(rates.exclusive[i]/rates.exclusiveSum);
            importance.push_back(guide->score('P',p.MassNumber(),p.AtomicNumber(),channels[i].code));
        }
        // A decreases strictly in PD: each initial A56 lineage splits at most once.
        const unsigned n=p.MassNumber()==56 ? copies : 1;
        *continuations+=n-1;
        for (unsigned k=0;k<n;++k) {
            const auto selected=weighted::selectDensity(physical,importance,fraction,rng.Rand());
            Particle branch=p;branch.Weight=weighted::reweight(p.Weight,selected.likelihood,n);
            photo.GetChannelSecondaries(branch,selected.index,products,rng);
        }
        return true;
    }
private:
    PBackgroundIntegral background;
    std::shared_ptr<const TalysData> data;
    std::shared_ptr<const weighted::Guide> guide;
    double fraction;unsigned copies;
    std::shared_ptr<std::atomic<unsigned long long>> continuations;
    bool forcePD;
    TalysPhotoDisintegration photo;
};
class WeightedDecay : public RandomInteraction {
public:
    WeightedDecay(std::shared_ptr<const NuclearDecayData> data,std::shared_ptr<const weighted::Guide> guide,double fraction)
        : data(data),guide(guide),fraction(fraction),decay(data,true) {}
    double Rate(const Particle& p) const override {return decay.Rate(p);}
    RandomInteraction* Clone() const override {return new WeightedDecay(data,guide,fraction);}
    bool GetSecondaries(Particle& p,std::vector<Particle>& products,Randomizer& rng) const override {
        if (fraction==0.) return decay.GetSecondaries(p,products,rng);
        if (decay.Rate(p)==0.) return false;
        const auto& channels=data->Get(p.MassNumber(),p.AtomicNumber());
        double sum=0.;for (const auto& c:channels) sum+=c.restRate;
        std::vector<double> physical,importance;
        for (const auto& c:channels) {
            physical.push_back(c.restRate/sum);
            importance.push_back(guide->score('D',p.MassNumber(),p.AtomicNumber(),c.code));
        }
        const auto selected=weighted::selectDensity(physical,importance,fraction,rng.Rand());
        Particle branch=p;branch.Weight=weighted::reweight(p.Weight,selected.likelihood);
        return decay.GetChannelSecondaries(branch,channels[selected.index].code,products,rng);
    }
private:
    std::shared_ptr<const NuclearDecayData> data;
    std::shared_ptr<const weighted::Guide> guide;
    double fraction;
    NuclearDecay decay;
};
class Output : public Utils::TSmartReferencedObj<IOutput> {
public:
    struct Row {std::string species;double weight,energy;};
    std::map<unsigned long long,std::vector<Row>> histories;
    void AddParticle(const Particle& p) override {
        if (!p.SourceParticle || !std::isfinite(p.Weight) || !(p.Weight>0.))
            throw std::runtime_error("Invalid genealogical output");
        histories[p.SourceParticle->id].push_back({p.Name(),p.Weight,p.Weight*p.Energy/units.eV});
    }
    void Flush(bool) override {}
};
}
int main(int argc,char** argv) {
    try {
        std::string dataPath,photons,guidePath,output;
        unsigned primaries=100,seed=3024000000u,copies=1,a=56,z=26;bool forcePD=false;
        double fraction=0.,step=.0025,eps=1e-6,epa=3.9810717055349857e18,length=1.;
        for (int i=1;i<argc;++i) {
            if (i+1==argc) throw std::runtime_error("Missing argument");
            const std::string key=argv[i],value=argv[++i];
            if(key=="--data") dataPath=value;else if(key=="--photon-table") photons=value;
            else if(key=="--guide") guidePath=value;else if(key=="--output") output=value;
            else if(key=="--primaries") primaries=std::stoul(value);else if(key=="--seed") seed=std::stoul(value);
            else if(key=="--copies") copies=std::stoul(value);else if(key=="--fraction") fraction=std::stod(value);
            else if(key=="--a") a=std::stoul(value);else if(key=="--charge") z=std::stoul(value);
            else if(key=="--background-step") step=std::stod(value);else if(key=="--eps") eps=std::stod(value);
            else if(key=="--energy-per-nucleon") epa=std::stod(value);else if(key=="--distance-mpc") length=std::stod(value);
            else if(key=="--force-pd-adapter" && (value=="0" || value=="1")) forcePD=value=="1";
            else throw std::runtime_error("Unknown weighted scenario option: "+key);
        }
        if (!primaries || a<2 || a>56 || !z || z>=a || !copies || copies>64 || !std::isfinite(fraction) || fraction<0. || fraction>=1. ||
            !std::isfinite(step) || step<=0. || step>.1 || !std::isfinite(eps) || eps<1e-8 || eps>.01 ||
            !std::isfinite(epa) || epa<=0. || !std::isfinite(length) || length<=0.)
            throw std::runtime_error("Invalid weighted scenario settings");
        auto guide=std::make_shared<weighted::Guide>(guidePath);
        auto nuclei=std::make_shared<TalysData>(dataPath,true,false,true);
        auto decays=std::make_shared<NuclearDecayData>(dataPath,nuclei,true);
        decays->ValidateNetwork(a,z,true);
        cosmology.Init(.1,.692,67.8);
        CompoundBackground field;field.AddComponent(new SharedPhotonTable(photons));
        PBackgroundIntegral background=new ContinuousBackgroundIntegral(field,.02,std::pow(10.,step),.02,eps);
        Particle primary(Nucleus,.01);nuclei->SetParticle(primary,a,z);primary.Energy=a*epa*units.eV;
        if (primary.Time.t()+length*units.Mpc>0.) throw std::runtime_error("Path exceeds fixed reference interval");
        ParticleStack stack;Result result;CosmoTime end;end.setT(primary.Time.t()+length*units.Mpc);result.SetEndTime(end);
        auto* records=new Output;result.AddOutput(records);
        PropagationEngine engine(stack,result,seed);
        auto continuations=std::make_shared<std::atomic<unsigned long long>>(0);
        engine.AddInteraction(new WeightedPhoto(background,nuclei,guide,fraction,copies,continuations,forcePD));
        engine.AddInteraction(new WeightedDecay(decays,guide,fraction));
        for(unsigned i=0;i<primaries;++i) {
            Particle p(Nucleus,.01);nuclei->SetParticle(p,a,z);p.Energy=primary.Energy;stack.AddPrimary(p);
        }
        engine.Run(0.);
        if (engine.FailureCount() || records->histories.size()!=primaries)
            throw std::runtime_error("Discarded or incomplete weighted histories");
        std::ofstream out(output);if(!out) throw std::runtime_error("Cannot create weighted output");
        out<<std::setprecision(17)<<"# backend=mcray process=talys-decay primary="<<primary.Name()<<" seed="<<seed
           <<" primaries="<<primaries<<" fraction="<<fraction<<" copies="<<copies<<" L_mpc="<<length
           <<" E_per_nucleon_eV="<<epa<<" background_step="<<step<<" eps="<<eps<<'\n'
           <<"# PD_photons=On PD_photon_yields=CRPropa decay_photons=On EM_transport=Off free_neutron_decay=Off nuclear_energy_debit=Off\n"
           <<"# weighted_energy=1 genealogy=primary batch history species weight weighted_energy_eV\n";
        unsigned history=0;
        for(const auto& h:records->histories) {
            for(const auto& r:h.second) out<<0<<'\t'<<history<<'\t'<<r.species<<'\t'<<r.weight<<'\t'<<r.energy<<'\n';
            ++history;
        }
        out<<"# extra_continuations="<<continuations->load()<<"\n# batch=0 primaries="<<primaries<<'\n';
        if (!out) throw std::runtime_error("Failed weighted output write");
    } catch(Utils::Exception* error) {std::cerr<<error->Message()<<'\n';delete error;return 1;}
      catch(const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}
