// Controlled sampler substitution: all transport physics uses the current tree.
#include "ICS.h"
#include "PropagationEngine.h"
#include "Output.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <chrono>
#include <cstdlib>
using namespace Utils;
using namespace mcray;
using namespace Interactions;
bool HistoricalSampleLogDistribution(const Function&, double, double&, double&, double, double, double);
namespace {
bool historical=false, diagnostic=false;
bool rejection=false, profiling=false, cdfMode=false;
using Clock=std::chrono::steady_clock;
double invariantSeconds=0.;
unsigned long long invariantCalls=0;
int sampleCount=20000, gridCount=4096;
std::vector<double> samples;
struct LogPDF : Function {
    const Function& source;
    mutable double maximum=0.;
    explicit LogPDF(const Function& source):source(source){}
    double f(double t) const override {const double x=std::exp(t), y=x*source(x);maximum=std::max(maximum,y);return y;}
};
// Independent composite fixed Gauss-Legendre quadrature; no adaptive ODE/CDF.
struct Reference {
    const Function& f;
    double lo,step,total=0.;
    std::vector<double> cdf;
    gsl_integration_glfixed_table* table;
    Reference(const Function& f,double xmin,double xmax):f(f),lo(std::log(xmin)),step(std::log(xmax/xmin)/gridCount),cdf(gridCount+1),table(gsl_integration_glfixed_table_alloc(16)) {
        gsl_function gf;gf.params=const_cast<Function*>(&f);
        gf.function=[](double t,void* p){double x=std::exp(t);return x*(*static_cast<Function*>(p))(x);};
        for(int i=0;i<gridCount;++i) cdf[i+1]=cdf[i]+gsl_integration_glfixed(&gf,lo+i*step,lo+(i+1)*step,table);
        total=cdf.back();if(!(total>0.)) throw std::runtime_error("Empty reference kernel");
        for(auto& c:cdf)c/=total;
    }
    ~Reference(){gsl_integration_glfixed_table_free(table);}
    double pit(double x) const {
        double t=std::clamp((std::log(x)-lo)/step,0.,double(gridCount));
        int i=std::min(int(t),gridCount-1);
        gsl_function gf;gf.params=const_cast<Function*>(&f);
        gf.function=[](double v,void* p){double y=std::exp(v);return y*(*static_cast<Function*>(p))(y);};
        return cdf[i]+gsl_integration_glfixed(&gf,lo+i*step,lo+t*step,table)/total;
    }
};
class HistoryOutput : public TSmartReferencedObj<IOutput> {
public:
    // Per-primary statistics retain correlations among descendants and zero bins.
    struct History {std::array<double,32> counts{},energy{};};
    std::map<unsigned long long,History> histories;
    void AddParticle(const Particle& p) override {
        if(!p.SourceParticle)throw std::runtime_error("Missing ancestry");
        auto& h=histories[p.SourceParticle->id];
        if(p.Energy<1e12*units.eV)return;
        int bin=std::clamp(int(std::log10(p.Energy/(1e12*units.eV))*4),0,15);
        int index=(p.Type==Photon?0:16)+bin;
        h.counts[index]+=p.Weight;h.energy[index]+=p.Weight*p.Energy/units.eV;
    }
    void Flush(bool) override {}
};
}
bool SamplingProbe(const Function& f,Randomizer& rng,double& x,double& integral,double xmin,double xmax,double tol) {
    static thread_local MathUtils math;
    if(!diagnostic) {
        const auto start=profiling?Clock::now():Clock::time_point{};
        const bool ok=historical ? HistoricalSampleLogDistribution(f,rng.Rand(),x,integral,xmin,xmax,tol)
            : rejection ? math.SampleLogDistribution(f,rng,x,integral,xmin,xmax,tol)
            : cdfMode ? math.SampleLogDistributionCDF(f,rng,x,integral,xmin,xmax,tol)
            : math.SampleLogDistributionBounded(f,rng,x,integral,xmin,xmax,tol);
        if(profiling){invariantSeconds+=std::chrono::duration<double>(Clock::now()-start).count();++invariantCalls;}
        return ok;
    }
    Reference ref(f,xmin,xmax);
    LogPDF logpdf(f);
    const double qag=math.Integration_qag(logpdf,std::log(xmin),std::log(xmax),0.,tol,1000,GSL_INTEG_GAUSS61);
    const double envelope=logpdf.maximum;
    for(int i=0;i<=gridCount*4;++i)logpdf(std::log(xmin)+std::log(xmax/xmin)*i/(gridCount*4.));
    double nrError=0.,nrIntegral=0.;
    for(int i=0;i<2048;++i) {
        const double u=(i+0.5)/2048.;
        if(!HistoricalSampleLogDistribution(f,u,x,nrIntegral,xmin,xmax,tol))throw std::runtime_error("NR failed");
        nrError=std::max(nrError,std::abs(ref.pit(x)-u));
    }
    samples.clear();std::vector<double> pits;
    for(int i=0;i<sampleCount;++i) {
        if(!math.SampleLogDistributionCDF(f,rng,x,integral,xmin,xmax,tol))throw std::runtime_error("Current sampler failed");
        samples.push_back(x);pits.push_back(ref.pit(x));
    }
    std::sort(pits.begin(),pits.end());double ks=0.;
    for(int i=0;i<sampleCount;++i)ks=std::max({ks,(i+1.)/sampleCount-pits[i],pits[i]-double(i)/sampleCount});
    std::cout<<'\t'<<nrIntegral/ref.total-1<<'\t'<<qag/ref.total-1<<'\t'<<nrError<<'\t'<<ks<<'\t'<<logpdf.maximum/envelope-1;
    return true;
}
int main(int argc,char** argv) {
    try {
        gsl_set_error_handler_off();
        if(argc<2)throw std::runtime_error("Usage: em_sampling_probe kernels [N grid tolerance] | cascade old|current|rejection N seed output.tsv distance_Mpc [max_rate_variation] | interactions current|rejection N seed output.tsv ICS|PP E_EV; MCRAY_PROFILE_S=1 enables timing");
        profiling=std::getenv("MCRAY_PROFILE_S")!=nullptr;
        cosmology.Init(0.1,0.692,67.8);
        const double tolerance=(std::string(argv[1])=="kernels" && argc>4)?std::stod(argv[4]):(std::string(argv[1])=="cascade"?1e-3:1e-4);
        if(!std::isfinite(tolerance) || tolerance<=0. || tolerance>0.01)throw std::runtime_error("Invalid tolerance");
        const double temperature=2.7255*units.K;
        CompoundBackground bg;bg.AddComponent(new PlankBackground(temperature,1e-3*temperature,1e3*temperature,0.,0.1));
        PBackgroundIntegral back=new ContinuousBackgroundIntegral(bg,0.01,std::pow(10.,0.01),0.1,tolerance);
        if(std::string(argv[1])=="kernels") {
            if(argc>2)sampleCount=std::stoi(argv[2]);if(argc>3)gridCount=std::stoi(argv[3]);
            if(sampleCount<1||gridCount<16)throw std::runtime_error("Invalid diagnostic size");
            diagnostic=true;
            std::cout<<std::setprecision(12)<<"# CMB_K=2.7255 z=0 tolerance="<<tolerance<<" seed=20261012 N="<<sampleCount<<" grid="<<gridCount<<"\n"
                <<"# process E_eV NR_rate_rel GSL_rate_rel NR_max_CDF_error current_KS envelope_shortfall secondary_fraction_mean\n";
            for(bool photon:{false,true})for(double e:{1e13,1e14,1e15,1e16,1e18,1e20}) {
                Particle p(photon?Photon:Electron,0.);p.Energy=e*units.eV;
                ICS ics(back,1e12*units.eV);GammaPP pp(back);Randomizer rng(20261012);
                if((photon?pp.Rate(p):ics.Rate(p))*units.Mpc<1e-20)continue;
                std::cout<<(photon?"PP":"ICS")<<'\t'<<e;double s;
                if(!(photon?pp.SampleS(p,s,rng):ics.SampleS(p,s,rng)))throw std::runtime_error("SampleS failed");
                double fraction=0.;
                for(double value:samples){std::vector<Particle> products;if(photon)pp.SampleSecondaries(p,products,value,rng);else ics.SampleSecondaries(p,products,value,rng);fraction+=products[0].Energy/p.Energy;}
                std::cout<<'\t'<<fraction/sampleCount<<'\n'<<std::flush;
            }
        } else if(std::string(argv[1])=="interactions" && argc==8) {
            const std::string mode=argv[2], process=argv[6];
            if((mode!="current"&&mode!="rejection"&&mode!="cdf")||(process!="ICS"&&process!="PP"))throw std::runtime_error("Invalid generation case");
            rejection=mode=="rejection";cdfMode=mode=="cdf";
            const auto n=std::stoull(argv[3]);const auto seed=std::stoul(argv[4]);const double energy=std::stod(argv[7]);
            if(!n||!seed||!std::isfinite(energy)||energy<=1e12)throw std::runtime_error("Invalid generation parameters");
            const bool pair=process=="PP";
            Particle p(pair?Photon:Electron,0.);p.Energy=energy*units.eV;
            ICS ics(back,1e12*units.eV);GammaPP pp(back);Randomizer rng(seed);
            double checksum=0.;
            const auto start=Clock::now();
            for(unsigned long long i=0;i<n;++i) {
                double s;std::vector<Particle> products;
                if(!(pair?pp.SampleS(p,s,rng):ics.SampleS(p,s,rng)))throw std::runtime_error("SampleS failed");
                if(pair)pp.SampleSecondaries(p,products,s,rng);else ics.SampleSecondaries(p,products,s,rng);
                if(products.size()!=2)throw std::runtime_error("Invalid multiplicity");
                checksum+=products[0].Energy/p.Energy;
            }
            const double seconds=std::chrono::duration<double>(Clock::now()-start).count();
            std::ofstream out(argv[5]);
            out<<std::setprecision(17)<<"# transport_seconds="<<seconds<<" invariant_seconds="<<invariantSeconds<<" invariant_calls="<<invariantCalls<<" profiling="<<profiling<<"\n"
                <<"# sampler="<<mode<<" process="<<process<<" E_eV="<<energy<<" N="<<n<<" seed="<<seed<<" tolerance="<<tolerance<<" checksum="<<checksum<<"\n";
            if(!out)throw std::runtime_error("Output failed");
        } else if(std::string(argv[1])=="cascade" && (argc==7 || argc==8)) {
            const std::string mode=argv[2];if(mode!="old"&&mode!="current"&&mode!="rejection"&&mode!="cdf")throw std::runtime_error("Invalid sampler");
            rejection=mode=="rejection";cdfMode=mode=="cdf";
            historical=mode=="old";int n=std::stoi(argv[3]);unsigned long seed=std::stoul(argv[4]);
            double distance=std::stod(argv[6]);
            const double rateVariation=argc==8?std::stod(argv[7]):0.;
            if(!std::isfinite(rateVariation)||rateVariation<0.||rateVariation>0.1)throw std::runtime_error("Invalid rate variation");
            if(n<1||!std::isfinite(distance)||distance<=0||distance>1.)throw std::runtime_error("Invalid cascade size/distance");
            ParticleStack stack;Result result;CosmoTime end;end.setT(cosmology.z2t0(0.01)+distance*units.Mpc);result.SetEndTime(end);
            auto* output=new HistoryOutput;result.AddOutput(output);
            PropagationEngine engine(stack,result,seed);
            engine.AddInteraction(PRandomInteraction(new GammaPP(back)));
            engine.AddInteraction(PRandomInteraction(new ICS(back,1e12*units.eV)));
            engine.AddInteraction(PCELInteraction(new IcsCEL(back,1e12*units.eV)));
            engine.SetPropagationFilter(new EnergyBasedFilter(1e12*units.eV));
            for(int i=0;i<n;++i){Particle p(Photon,0.01);p.Energy=1e15*units.eV;stack.AddPrimary(p);}
            const auto transportStart=Clock::now();
            engine.Run(rateVariation);
            const double transportSeconds=std::chrono::duration<double>(Clock::now()-transportStart).count();
            if(engine.FailureCount())throw std::runtime_error("Discarded cascade histories");
            std::ofstream out(argv[5]);if(!out)throw std::runtime_error("Cannot write output");
            out<<std::setprecision(17)<<"# sampler="<<mode<<" N="<<n<<" seed="<<seed<<" E0_eV=1e15 Emin_eV=1e12 L_Mpc="<<distance<<" source_z=0.01 CMB_K=2.7255 no_thinning tolerance="<<tolerance<<" max_rate_variation="<<rateVariation<<"\n# bin count_sum count_sumsq energy_sum energy_sumsq\n";
            out<<"# transport_seconds="<<transportSeconds<<" invariant_seconds="<<invariantSeconds<<" invariant_calls="<<invariantCalls<<" profiling="<<profiling<<"\n";
            for(int j=0;j<32;++j){double c=0,c2=0,e=0,e2=0;for(const auto& h:output->histories){c+=h.second.counts[j];c2+=std::pow(h.second.counts[j],2);e+=h.second.energy[j];e2+=std::pow(h.second.energy[j],2);}out<<j<<'\t'<<c<<'\t'<<c2<<'\t'<<e<<'\t'<<e2<<'\n';}
            if(!out)throw std::runtime_error("Output failed");
        } else throw std::runtime_error("Invalid arguments");
        return 0;
    }catch(Exception* e){std::cerr<<e->Message()<<'\n';delete e;}
     catch(const std::exception& e){std::cerr<<e.what()<<'\n';}
     catch(const char* e){std::cerr<<e<<'\n';}
    return 1;
}
