// Opt-in single-interaction spectra. Production files are only copied/instrumented.
#include "ICS.h"
#include "secondary_oracle.h"
#include <gsl/gsl_errno.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cstring>
using namespace Utils;
using namespace mcray;
using namespace Interactions;
namespace secondary {Audit audit;}
namespace {
using secondary::Real;
bool pairProcess=false,oracleMode=false,oldMode=false;
int cells=512,order=16;
double primaryEnergy=0,threshold=0,lastS=0,sMin=0,sMax=0;
std::string destination;
std::vector<secondary::Region> observables;
MathUtils currentMath;
void oracle(const Function& f,double lo,double hi) {
    secondary::Quadrature q(order);
    std::vector<std::array<Real,3>> sums(observables.size());
    const Real m2=std::pow(Real(Particle::Mass(Electron)),2),a=std::log(Real(lo)),b=std::log(Real(hi));
    Real total=0;
    // Split at every moving kinematic boundary as well as a uniform log(s) grid.
    std::vector<Real> knots;
    for(int i=0;i<=cells;++i)knots.push_back(a+(b-a)*i/cells);
    for(const auto& o:observables)for(Real x:{o.lo,o.hi,1-o.lo,1-o.hi})if(x>0 && x<1) {
        const Real s=pairProcess?1/(x*(1-x)):1/x;
        const Real t=std::log(s*m2);if(t>a && t<b)knots.push_back(t);
    }
    std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
    for(size_t k=1;k<knots.size();++k)for(size_t i=0;i<q.table->n;++i) {
        double t,w;gsl_integration_glfixed_point(double(knots[k-1]),double(knots[k]),i,&t,&w,q.table);
        const Real s=std::exp(Real(t)),weight=Real(w)*s*f(double(s));
        if(!std::isfinite(weight)||weight<0)throw std::runtime_error("Invalid production s kernel");
        if(weight==0)continue;
        const auto conditional=secondary::conditional(pairProcess,s/m2,threshold,observables,q);
        total+=weight;
        for(size_t j=0;j<sums.size();++j)for(int c=0;c<3;++c)sums[j][c]+=weight*conditional[j][c];
    }
    if(!(total>0))throw std::runtime_error("Empty reference s kernel");
    const auto possible=secondary::conditional(pairProcess,Real(hi)/m2,threshold,observables,q,true);
    std::ofstream out(destination);out<<std::setprecision(21);
    out<<"# process="<<(pairProcess?"PP":"ICS")<<" E_eV="<<primaryEnergy<<" Emin_eV="<<threshold*primaryEnergy<<" s_min_over_m2="<<Real(lo)/m2<<" s_max_over_m2="<<Real(hi)/m2<<" cells="<<cells<<" order="<<order<<" kernel_integral="<<total<<'\n';
    out<<"# id view kind lo hi p0 p1 p2 possible0 possible1 possible2\n";
    for(size_t j=0;j<sums.size();++j) {
        const auto& o=observables[j];out<<o.id<<'\t'<<o.view<<'\t'<<o.kind<<'\t'<<o.lo<<'\t'<<o.hi;
        for(int c=0;c<3;++c)out<<'\t'<<sums[j][c]/total;
        for(int c=0;c<3;++c)out<<'\t'<<(possible[j][c]>0);
        out<<'\n';
    }
    if(!out)throw std::runtime_error("Cannot save independent reference");
}
void binaryDouble(std::ostream& out,double x){out.write(reinterpret_cast<const char*>(&x),sizeof(x));}
}
bool HistoricalSampleLogDistribution(const Function&,double,double&,double&,double,double,double);
bool SamplingProbe(const Function& f,Randomizer& rng,double& x,double& integral,double lo,double hi,double tol) {
    sMin=lo;sMax=hi;
    if(oracleMode){oracle(f,lo,hi);x=std::sqrt(lo*hi);integral=1.;return true;}
    const bool ok=oldMode?HistoricalSampleLogDistribution(f,rng.Rand(),x,integral,lo,hi,tol):currentMath.SampleLogDistribution(f,rng,x,integral,lo,hi,tol);
    lastS=x;
    if(!ok||!std::isfinite(integral)||integral<=0||!std::isfinite(x)||x<lo||x>hi)throw std::runtime_error("Invalid sampled invariant");
    return true;
}
int main(int argc,char** argv) {
    try {
        if(argc!=10 && argc!=12)throw std::runtime_error("Usage: em_secondary_probe oracle PROCESS E_EV EMIN_EV TOL EDGES CELLS ORDER OUTPUT | sample PROCESS E_EV EMIN_EV TOL EDGES old|current N SEED EVENTS STATS");
        gsl_set_error_handler_off();
        oracleMode=std::string(argv[1])=="oracle";
        if(!oracleMode && std::string(argv[1])!="sample")throw std::runtime_error("Invalid mode");
        if((oracleMode&&argc!=10)||(!oracleMode&&argc!=12))throw std::runtime_error("Invalid argument count");
        const std::string process=argv[2];if(process!="PP"&&process!="ICS")throw std::runtime_error("Invalid process");
        pairProcess=process=="PP";primaryEnergy=std::stod(argv[3]);const double emin=std::stod(argv[4]),tol=std::stod(argv[5]);
        if(!std::isfinite(primaryEnergy)||!std::isfinite(emin)||!std::isfinite(tol)||primaryEnergy<=emin||emin<=0||tol<=0||tol>.01)throw std::runtime_error("Invalid physics parameters");
        threshold=emin/primaryEnergy;
        std::ifstream edgesIn(argv[6]);std::vector<Real> edges;Real edge;while(edgesIn>>edge)edges.push_back(edge);
        if(!edgesIn.eof()||edges.size()<3||edges.front()!=0||edges.back()!=1||!std::is_sorted(edges.begin(),edges.end())||std::adjacent_find(edges.begin(),edges.end())!=edges.end())throw std::runtime_error("Invalid frozen energy bin edges");
        observables=secondary::regions(pairProcess,edges);
        cosmology.Init(.1,.692,67.8);
        const double t=2.7255*units.K;CompoundBackground bg;bg.AddComponent(new PlankBackground(t,1e-3*t,1e3*t,0.,.1));
        PBackgroundIntegral background=new ContinuousBackgroundIntegral(bg,.01,std::pow(10.,.01),.1,tol);
        ICS ics(background,emin*units.eV);GammaPP pp(background);Particle primary(pairProcess?Photon:Electron,0.);primary.Energy=primaryEnergy*units.eV;
        if(oracleMode) {
            cells=std::stoi(argv[7]);order=std::stoi(argv[8]);destination=argv[9];
            if(cells<8||order<8||order>64)throw std::runtime_error("Invalid oracle quadrature");
            Randomizer rng(1);double s;
            if(!(pairProcess?pp.SampleS(primary,s,rng):ics.SampleS(primary,s,rng)))throw std::runtime_error("Cannot construct reference kernel");
            return 0;
        }
        const std::string mode=argv[7];if(mode!="old"&&mode!="current")throw std::runtime_error("Invalid sampler");oldMode=mode=="old";
        const auto count=std::stoull(argv[8]);const auto seed=std::stoul(argv[9]);if(!count||!seed||seed>0xffffffffUL)throw std::runtime_error("Invalid count/seed");
        const unsigned short endian=1;if(*reinterpret_cast<const unsigned char*>(&endian)!=1||sizeof(double)!=8||!std::numeric_limits<double>::is_iec559)throw std::runtime_error("Raw format requires little-endian IEEE754 binary64");
        Randomizer rng(seed);std::ofstream raw(argv[10],std::ios::binary);raw.write("MCRAYEM1",8);
        std::vector<std::array<unsigned long long,3>> counts(observables.size());
        secondary::Audit aggregate;unsigned long long envelopeEvents=0,narrowingEvents=0,argmaxEvents=0;
        double maxEnergyResidual=0.,maxSupportResidual=0.;
        const Real m2=std::pow(Real(Particle::Mass(Electron)),2);
        for(unsigned long long i=0;i<count;++i) {
            secondary::audit={};double s=std::numeric_limits<double>::quiet_NaN();
            std::vector<Particle> products;
            try {
            if(!(pairProcess?pp.SampleS(primary,s,rng):ics.SampleS(primary,s,rng)))throw std::runtime_error("Secondary sampling failed");
            if(pairProcess)pp.SampleSecondaries(primary,products,s,rng);else ics.SampleSecondaries(primary,products,s,rng);
            if(products.size()!=2)throw std::runtime_error("Wrong secondary multiplicity");
            const double first=products[0].Energy/primary.Energy,second=products[1].Energy/primary.Energy;
            const auto& a=secondary::audit;
            const Real sd=Real(s)/m2;
            binaryDouble(raw,double(sd));binaryDouble(raw,first);binaryDouble(raw,second);
            const unsigned char type=static_cast<unsigned char>(products[0].Type),flags=(a.envelope?1:0)|(a.narrowing?2:0)|(a.argmax?4:0);
            raw.put(type);raw.put(flags);
            if(!std::isfinite(first)||!std::isfinite(second)||first<=0||second<=0||first>=1||second>=1)throw std::runtime_error("Nonphysical secondary energy");
            const Real lower=pairProcess?secondary::pairMinimum(sd):threshold,upper=pairProcess?.5L:1-1/sd;
            const Real violation=std::max({Real(0),lower-first,Real(first)-upper});
            maxSupportResidual=std::max(maxSupportResidual,double(violation));
            maxEnergyResidual=std::max(maxEnergyResidual,std::abs(first+second-1));
            if(violation>32*std::numeric_limits<double>::epsilon()||std::abs(first+second-1)>32*std::numeric_limits<double>::epsilon())throw std::runtime_error("Kinematic or energy-conservation violation");
            const bool firstElectron=products[0].Type==Electron;
            if(pairProcess && !((products[0].Type==Electron&&products[1].Type==Positron)||(products[0].Type==Positron&&products[1].Type==Electron)))throw std::runtime_error("Invalid pair types");
            if(!pairProcess && !(products[0].Type==Photon&&products[1].Type==Electron))throw std::runtime_error("Invalid ICS types");
            for(size_t j=0;j<observables.size();++j)++counts[j][observables[j].multiplicity(first,second,firstElectron)];
            aggregate.proposals+=a.proposals;aggregate.envelope+=a.envelope;aggregate.narrowing+=a.narrowing;aggregate.argmax+=a.argmax;aggregate.largestEnvelopeRatio=std::max(aggregate.largestEnvelopeRatio,a.largestEnvelopeRatio);
            envelopeEvents+=a.envelope>0;narrowingEvents+=a.narrowing>0;argmaxEvents+=a.argmax>0;
            if((i+1)%4096==0)raw.flush();
            if(!raw)throw std::runtime_error("Cannot preserve raw events");
            } catch(...) {
                raw.flush();std::cerr<<std::setprecision(17)<<"Failure at event="<<i<<" seed="<<seed<<" s="<<s<<" products="<<products.size();
                for(const auto& product:products)std::cerr<<" type="<<product.Type<<" E_fraction="<<product.Energy/primary.Energy;
                std::cerr<<'\n';throw;
            }
        }
        raw.close();if(!raw)throw std::runtime_error("Cannot finalize raw events");
        std::ofstream out(argv[11]);out<<std::setprecision(17)<<"# N="<<count<<" seed="<<seed<<" sampler="<<mode<<" process="<<process<<" E_eV="<<primaryEnergy<<" Emin_eV="<<emin<<" tolerance="<<tol<<" proposals="<<aggregate.proposals<<" envelope_exceedances="<<aggregate.envelope<<" envelope_events="<<envelopeEvents<<" max_envelope_ratio="<<aggregate.largestEnvelopeRatio<<" range_narrowings="<<aggregate.narrowing<<" narrowing_events="<<narrowingEvents<<" argmax_returns="<<aggregate.argmax<<" argmax_events="<<argmaxEvents<<" max_energy_residual="<<maxEnergyResidual<<" max_support_residual="<<maxSupportResidual<<"\n# id n0 n1 n2\n";
        for(size_t j=0;j<counts.size();++j)out<<observables[j].id<<'\t'<<counts[j][0]<<'\t'<<counts[j][1]<<'\t'<<counts[j][2]<<'\n';
        if(!out)throw std::runtime_error("Cannot finalize spectrum statistics");
        return 0;
    }catch(Exception* e){std::cerr<<e->Message()<<'\n';delete e;}
     catch(const std::exception& e){std::cerr<<e.what()<<'\n';}
     catch(const char* e){std::cerr<<e<<'\n';}
    return 1;
}
