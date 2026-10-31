// Independent processes retain one uninterrupted source stream per case/path.
#include "fork_api.h"
#include "observer_modern.h"
#include "SpatialGeometry.h"
#include "Utils.h"
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace mcray;
using namespace observer_validation;
namespace {
constexpr std::array<long double,4> thresholds{.092L,.094L,.096L,.098L};
struct Sum {
    long double value=0, correction=0;
    void Add(long double x) {
        const auto adjusted=x-correction;
        const auto next=value+adjusted;
        correction=(next-value)-adjusted;
        value=next;
    }
};
struct Statistics {
    unsigned long count=0,hits=0,pastEnd=0,maxSteps=0;
    long double positiveEndError=0;
    std::array<Sum,11> sums;
    void Add(unsigned long index,const ObserverInput& state,const ObserverFlight& flight,double weight) {
        if(index!=count || weight!=1) throw std::runtime_error("Invalid primary sequence or weight");
        if(!std::isfinite(flight.elapsed) || flight.elapsed<0 ||
           !std::isfinite(flight.endError) || !std::isfinite(flight.endErrorUlps))
            throw std::runtime_error("Invalid flight time");
        sums[0].Add(flight.reached);
        const auto elapsed=flight.elapsed/units.Mpc;
        for(size_t k=0;k<thresholds.size();++k)
            sums[1+k].Add(flight.reached && elapsed<=thresholds[k]);
        for(size_t j=0;j<3;++j) {
            const auto position=flight.position[j]/(.25L*units.Mpc);
            if(!std::isfinite(position) || std::abs(position)>1 ||
               !std::isfinite(state.direction[j]) || std::abs(state.direction[j])>1)
                throw std::runtime_error("Observable outside frozen bounds");
            sums[5+j].Add(position);
            sums[8+j].Add(state.direction[j]);
        }
        ++count; hits+=flight.reached; pastEnd+=flight.pastEnd;
        maxSteps=std::max(maxSteps,flight.steps);
        positiveEndError=std::max(positiveEndError,flight.endError/units.Mpc);
    }
};
void Run(const char* file,const std::string& selected,bool modern,const char* summary) {
    cosmology.Init(10.01,.7,70);
    ForkUnits(units.Mpc);
    std::ifstream input(file);
    if(!input) throw std::runtime_error("Missing ensemble design");
    std::string name;
    double degrees,theta,phi;
    unsigned long forkSeed,modernSeed,count;
    bool found=false;
    while(input >> name >> degrees >> theta >> phi >> forkSeed >> modernSeed >> count) {
        if(name!=selected) continue;
        if(found || count==0) throw std::runtime_error("Invalid case/count");
        found=true;
        const auto axis=SourceGeometry::NormalizeAxis(
            {std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)});
        SourceGeometry source({},axis,degrees==0?SourceGeometry::AngularModel::Pencil:
                                                 SourceGeometry::AngularModel::Cone,degrees);
        Randomizer random(modernSeed);
        std::unique_ptr<void,void(*)(void*)> reference(
            ForkCone(degrees*std::acos(-1.)/180,theta,phi,forkSeed),ForkDeleteCone);
        ObserverInput state{};
        state.mode=1; state.radius=.01*units.Mpc;
        state.initialTime=cosmology.z2t0(.01); state.duration=.2*units.Mpc;
        for(int j=0;j<3;++j) state.center[j]=static_cast<double>(.1*units.Mpc*axis[j]);
        Statistics stats;
        const auto before=std::chrono::steady_clock::now();
        for(unsigned long i=0;i<count;++i) {
            double weight;
            if(modern) {
                auto particle=Initial(state,1);
                source.Initialize(particle,random);
                weight=particle.Weight;
                for(int j=0;j<3;++j) {
                    state.position[j]=particle.X[j]; state.direction[j]=particle.Pdir[j];
                }
            } else {
                const auto sample=ForkNext(reference.get()); weight=sample.weight;
                for(int j=0;j<3;++j) state.direction[j]=sample.direction[j];
            }
            const auto flight=modern?ModernFlight(state,false,weight):HistoricalObserverFlight(state);
            stats.Add(i,state,flight,weight);
            std::cout << "{\"record\":\"history\",\"case\":\"" << name
                      << "\",\"implementation\":\"" << (modern?"modern":"fork")
                      << "\",\"primary\":" << i << ",\"weight\":" << weight << ",\"direction\":[";
            for(int j=0;j<3;++j) std::cout << (j?",":"") << state.direction[j];
            std::cout << "],\"flight\":"; PrintFlight(flight);
            std::cout << ",\"historical_steps\":" << (modern?0:flight.steps) << "}\n";
            if((i+1)%100000==0) std::cerr << "processed " << (i+1) << '/' << count << '\n';
        }
        std::cout.flush();
        if(!std::cout) throw std::runtime_error("Raw history output failed");
        std::ofstream out(summary);
        if(!out) throw std::runtime_error("Cannot create summary");
        out << std::setprecision(21) << "{\"case\":\"" << name << "\",\"implementation\":\""
            << (modern?"modern":"fork") << "\",\"histories\":" << stats.count
            << ",\"hits\":" << stats.hits << ",\"failures\":0,\"emitted_weight\":" << stats.count
            << ",\"past_end_histories\":" << stats.pastEnd
            << ",\"max_positive_end_error_Mpc\":" << stats.positiveEndError
            << ",\"max_historical_steps\":" << stats.maxSteps << ",\"sums\":[";
        for(size_t k=0;k<stats.sums.size();++k) out << (k?",":"") << stats.sums[k].value;
        out << "],\"wall_seconds\":"
            << std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count() << "}\n";
        out.flush();
        if(!out) throw std::runtime_error("Summary output failed");
    }
    if(!input.eof() || !found) throw std::runtime_error("Malformed or missing ensemble case");
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=5 || (std::string(argv[3])!="fork" && std::string(argv[3])!="modern"))
            throw std::runtime_error("Usage: source-observer-ensemble CASES CASE fork|modern SUMMARY");
        std::ios::sync_with_stdio(false);
        std::cout << std::setprecision(21) << std::boolalpha;
        Run(argv[1],argv[2],std::string(argv[3])=="modern",argv[4]); return 0;
    } catch(Utils::Exception* e) {
        std::unique_ptr<Utils::Exception> owner(e); std::cerr << e->Message() << '\n';
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; }
    return 1;
}
