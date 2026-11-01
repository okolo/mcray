// Cost/failure pilot only: independent source streams, no acceptance verdict.
#include "fork_api.h"
#include "observer_modern.h"
#include "SpatialGeometry.h"
#include "Utils.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

using namespace mcray;
using namespace observer_validation;
namespace {
void Run(const char* file) {
    cosmology.Init(10.01,.7,70);
    ForkUnits(units.Mpc);
    std::ifstream input(file);
    if(!input) throw std::runtime_error("Missing pilot design");
    std::string name;
    double degrees,theta,phi;
    unsigned long forkSeed,modernSeed;
    size_t count;
    while(input >> name >> degrees >> theta >> phi >> forkSeed >> modernSeed >> count) {
        if(count!=1000) throw std::runtime_error("Cost pilot requires exactly 1000 histories per case");
        const auto axis=SourceGeometry::NormalizeAxis(
            {std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)});
        SourceGeometry source({},axis,degrees==0?SourceGeometry::AngularModel::Pencil:
                                                 SourceGeometry::AngularModel::Cone,degrees);
        Randomizer random(modernSeed);
        // One mutable historical cone/RNG per case, as in the original source loop.
        std::unique_ptr<void,void(*)(void*)> reference(
            ForkCone(degrees*std::acos(-1.)/180,theta,phi,forkSeed),ForkDeleteCone);
        ObserverInput state{};
        state.mode=1; state.radius=.01*units.Mpc;
        state.initialTime=cosmology.z2t0(.01); state.duration=.2*units.Mpc;
        // Share the represented center, rather than independently rotating it.
        for(int j=0;j<3;++j) state.center[j]=static_cast<double>(.1*units.Mpc*axis[j]);
        for(bool modern:{false,true}) {
            double seconds=0;
            for(size_t i=0;i<count;++i) {
                const auto before=std::chrono::steady_clock::now();
                double weight;
                if(modern) {
                    auto particle=Initial(state,1);
                    source.Initialize(particle,random);
                    weight=particle.Weight;
                    for(int j=0;j<3;++j) {
                        state.position[j]=particle.X[j]; state.direction[j]=particle.Pdir[j];
                    }
                } else {
                    const auto sample=ForkNext(reference.get());
                    weight=sample.weight;
                    for(int j=0;j<3;++j) state.direction[j]=sample.direction[j];
                }
                const auto flight=modern?ModernFlight(state,false,weight):HistoricalObserverFlight(state);
                seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
                if(!std::isfinite(weight)) throw std::runtime_error("Nonfinite emitted weight");
                std::cout << "{\"record\":\"history\",\"case\":\"" << name
                          << "\",\"implementation\":\"" << (modern?"modern":"fork")
                          << "\",\"primary\":" << i << ",\"weight\":" << weight << ",\"direction\":[";
                for(int j=0;j<3;++j) std::cout << (j?",":"") << state.direction[j];
                std::cout << "],\"flight\":"; PrintFlight(flight);
                std::cout << ",\"historical_steps\":" << (modern?0:flight.steps) << "}\n";
            }
            // Exclude JSON output/compilation/startup from this compute timer.
            // Retain an independent whole-process wall clock in the Python driver.
            std::cout << "{\"record\":\"timing\",\"case\":\"" << name
                      << "\",\"implementation\":\"" << (modern?"modern":"fork")
                      << "\",\"compute_seconds\":" << seconds << "}\n";
        }
    }
    if(!input.eof()) throw std::runtime_error("Malformed pilot design");
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::runtime_error("Usage: source-observer-pilot CASES");
        std::cout << std::setprecision(21) << std::boolalpha;
        Run(argv[1]); return 0;
    } catch(Utils::Exception* e) {
        std::unique_ptr<Utils::Exception> owner(e); std::cerr << e->Message() << '\n';
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; }
    return 1;
}
