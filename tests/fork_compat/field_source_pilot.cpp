// Source cost/failure diagnostic only. No intervals or equivalence decision.
#include "coupled_modern.h"
#include "fork_api.h"
#include "SpatialGeometry.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>

using namespace mcray;
namespace {
void Run(const char* cases,const char* grid,bool eev) {
    cosmology.Init(10.01,.7,70);
    ForkUnits(units.Mpc);
    const double unusedField[]{0,0,0};
    // Exclude grid parsing from sampling/propagation timers. Both fields are
    // immutable and shared sequentially; per-history counters remain separate.
    auto field=ModernCoupledField(2,grid,unusedField,.1*units.Mpc);
    std::unique_ptr<void,void(*)(void*)> referenceField(
        HistoricalCoupledField(2,grid,unusedField,.1*units.Mpc),DeleteHistoricalCoupledField);
    std::ifstream input(cases);
    if(!input) throw std::runtime_error("Missing field/source pilot design");
    std::string name;
    int species,accuracy;
    double degrees,theta,phi;
    unsigned long oldSeed,newSeed,count;
    while(input >> name >> species >> degrees >> theta >> phi >> oldSeed >> newSeed >> count >> accuracy) {
        if((species!=0 && species!=3 && (!eev || species!=2)) || (count!=2 && count!=1000) || accuracy!=20)
            throw std::runtime_error("Unsupported field/source pilot settings");
        const auto axis=SourceGeometry::NormalizeAxis(
            {std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)});
        SourceGeometry source({},axis,degrees==0?SourceGeometry::AngularModel::Pencil:
                                                 SourceGeometry::AngularModel::Cone,degrees);
        CoupledInput state{};
        state.particle=species; state.accuracy=accuracy; state.fieldKind=2;
        state.energy=(eev?1e18:1e15)*units.eV;
        auto& observer=state.observer;
        observer.mode=1; observer.radius=.01*units.Mpc;
        observer.initialTime=cosmology.z2t0(.01); observer.duration=.2*units.Mpc;
        for(int j=0;j<3;++j) observer.center[j]=static_cast<double>(.1*units.Mpc*axis[j]);
        for(bool modern:{false,true}) {
            Randomizer random(newSeed);
            std::unique_ptr<void,void(*)(void*)> cone(
                ForkCone(degrees*std::acos(-1.)/180,theta,phi,oldSeed),ForkDeleteCone);
            double seconds=0;
            for(unsigned long i=0;i<count;++i) {
                const auto before=std::chrono::steady_clock::now();
                double weight;
                if(modern) {
                    const ParticleType types[]{Electron,Positron,Proton,Photon};
                    Particle particle(types[species],.01);
                    particle.Energy=state.energy; particle.Weight=1;
                    source.Initialize(particle,random);
                    weight=particle.Weight;
                    for(int j=0;j<3;++j) {
                        observer.position[j]=particle.X[j]; observer.direction[j]=particle.Pdir[j];
                    }
                } else {
                    const auto sample=ForkNext(cone.get()); weight=sample.weight;
                    for(int j=0;j<3;++j) { observer.position[j]=0; observer.direction[j]=sample.direction[j]; }
                }
                if(weight!=1) throw std::runtime_error("Unexpected ordinary-cone weight");
                const auto flight=modern?ModernCoupled(state,*field,weight):
                                         HistoricalCoupledWithField(state,referenceField.get());
                seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
                std::cout << "{\"record\":\"history\",\"case\":\"" << name
                          << "\",\"implementation\":\"" << (modern?"modern":"fork")
                          << "\",\"primary\":" << i << ",\"weight\":" << weight << ",\"direction\":[";
                for(int j=0;j<3;++j) std::cout << (j?",":"") << observer.direction[j];
                std::cout << "],\"flight\":"; PrintCoupled(flight);
                std::cout << ",\"historical_steps\":" << (modern?0:flight.moves) << "}\n";
            }
            std::cout << "{\"record\":\"timing\",\"case\":\"" << name
                      << "\",\"implementation\":\"" << (modern?"modern":"fork")
                      << "\",\"compute_seconds\":" << seconds << "}\n";
            std::cerr << "completed " << name << ' ' << (modern?"modern":"fork") << ' ' << count << '\n';
        }
    }
    if(!input.eof()) throw std::runtime_error("Malformed field/source design");
    std::cout.flush();
    if(!std::cout) throw std::runtime_error("Field/source raw output failed");
}
}
int main(int argc,char** argv) {
    try {
        if((argc!=3 && argc!=4) || (argc==4 && std::string(argv[3])!="pev" && std::string(argv[3])!="eev"))
            throw std::runtime_error("Usage: field-source-pilot CASES GRID [pev|eev]");
        std::cout << std::setprecision(21) << std::boolalpha;
        Run(argv[1],argv[2],argc==4 && std::string(argv[3])=="eev"); return 0;
    } catch(Utils::Exception* error) {
        std::unique_ptr<Utils::Exception> owner(error); std::cerr << error->Message() << '\n';
    } catch(const std::exception& error) { std::cerr << error.what() << '\n'; }
    return 1;
}
