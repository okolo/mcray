// One uninterrupted source stream; Python computes the frozen confidence family.
#include "coupled_modern.h"
#include "fork_api.h"
#include "SpatialGeometry.h"
#include <chrono>
#include <array>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>

using namespace mcray;
namespace {
struct Sum {
    long double value=0,correction=0;
    void Add(long double x) {
        const auto adjusted=x-correction, next=value+adjusted;
        correction=(next-value)-adjusted; value=next;
    }
};
struct Statistics {
    unsigned long count=0,hits=0,pastEnd=0,maxMoves=0;
    long double positiveEndError=0;
    std::array<Sum,14> sums;
    void Add(unsigned long id,const ObserverInput& input,const CoupledOutput& out,double weight) {
        if(id!=count || weight!=1 || !std::isfinite(out.transport.elapsed) ||
           out.transport.elapsed<0 || !std::isfinite(out.endError))
            throw std::runtime_error("Invalid primary sequence, weight or time");
        const auto elapsed=out.transport.elapsed/units.Mpc;
        sums[0].Add(out.hit);
        const std::array<long double,4> thresholds{.092L,.094L,.096L,.098L};
        for(size_t k=0;k<4;++k) sums[1+k].Add(out.hit && elapsed<=thresholds[k]);
        for(size_t j=0;j<3;++j) {
            const std::array<long double,3> values{out.transport.position[j]/(.25L*units.Mpc),
                out.transport.direction[j],input.direction[j]};
            for(size_t k=0;k<3;++k) {
                if(!std::isfinite(values[k]) || std::abs(values[k])>1)
                    throw std::runtime_error("Observable outside frozen bounds");
                sums[5+3*k+j].Add(values[k]);
            }
        }
        ++count; hits+=out.hit; pastEnd+=out.pastEnd;
        maxMoves=std::max(maxMoves,out.moves);
        positiveEndError=std::max(positiveEndError,out.endError/units.Mpc);
    }
};
void Run(const char* cases,const std::string& selected,bool modern,const char* summary) {
    const auto gridPath=std::string(cases).substr(0,std::string(cases).find_last_of('/'))+"/transport-fixture/transport.grid";
    const char* grid=gridPath.c_str();
    bool found=false;
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
        if(name!=selected) continue;
        if(found) throw std::runtime_error("Duplicate ensemble case");
        found=true;
        if((species!=0 && species!=3 && species!=2) || (count!=1000 && count!=9000000) || accuracy!=40)
            throw std::runtime_error("Unsupported field/source ensemble settings");
        const auto axis=SourceGeometry::NormalizeAxis(
            {std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)});
        SourceGeometry source({},axis,degrees==0?SourceGeometry::AngularModel::Pencil:
                                                 SourceGeometry::AngularModel::Cone,degrees);
        CoupledInput state{};
        state.particle=species; state.accuracy=accuracy; state.fieldKind=2;
        state.energy=1e18*units.eV;
        auto& observer=state.observer;
        observer.mode=1; observer.radius=.01*units.Mpc;
        observer.initialTime=cosmology.z2t0(.01); observer.duration=.2*units.Mpc;
        for(int j=0;j<3;++j) observer.center[j]=static_cast<double>(.1*units.Mpc*axis[j]);
        {
            Randomizer random(newSeed);
            std::unique_ptr<void,void(*)(void*)> cone(
                ForkCone(degrees*std::acos(-1.)/180,theta,phi,oldSeed),ForkDeleteCone);
            const auto started=std::chrono::steady_clock::now();
            Statistics statistics;
            for(unsigned long i=0;i<count;++i) {
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
                statistics.Add(i,observer,flight,weight);
                std::cout << "{\"record\":\"history\",\"case\":\"" << name
                          << "\",\"implementation\":\"" << (modern?"modern":"fork")
                          << "\",\"primary\":" << i << ",\"weight\":" << weight << ",\"direction\":[";
                for(int j=0;j<3;++j) std::cout << (j?",":"") << observer.direction[j];
                std::cout << "],\"flight\":"; PrintCoupled(flight);
                std::cout << ",\"historical_steps\":" << (modern?0:flight.moves) << "}\n";
                if((i+1)%100000==0) std::cerr << "processed " << (i+1) << "/" << count << '\n';
            }
            std::cout.flush();
            if(!std::cout) throw std::runtime_error("Raw output failed");
            std::ofstream out(summary);
            out << std::setprecision(21) << "{\"case\":\"" << name << "\",\"implementation\":\""
                << (modern?"modern":"fork") << "\",\"histories\":" << statistics.count
                << ",\"hits\":" << statistics.hits << ",\"failures\":0,\"emitted_weight\":" << statistics.count
                << ",\"past_end_histories\":" << statistics.pastEnd
                << ",\"max_positive_end_error_Mpc\":" << statistics.positiveEndError
                << ",\"max_engine_moves\":" << statistics.maxMoves << ",\"sums\":[";
            for(size_t k=0;k<14;++k) out << (k?",":"") << statistics.sums[k].value;
            out << "],\"wall_seconds\":" << std::chrono::duration<double>(
                std::chrono::steady_clock::now()-started).count() << "}\n";
            out.flush();
            if(!out) throw std::runtime_error("Summary output failed");
        }
    }
    if(!input.eof() || !found) throw std::runtime_error("Malformed field/source design");
    std::cout.flush();
    if(!std::cout) throw std::runtime_error("Field/source raw output failed");
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=5 || (std::string(argv[3])!="fork" && std::string(argv[3])!="modern"))
            throw std::runtime_error("Usage: eev-field-ensemble CASES CASE fork|modern SUMMARY");
        std::ios::sync_with_stdio(false);
        std::cout << std::setprecision(21) << std::boolalpha;
        Run(argv[1],argv[2],std::string(argv[3])=="modern",argv[4]); return 0;
    } catch(Utils::Exception* error) {
        std::unique_ptr<Utils::Exception> owner(error); std::cerr << error->Message() << '\n';
    } catch(const std::exception& error) { std::cerr << error.what() << '\n'; }
    return 1;
}
