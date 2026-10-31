#include "observer_api.h"
#include "SpatialGeometry.h"
#include "ParticleStack.h"
#include "PropagationEngine.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
using namespace mcray;
namespace {
class Capture : public Utils::TSmartReferencedObj<IOutput> {
public:
    std::vector<Particle> particles;
    void AddParticle(const Particle& p) override { particles.push_back(p); }
    void Flush(bool) override {}
};
std::shared_ptr<SphereObserver> Observer(const ObserverInput& input) {
    if(input.mode==-1) return {};
    return std::make_shared<SphereObserver>(
        SpatialVector{input.center[0],input.center[1],input.center[2]},input.radius,
        input.mode==1?SphereObserver::StopWhen::Inside:SphereObserver::StopWhen::Outside);
}
Particle Initial(const ObserverInput& input) {
    Particle p(Photon,.01); p.Energy=1e12*units.eV; p.Weight=2.5;
    p.Time.setT(input.initialTime);
    for(int i=0;i<3;++i) { p.X[i]=input.position[i]; p.Pdir[i]=input.direction[i]; }
    return p;
}
ObserverFlight ModernFlight(const ObserverInput& input,bool variable) {
    auto observer=Observer(input);
    Result result;
    result.SetObserver(observer);
    CosmoTime end; end.setT(input.initialTime+input.duration);
    result.SetEndTime(end);
    Utils::SmartPtr<Capture> capture=new Capture;
    result.AddOutput(capture);
    auto p=Initial(input);
    ParticleStack stack; stack.AddPrimary(p);
    PropagationEngine engine(stack,result,17);
    engine.Run(variable?.1:0);
    if(engine.FailureCount()) throw std::runtime_error("Ballistic engine reported failed histories");
    result.Flush();
    if(capture->particles.size()!=1) throw std::runtime_error("Expected one final particle");
    const auto& final=capture->particles.front();
    if(final.Weight!=p.Weight || final.Energy!=p.Energy || final.Type!=Photon)
        throw std::runtime_error("Ballistic metadata changed");
    ObserverFlight out{};
    out.reached=observer && observer->Reached(final);
    out.elapsed=final.Time.t()-input.initialTime;
    out.pastEnd=final.Time.t()>input.initialTime+input.duration;
    const auto target=input.initialTime+input.duration;
    out.endError=final.Time.t()-target;
    out.endErrorUlps=out.endError/(std::nextafter(target,std::numeric_limits<cosmo_time>::infinity())-target);
    for(int i=0;i<3;++i) out.position[i]=final.X[i];
    return out;
}
void PrintFlight(const ObserverFlight& out) {
    std::cout << "{\"hit\":" << out.reached << ",\"past_end\":" << out.pastEnd
              << ",\"flight_Mpc\":" << out.elapsed/units.Mpc
              << ",\"end_error_Mpc\":" << out.endError/units.Mpc
              << ",\"end_error_ulps\":" << out.endErrorUlps
              << ",\"position_Mpc\":[";
    for(int i=0;i<3;++i) {
        if(!std::isfinite(out.position[i])) throw std::runtime_error("Nonfinite endpoint");
        std::cout << (i?",":"") << out.position[i]/units.Mpc;
    }
    std::cout << "]}";
}
void Run(const char* file) {
    cosmology.Init(10.01,.7,70);
    std::ifstream stream(file);
    if(!stream) throw std::runtime_error("Missing observer cases");
    std::string name;
    int kind,mode,point;
    double radiusMpc,centerMpc[3],durationR;
    long double start[3],direction[3];
    while(stream >> name >> kind >> mode >> point >> radiusMpc
                 >> centerMpc[0] >> centerMpc[1] >> centerMpc[2]
                 >> start[0] >> start[1] >> start[2]
                 >> direction[0] >> direction[1] >> direction[2] >> durationR) {
        ObserverInput input{}; input.mode=mode;
        input.radius=radiusMpc*units.Mpc;
        input.duration=durationR*input.radius;
        input.initialTime=kind==0?-32*static_cast<long double>(input.radius):cosmology.z2t0(.01);
        for(int i=0;i<3;++i) {
            input.center[i]=centerMpc[i]*units.Mpc;
            input.position[i]=input.center[i]+start[i]*input.radius;
            input.direction[i]=direction[i];
        }
        if(point==1 || point==2) {
            const auto toward=point==1?-std::numeric_limits<long double>::infinity():
                                      std::numeric_limits<long double>::infinity();
            input.position[0]=std::nextafter(input.position[0],toward);
        }
        std::cout << "{\"case\":\"" << name << "\",\"kind\":" << kind
                  << ",\"mode\":" << mode << ",\"radius_Mpc\":" << radiusMpc
                  << ",\"center_Mpc\":[" << centerMpc[0] << ',' << centerMpc[1] << ',' << centerMpc[2]
                  << "],\"requested_Mpc\":" << input.duration/units.Mpc;
        std::cout << ",\"initial_position_Mpc\":[";
        for(int i=0;i<3;++i) std::cout << (i?",":"") << input.position[i]/units.Mpc;
        std::cout << ']';
        if(kind==0) {
            const auto old=HistoricalObserverPoint(input);
            auto observer=Observer(input); auto p=Initial(input);
            const bool hit=observer && observer->Reached(p);
            const auto step=observer?observer->LimitStep(p,input.duration):input.duration;
            if(!std::isfinite(old.step) || !std::isfinite(step)) throw std::runtime_error("Nonfinite observer step");
            std::cout << ",\"fork_hit\":" << old.reached << ",\"modern_hit\":" << hit
                      << ",\"fork_step_Mpc\":" << old.step/units.Mpc
                      << ",\"modern_step_Mpc\":" << step/units.Mpc
                      << ",\"fork_exceeds_remaining\":" << (old.step>input.duration)
                      << ",\"modern_exceeds_remaining\":" << (step>input.duration);
        } else {
            const auto old=HistoricalObserverFlight(input);
            const auto modern=ModernFlight(input,false),variable=ModernFlight(input,true);
            bool same=modern.reached==variable.reached && modern.elapsed==variable.elapsed;
            for(int i=0;i<3;++i) same=same && modern.position[i]==variable.position[i];
            std::cout << ",\"fork\":"; PrintFlight(old);
            std::cout << ",\"modern\":"; PrintFlight(modern);
            std::cout << ",\"variable_rate\":"; PrintFlight(variable);
            std::cout << ",\"fork_steps\":" << old.steps;
            std::cout << ",\"modern_paths_equal\":" << same;
            long double squared=0;
            for(int i=0;i<3;++i) squared+=(old.position[i]-modern.position[i])*(old.position[i]-modern.position[i]);
            std::cout << ",\"difference_position_Mpc\":" << std::sqrt(squared)/units.Mpc
                      << ",\"difference_flight_Mpc\":" << std::abs(old.elapsed-modern.elapsed)/units.Mpc
                      << ",\"difference_rate_path_time_Mpc\":" << std::abs(modern.elapsed-variable.elapsed)/units.Mpc;
        }
        std::cout << "}\n";
    }
    if(!stream.eof()) throw std::runtime_error("Malformed observer input");
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=2) throw std::runtime_error("Usage: fork-observer-probe CASES");
        std::cout << std::setprecision(21) << std::boolalpha;
        Run(argv[1]); return 0;
    } catch(Utils::Exception* e) {
        std::unique_ptr<Utils::Exception> owner(e); std::cerr << e->Message() << '\n';
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; }
    return 1;
}
