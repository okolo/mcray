#include "observer_modern.h"
#include "SpatialGeometry.h"
#include "ParticleStack.h"
#include "PropagationEngine.h"
#include <limits>
#include <memory>
using namespace mcray;
namespace observer_validation {
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
Particle Initial(const ObserverInput& input,double weight) {
    Particle p(Photon,.01); p.Energy=1e12*units.eV; p.Weight=weight;
    p.Time.setT(input.initialTime);
    for(int i=0;i<3;++i) { p.X[i]=input.position[i]; p.Pdir[i]=input.direction[i]; }
    return p;
}
ObserverFlight ModernFlight(const ObserverInput& input,bool variable,double weight) {
    auto observer=Observer(input);
    Result result;
    result.SetObserver(observer);
    CosmoTime end; end.setT(input.initialTime+input.duration);
    result.SetEndTime(end);
    Utils::SmartPtr<Capture> capture=new Capture;
    result.AddOutput(capture);
    auto p=Initial(input,weight);
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
}
