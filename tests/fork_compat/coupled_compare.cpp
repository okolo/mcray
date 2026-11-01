#include "coupled_api.h"
#include "GridMagneticField.h"
#include "ParticleStack.h"
#include "PropagationEngine.h"
#include "SpatialGeometry.h"
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
using namespace mcray;
using namespace Interactions;
namespace {
class ConstantField : public MagneticField {
public:
    ConstantField(const double* field,double scale): field_(field,field+3),scale_(scale) {}
    void GetValueGauss(const coord_type*,const CosmoTime&,std::vector<double>& b) const override { b=field_; }
    double MinVariabilityScale(const CosmoTime&) const override { return scale_; }
private:
    std::vector<double> field_;
    double scale_;
};
class CountedField : public MagneticField {
public:
    explicit CountedField(MagneticField& field): field_(field) {}
    void GetValueGauss(const coord_type* x,const CosmoTime& t,std::vector<double>& b) const override {
        ++values; field_.GetValueGauss(x,t,b);
    }
    double MinVariabilityScale(const CosmoTime& t) const override { ++scales; return field_.MinVariabilityScale(t); }
    mutable unsigned long values=0,scales=0;
private:
    MagneticField& field_;
};
class Capture : public Utils::TSmartReferencedObj<IOutput> {
public:
    std::vector<Particle> particles;
    void AddParticle(const Particle& p) override { particles.push_back(p); }
    void Flush(bool) override {}
};
struct MoveLogger : ILogger {
    unsigned long calls=0;
    void log(const Particle&) override {
        // Run logs the initial state once, then each complete engine move.
        if(++calls>100001) throw std::runtime_error("Modern coupled move budget exceeded");
    }
};
CoupledOutput Modern(const CoupledInput& input,const char* grid,const double* field,double scale) {
    const ParticleType types[]{Electron,Positron,Proton};
    Particle p(types[input.particle],.01); p.Weight=1; p.Energy=input.energy;
    p.Time.setT(input.observer.initialTime);
    for(int j=0;j<3;++j) { p.X[j]=input.observer.position[j]; p.Pdir[j]=input.observer.direction[j]; }
    const auto& o=input.observer;
    auto observer=std::make_shared<SphereObserver>(SpatialVector{o.center[0],o.center[1],o.center[2]},o.radius);
    Result result; result.SetObserver(observer);
    CosmoTime end; end.setT(o.initialTime+o.duration); result.SetEndTime(end);
    Utils::SmartPtr<Capture> capture=new Capture; result.AddOutput(capture);
    Utils::SmartPtr<MagneticField> owned=input.fieldKind==2
        ? static_cast<MagneticField*>(new SavedMF(grid,SavedMF::FILLING_ZERO))
        : static_cast<MagneticField*>(new ConstantField(field,scale));
    Utils::SmartPtr<CountedField> counted=new CountedField(*owned);
    Utils::SmartPtr<Deflection3D> deflector=new Deflection3D(counted,input.accuracy);
    MoveLogger logger;
    ParticleStack stack; stack.AddPrimary(p);
    PropagationEngine engine(stack,result,17);
    engine.SetLogger(&logger);
    if(input.fieldKind) engine.AddInteraction(PDeflectionInteraction(deflector.operator->()));
    CoupledOutput out{};
    out.transport.beta=p.beta(); out.transport.mass=p.Mass(); out.transport.charge=p.ElectricCharge();
    out.transport.globalScale=owned->MinVariabilityScale(p.Time);
    out.transport.initialScale=out.transport.globalScale;
    out.transport.initialRate=input.fieldKind?deflector->Rate(p):0;
    engine.Run(0); result.Flush();
    if(engine.FailureCount() || capture->particles.size()!=1)
        throw std::runtime_error("Modern coupled trajectory failed or returned incomplete output");
    const auto& final=capture->particles.front();
    if(final.Weight!=1 || final.Energy!=p.Energy || final.Type!=p.Type)
        throw std::runtime_error("Modern metadata changed without interactions");
    out.hit=observer->Reached(final); out.pastEnd=final.Time.t()>end.t();
    out.endError=final.Time.t()-end.t(); out.moves=logger.calls-1;
    out.transport.elapsed=final.Time.t()-o.initialTime; out.transport.energy=final.Energy;
    out.transport.fieldCalls=counted->values; out.transport.scaleCalls=counted->scales;
    for(int j=0;j<3;++j) { out.transport.position[j]=final.X[j]; out.transport.direction[j]=final.Pdir[j]; }
    return out;
}
long double Difference(const long double* a,const long double* b) {
    return std::hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
}
void Print(const CoupledOutput& out) {
    const auto& p=out.transport;
    std::cout << "{\"hit\":" << out.hit << ",\"past_end\":" << out.pastEnd
              << ",\"end_error_Mpc\":" << out.endError/units.Mpc << ",\"moves\":" << out.moves
              << ",\"flight_Mpc\":" << p.elapsed/units.Mpc << ",\"beta\":" << p.beta
              << ",\"mass\":" << p.mass << ",\"charge\":" << p.charge << ",\"energy\":" << p.energy
              << ",\"initial_rate\":" << p.initialRate << ",\"scale_Mpc\":" << p.globalScale/units.Mpc
              << ",\"field_calls\":" << p.fieldCalls << ",\"scale_calls\":" << p.scaleCalls
              << ",\"zero_scale_calls\":" << p.zeroScaleCalls << ",\"position_Mpc\":[";
    for(int j=0;j<3;++j) std::cout << (j?",":"") << p.position[j]/units.Mpc;
    std::cout << "],\"direction\":[";
    for(int j=0;j<3;++j) std::cout << (j?",":"") << p.direction[j];
    std::cout << "],\"direction_norm_error\":" << std::abs(std::hypot(p.direction[0],p.direction[1],p.direction[2])-1) << '}';
}
void Run(const char* cases,const char* grid) {
    cosmology.Init(10.01,.7,70);
    std::ifstream input(cases); if(!input) throw std::runtime_error("Missing coupled input");
    std::string profile,ray; int fieldKind,particle,accuracy;
    double center[3]; long double direction[3];
    while(input >> profile >> ray >> fieldKind >> particle >> accuracy
                >> center[0] >> center[1] >> center[2] >> direction[0] >> direction[1] >> direction[2]) {
        CoupledInput c{}; c.particle=particle; c.accuracy=accuracy; c.fieldKind=fieldKind; c.energy=1e18*units.eV;
        c.observer.mode=1; c.observer.radius=.01*units.Mpc;
        c.observer.initialTime=cosmology.z2t0(.01); c.observer.duration=.2*units.Mpc;
        for(int j=0;j<3;++j) { c.observer.center[j]=center[j]*units.Mpc; c.observer.direction[j]=direction[j]; }
        const double field[]{2e-10,3e-10,5e-10};
        const auto old=HistoricalCoupled(c,grid,field,.1*units.Mpc);
        const auto modern=Modern(c,grid,field,.1*units.Mpc);
        std::cout << "{\"profile\":\"" << profile << "\",\"ray\":\"" << ray << "\",\"field\":" << fieldKind
                  << ",\"particle\":" << particle << ",\"accuracy\":" << accuracy << ",\"fork\":";
        Print(old); std::cout << ",\"modern\":"; Print(modern);
        std::cout << ",\"same_hit\":" << (old.hit==modern.hit)
                  << ",\"difference_position_Mpc\":" << Difference(old.transport.position,modern.transport.position)/units.Mpc
                  << ",\"difference_direction\":" << Difference(old.transport.direction,modern.transport.direction)
                  << ",\"difference_time_Mpc\":" << std::abs(old.transport.elapsed-modern.transport.elapsed)/units.Mpc << "}\n";
    }
    if(!input.eof()) throw std::runtime_error("Malformed coupled input");
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=3) throw std::runtime_error("Usage: coupled-probe CASES GRID");
        std::cout << std::setprecision(21) << std::boolalpha;
        Run(argv[1],argv[2]); return 0;
    } catch(Utils::Exception* error) {
        std::unique_ptr<Utils::Exception> owner(error); std::cerr << error->Message() << '\n';
    } catch(const std::exception& error) { std::cerr << error.what() << '\n'; }
    return 1;
}
