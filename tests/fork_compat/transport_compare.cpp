// Deterministic driver: no interactions, observers or sampled source histories.
#include "GridMagneticField.h"
#include "transport_api.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
using namespace mcray;
using namespace Interactions;
using Vec3 = std::array<long double,3>;
namespace {
class CountedField : public MagneticField {
public:
    explicit CountedField(MagneticField& field): field_(field) {}
    void GetValueGauss(const coord_type* x, const CosmoTime& t, std::vector<double>& out) const override {
        ++values; field_.GetValueGauss(x,t,out);
    }
    double MinVariabilityScale(const CosmoTime& t) const override {
        ++scales; return field_.MinVariabilityScale(t);
    }
    mutable unsigned long values=0,scales=0;
private:
    MagneticField& field_;
};
class ConstantField : public MagneticField {
public:
    ConstantField(const double* b,double scale): value_(b,b+3),scale_(scale) {}
    void GetValueGauss(const coord_type*,const CosmoTime&,std::vector<double>& out) const override { out=value_; }
    double MinVariabilityScale(const CosmoTime&) const override { return scale_; }
private:
    std::vector<double> value_;
    double scale_;
};
long double Norm(const long double* x) { return std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); }
long double Difference(const long double* x,const long double* y) {
    Vec3 d{}; for(int i=0;i<3;++i) d[i]=x[i]-y[i]; return Norm(d.data());
}
TransportOutput Modern(const TransportInput& input,const char* grid,int mode,const double* b,double scale) {
    const ParticleType types[]={Electron,Positron,Proton};
    Particle particle(types[input.particle],0.01);
    particle.Energy=input.energy; particle.Time.setT(input.initialTime);
    for(int i=0;i<3;++i) { particle.X[i]=input.position[i]; particle.Pdir[i]=input.direction[i]; }
    // SmartPtr owns both objects; the counter borrows the underlying field only
    // until the deflector has finished. No stack object enters reference ownership.
    Utils::SmartPtr<MagneticField> owned=grid
        ? static_cast<MagneticField*>(new SavedMF(grid,static_cast<SavedMF::FillingType>(mode)))
        : static_cast<MagneticField*>(new ConstantField(b,scale));
    auto* counted=new CountedField(*owned);
    Deflection3D deflector(counted,input.accuracy);
    Randomizer random(91827);
    TransportOutput result{};
    result.beta=particle.beta(); result.mass=particle.Mass(); result.charge=particle.ElectricCharge();
    result.globalScale=owned->MinVariabilityScale(particle.Time);
    result.initialScale=result.globalScale;
    result.initialRate=deflector.Rate(particle);
    for(int step=0;step<input.outerSteps;++step)
        if(!deflector.Propagate(input.duration/input.outerSteps,particle,random))
            throw std::runtime_error("Charged modern propagation refused");
    for(int i=0;i<3;++i) { result.position[i]=particle.X[i]; result.direction[i]=particle.Pdir[i]; }
    result.elapsed=particle.Time.t()-input.initialTime;
    result.energy=particle.Energy;
    result.fieldCalls=counted->values; result.scaleCalls=counted->scales;
    return result;
}
// Exact constant-field helix for dn/dt=(q e/E) n cross B. Each implementation's
// own beta is used, so beta rounding is not misattributed to its integrator.
std::pair<Vec3,Vec3> Analytic(const TransportInput& input,const double* field,
                                long double charge,long double beta) {
    const long double magnitude=std::sqrt(field[0]*field[0]+field[1]*field[1]+field[2]*field[2]);
    Vec3 position{},direction{},axis{},parallel{},perpendicular{},cross{};
    if(magnitude==0) {
        for(int i=0;i<3;++i) {
            direction[i]=input.direction[i];
            position[i]=input.position[i]+beta*input.duration*direction[i];
        }
        return {position,direction};
    }
    long double dot=0;
    for(int i=0;i<3;++i) { axis[i]=field[i]/magnitude; dot+=input.direction[i]*axis[i]; }
    for(int i=0;i<3;++i) { parallel[i]=dot*axis[i]; perpendicular[i]=input.direction[i]-parallel[i]; }
    for(int i=0;i<3;++i)
        cross[i]=input.direction[(i+1)%3]*axis[(i+2)%3]-input.direction[(i+2)%3]*axis[(i+1)%3];
    const long double omega=charge*units.Gauss*units.e*magnitude/input.energy;
    const auto angle=omega*input.duration;
    const auto sine=std::sin(angle),cosine=std::cos(angle);
    const auto oneMinusCosine=2*std::sin(angle/2)*std::sin(angle/2);
    for(int i=0;i<3;++i) {
        direction[i]=parallel[i]+perpendicular[i]*cosine+cross[i]*sine;
        position[i]=input.position[i]+beta*(parallel[i]*input.duration+
            perpendicular[i]*sine/omega+cross[i]*oneMinusCosine/omega);
    }
    return {position,direction};
}
void CheckAnalyticReference() {
    // Anchor the independent helix formula with simple known solutions before
    // using it to diagnose either numerical propagator. This is a floating-point
    // implementation check, not a trajectory-equivalence acceptance threshold.
    TransportInput input{};
    input.direction[0]=1;
    input.energy=units.Gauss*units.e;
    input.duration=std::acos(-1.L)/2;
    const double field[]={0,0,1}, zero[]={0,0,0};
    const auto tolerance=128*std::numeric_limits<double>::epsilon();
    for(const int charge : {-1,1}) {
        const auto solution=Analytic(input,field,charge,.5L);
        const Vec3 position={.5L,-charge*.5L,0}, direction={0,static_cast<long double>(-charge),0};
        const auto error=Difference(solution.first.data(),position.data())+
                         Difference(solution.second.data(),direction.data());
        if(!std::isfinite(error) || error>tolerance)
            throw std::runtime_error("Analytic quarter-turn check failed");
    }
    const auto free=Analytic(input,zero,1,.5L);
    const Vec3 position={.5L*input.duration,0,0};
    if(Difference(free.first.data(),position.data())>tolerance ||
       Difference(free.second.data(),input.direction)>tolerance)
        throw std::runtime_error("Analytic free-motion check failed");
}
void PrintArray(const long double* values,long double scale) {
    std::cout << '[';
    for(int i=0;i<3;++i) {
        if(!std::isfinite(values[i])) throw std::runtime_error("Nonfinite trajectory");
        std::cout << (i?",":"") << values[i]/scale;
    }
    std::cout << ']';
}
void PrintOutput(const TransportOutput& out,const TransportInput& input,const double* field,bool analytic) {
    if(!std::isfinite(out.elapsed) || !std::isfinite(out.beta)) throw std::runtime_error("Invalid time/speed");
    std::cout << "{\"position_Mpc\":"; PrintArray(out.position,units.Mpc);
    std::cout << ",\"direction\":"; PrintArray(out.direction,1);
    std::cout << ",\"flight_Mpc\":" << out.elapsed/units.Mpc << ",\"beta\":" << out.beta
              << ",\"mass_MeV\":" << out.mass*units.Eunit << ",\"charge\":" << out.charge
              << ",\"energy_eV\":" << out.energy/units.eV
              << ",\"initial_scale_Mpc\":" << out.initialScale/units.Mpc
              << ",\"global_scale_Mpc\":" << out.globalScale/units.Mpc
              << ",\"initial_rate_per_Mpc\":" << out.initialRate*units.Mpc
              << ",\"field_calls\":" << out.fieldCalls << ",\"scale_calls\":" << out.scaleCalls
              << ",\"zero_scale_calls\":" << out.zeroScaleCalls
              << ",\"direction_norm_error\":" << std::abs(Norm(out.direction)-1)
              << ",\"time_error_Mpc\":" << std::abs(out.elapsed-input.duration)/units.Mpc;
    if(analytic) {
        const auto reference=Analytic(input,field,out.charge,out.beta);
        std::cout << ",\"analytic_position_error_Mpc\":" << Difference(out.position,reference.first.data())/units.Mpc
                  << ",\"analytic_direction_error\":" << Difference(out.direction,reference.second.data());
    }
    std::cout << '}';
}
void Compare(const char* cases,const char* grid) {
    cosmology.Init(10.01,0.7,70);
    std::ifstream stream(cases);
    if(!stream) throw std::runtime_error("Cannot open frozen cases");
    std::string name;
    int kind,mode;
    TransportInput input{};
    double field[3],energyEV,durationMpc,scaleMpc;
    while(stream >> name >> kind >> mode >> input.particle >> input.accuracy >> input.outerSteps
                 >> energyEV >> field[0] >> field[1] >> field[2]
                 >> input.position[0] >> input.position[1] >> input.position[2]
                 >> input.direction[0] >> input.direction[1] >> input.direction[2] >> durationMpc >> scaleMpc) {
        if(input.particle<0 || input.particle>2 || input.accuracy<1 || input.accuracy>80 ||
           input.outerSteps<1 || input.outerSteps>64 || durationMpc<=0 || durationMpc>1 || scaleMpc<=0)
            throw std::runtime_error("Case outside bounded diagnostic design");
        input.energy=energyEV*units.eV;
        // The mildly relativistic control keeps the same input energy in both
        // programs; historical and modern beta formulas are left untouched.
        if(kind==2) {
            const ParticleType types[]={Electron,Positron,Proton};
            input.energy=10*Particle::Mass(types[input.particle]);
            const double amplitude=0.5*input.energy/(units.Gauss*units.e*units.Mpc);
            for(auto& value:field) value*=amplitude;
        }
        input.duration=durationMpc*units.Mpc;
        input.initialTime=cosmology.z2t0(0.01);
        const auto norm=Norm(input.direction);
        for(int i=0;i<3;++i) { input.position[i]*=units.Mpc; input.direction[i]/=norm; }
        const auto old=HistoricalTransport(input,kind==3?grid:nullptr,mode,field,scaleMpc*units.Mpc);
        const auto modern=Modern(input,kind==3?grid:nullptr,mode,field,scaleMpc*units.Mpc);
        if(old.charge!=modern.charge || old.mass!=modern.mass) throw std::runtime_error("Particle identity mismatch");
        if(old.energy!=input.energy || modern.energy!=input.energy) throw std::runtime_error("Unexpected energy change");
        std::cout << "{\"case\":\"" << name << "\",\"kind\":" << kind << ",\"mode\":" << mode
                  << ",\"particle\":" << input.particle << ",\"accuracy\":" << input.accuracy
                  << ",\"outer_steps\":" << input.outerSteps << ",\"energy_eV\":" << input.energy/units.eV
                  << ",\"field_G\":[" << field[0] << ',' << field[1] << ',' << field[2] << ']'
                  << ",\"duration_Mpc\":" << durationMpc << ",\"initial_time_Mpc\":" << input.initialTime/units.Mpc
                  << ",\"initial_position_Mpc\":"; PrintArray(input.position,units.Mpc);
        std::cout << ",\"initial_direction\":"; PrintArray(input.direction,1);
        std::cout << ",\"fork\":"; PrintOutput(old,input,field,kind!=3);
        std::cout << ",\"modern\":"; PrintOutput(modern,input,field,kind!=3);
        std::cout << ",\"difference_position_Mpc\":" << Difference(old.position,modern.position)/units.Mpc
                  << ",\"difference_direction\":" << Difference(old.direction,modern.direction)
                  << ",\"difference_beta\":" << std::abs(old.beta-modern.beta) << "}\n";
    }
    if(!stream.eof()) throw std::runtime_error("Malformed frozen case");
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=3) throw std::runtime_error("Usage: fork-transport-probe CASES GRID");
        std::cout << std::setprecision(21);
        CheckAnalyticReference();
        Compare(argv[1],argv[2]);
        return 0;
    } catch(Utils::Exception* error) {
        std::unique_ptr<Utils::Exception> owner(error); std::cerr << error->Message() << '\n';
    } catch(const std::exception& error) { std::cerr << error.what() << '\n'; }
    return 1;
}
