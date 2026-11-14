#include "coupled_modern.h"
#include "SpatialGeometry.h"
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
using namespace mcray;
using namespace Interactions;
namespace {
long double Difference(const long double* a,const long double* b) {
    return std::hypot(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
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
        auto modernField=ModernCoupledField(fieldKind,grid,field,.1*units.Mpc);
        const auto modern=ModernCoupled(c,*modernField);
        std::cout << "{\"profile\":\"" << profile << "\",\"ray\":\"" << ray << "\",\"field\":" << fieldKind
                  << ",\"particle\":" << particle << ",\"accuracy\":" << accuracy << ",\"fork\":";
        PrintCoupled(old); std::cout << ",\"modern\":"; PrintCoupled(modern);
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
