// Fixed-s diagnostic; production sources are linked without modification.
#include "ICS.h"
#include "MathUtils.h"
#include "secondary_oracle.h"
#include <gsl/gsl_errno.h>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace Interactions {
// Existing diagnostic primitives in production ICS.cpp, without prefactors.
double PartialSigmaSeries(double r, double s);
double PartialSigmaExact(double r, double s);
}
using namespace Interactions;
using namespace mcray;
using secondary::Real;

namespace {
constexpr double rmin=1e-8;
struct Equation {double s,u; bool exact;};
double primitive(double r,double s,bool exact) {
    if(1-r<=1/s)return 0;
    return exact || 1-r-1/s>1e-3 ? PartialSigmaExact(r,s) : PartialSigmaSeries(r,s);
}
double equation(double r,void* opaque) {
    const auto& p=*static_cast<Equation*>(opaque);
    return primitive(r,p.s,p.exact)/primitive(rmin,p.s,p.exact)-p.u;
}
double invert(double s,double u,bool exact,double tolerance) {
    Equation p{s,u,exact};
    return Utils::MathUtils::SolveEquation(equation,rmin,1-1/s,&p,tolerance,100,gsl_root_fsolver_bisection);
}
Real cdf(Real y,Real s,const secondary::Quadrature& q) {
    auto density=[&](Real t){return secondary::comptonDensity(t,s);};
    return q.logarithmic(density,1/s,y)/q.logarithmic(density,1/s,1-rmin);
}
Real referenceY(Real s,Real u,const secondary::Quadrature& q) {
    Real lo=std::log(1/s),hi=std::log(1-rmin);
    for(int i=0;i<65;++i) {
        const Real mid=(lo+hi)/2;
        if(cdf(std::exp(mid),s,q)<u)lo=mid;else hi=mid;
    }
    return std::exp((lo+hi)/2);
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=2)throw std::runtime_error("Usage: ics_endpoint_diagnostic OUTPUT_PREFIX");
        gsl_set_error_handler_off();
        secondary::Quadrature q(32),check(16);
        std::ofstream primitives(std::string(argv[1])+"-primitives.tsv");
        std::ofstream quantiles(std::string(argv[1])+"-quantiles.tsv");
        primitives<<std::setprecision(18)<<"S\ty\ta_over_ymin\tproduction_cdf\tdirect_cdf\texact_cdf\n";
        quantiles<<std::setprecision(18)<<"S\tseed\tu\ty_current\ty_tight_only\ty_exact_only\ty_both\ty_reference\tcdf_current\tcdf_tight_only\tcdf_exact_only\tcdf_both\n";
        const double mass=Particle::Mass(Electron);
        double maximumReplicaDifference=0,maximumQuadratureDifference=0,maximumCurrentCDFResidual=0;
        for(double s:{10.,100.,1e3,1e4,1e6,1e8}) {
            for(double y:{1.0001/s,2/s,10/s,1e-4,5e-4,9e-4,1.1e-3,.01,.1}) {
                if(y<=1/s || y>=1-rmin)continue;
                const Real direct=cdf(y,s,q);
                maximumQuadratureDifference=std::max(maximumQuadratureDifference,double(std::abs(direct-cdf(y,s,check))));
                primitives<<s<<'\t'<<y<<'\t'<<(y-1/s)*s<<'\t'
                    <<primitive(1-y,s,false)/primitive(rmin,s,false)<<'\t'<<direct<<'\t'
                    <<primitive(1-y,s,true)/primitive(rmin,s,true)<<'\n';
            }
            Particle primary(Electron,0.);primary.Energy=1e20*units.eV;
            ICS ics(nullptr,rmin*primary.Energy);
            for(unsigned seed=1;seed<=12;++seed) {
                Randomizer draw(seed),sample(seed);
                const double u=draw.Rand();
                std::vector<Particle> products;
                ics.SampleSecondaries(primary,products,s*mass*mass,sample);
                const double production=products.at(1).Energy/primary.Energy;
                maximumCurrentCDFResidual=std::max(maximumCurrentCDFResidual,double(std::abs(cdf(production,s,q)-u)));
                const double replica=1-invert(s,u,false,1e-3);
                maximumReplicaDifference=std::max(maximumReplicaDifference,std::abs(production-replica));
                const double tight=1-invert(s,u,false,1e-12),exact=1-invert(s,u,true,1e-3),both=1-invert(s,u,true,1e-12);
                quantiles<<s<<'\t'<<seed<<'\t'<<u<<'\t'<<production<<'\t'<<tight<<'\t'<<exact<<'\t'<<both<<'\t'<<referenceY(s,u,q)
                    <<'\t'<<cdf(production,s,q)<<'\t'<<cdf(tight,s,q)<<'\t'<<cdf(exact,s,q)<<'\t'<<cdf(both,s,q)<<'\n';
            }
        }
        if(!primitives || !quantiles)throw std::runtime_error("Cannot save diagnostics");
        std::cout<<std::setprecision(18)<<"maximum_legacy_replica_absolute_energy_difference="<<maximumReplicaDifference
            <<"\nmaximum_current_cdf_residual="<<maximumCurrentCDFResidual
            <<"\nmaximum_16_32_quadrature_cdf_difference="<<maximumQuadratureDifference<<'\n';
        if(maximumCurrentCDFResidual>1e-9 || maximumQuadratureDifference>1e-12)
            throw std::runtime_error("Diagnostic control failed");
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
