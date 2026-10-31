// SPDX-License-Identifier: MIT
#include "Jet.h"
#include <algorithm>
#include <cmath>
#include <limits>
namespace mcray {
namespace {
[[noreturn]] void JetError(const char* text) {
    Utils::Exception::ThrowFatal(std::string("Jet: ")+text);
    std::terminate();
}
double Weight(double scal,double kappa) {
    if(!std::isfinite(scal) || !std::isfinite(kappa)) JetError("nonfinite vMF input");
    // Preserve the original contract, including its discontinuity at zero.
    if(std::abs(kappa)<1.e-14) return 1.;
    if(std::abs(scal)>1+64*std::numeric_limits<double>::epsilon())
        JetError("vMF scalar product is outside the unit-vector range");
    scal=std::clamp(scal,-1.,1.); // Unit-vector dot products can round beyond a pole.
    // Negative kappa uses the same original density around the opposite pole.
    // This algebraic rewrite avoids exp(+2*abs(kappa)) overflow.
    if(kappa<0) { kappa=-kappa; scal=-scal; }
    const double offset=1-scal;
    if(offset>1 && kappa>std::numeric_limits<double>::max()/offset) return 0.;
    const double product=kappa*offset;
    const double mass=kappa>std::numeric_limits<double>::max()/2 ? 1. : -std::expm1(-2*kappa);
    const double exponential=std::exp(-product);
    const double pi=std::acos(-1.);
    double result=kappa*exponential/(2*pi*mass);
    // An underflowed exponential can still have a representable weighted
    // density after multiplication by kappa. Preserve that tail in log space.
    if(exponential==0 && std::isfinite(product))
        result=std::exp(std::log(kappa)-product-std::log(2*pi*mass));
    if(!std::isfinite(result) || result<0) JetError("unrepresentable vMF multiplier");
    return result;
}
}
double vMF_weight(double setup_theta,double setup_phi,const coord_type* pdir,double kappa) {
    // Keep the fork's angle-based entry point and axis construction.
    const double Px=std::cos(setup_phi)*std::sin(setup_theta);
    const double Py=std::sin(setup_phi)*std::sin(setup_theta);
    const double Pz=std::cos(setup_theta);
    const double scal=Px*pdir[0]+Py*pdir[1]+Pz*pdir[2];
    return Weight(scal,kappa);
}
double vMF_weight(const coord_type* axis,const coord_type* pdir,double kappa) {
    const double scal=axis[0]*pdir[0]+axis[1]*pdir[1]+axis[2]*pdir[2];
    return Weight(scal,kappa);
}
}
