#ifndef MCRAY_EM_SECONDARY_SAMPLING_H
#define MCRAY_EM_SECONDARY_SAMPLING_H
#include <algorithm>
#include <cmath>
#include "Utils.h"

namespace Interactions { namespace detail {
inline long double SecondaryLog1p(long double value) {
    if(std::abs(value)>=.01L)return std::log1p(value);
    // Also avoids inaccurate x87 log1pl implementations under emulation.
    long double term=value,sum=value;
    for(int n=2;n<=18;++n) {term*=-value;sum+=term/n;}
    return sum;
}
// ELMAG Eq. (9), with the cancellation-prone logarithmic remainder evaluated
// by its convergent series. Arguments are recoil fractions, not photon ones.
inline long double ComptonIntegral(long double x, long double y) {
    if(y<=x)return 0;
    const long double b=1-x, delta=y-x, q=delta/x;
    const long double logarithm=SecondaryLog1p(q);
    long double remainder;
    if(q<.01L) {
        remainder=-b*q*q/2;
        long double power=q*q;
        for(int n=3;n<=18;++n) {
            power*=q;
            remainder+=(n%2?1:-1)*power*(n-1-x)/n;
        }
    } else {
        remainder=q/(1+q)+delta-(1+x)*logarithm;
    }
    return x/b*(logarithm+delta*(y+x)/2+4*x/(b*b)*remainder);
}

// Integral of ELMAG Eq. (11) from the stable pair-production lower bound.
inline long double PairMinimum(long double s) {
    return 2/(s*(1+std::sqrt(1-4/s)));
}
inline long double PairIntegral(long double r,long double s) {
    const long double r0=PairMinimum(s), delta=r-r0;
    if(delta<=0)return 0;
    const long double a=4/s, u=1-r, u0=1-r0;
    const long double logarithm=SecondaryLog1p(delta/r0)-SecondaryLog1p(-delta/u0);
    return ((1+a-a*a/2)*logarithm-2*delta
        -a*a/4*delta*(1/(r*r0)+1/(u*u0)))/(1+2*a*(1-a));
}

// Invert in a logarithmic energy coordinate, resolving probabilities relative
// to BOTH tails. Do not change the shared MathUtils root-solver contract.
template<class CDF> long double InvertSecondaryCDF(CDF cdf,long double lo,
        long double hi,long double uniform) {
    if(uniform<=0)return lo;
    if(uniform>=1)return hi;
    const long double tolerance=1e-12L*std::min(uniform,1-uniform);
    for(int iteration=0;iteration<96;++iteration) {
        const long double middle=(lo+hi)/2;
        if(middle==lo || middle==hi)return middle;
        const long double probability=cdf(middle);
        if(!std::isfinite(probability) || probability< -1e-14L || probability>1+1e-14L)
            throw new Utils::Exception("Invalid secondary sampling CDF",false);
        if(std::abs(probability-uniform)<=tolerance)return middle;
        if(probability<uniform)lo=middle;else hi=middle;
    }
    throw new Utils::Exception("Secondary inverse CDF did not converge",false);
}
}}
#endif
