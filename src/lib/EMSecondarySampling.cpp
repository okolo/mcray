#include "EMSecondarySampling.h"
#include "Utils.h"
#include <algorithm>
#include <cmath>

namespace Interactions { namespace detail {

long double SecondaryLog1p(long double value) {
    if(std::abs(value)>=.01L)return std::log1p(value);
    // Around zero, sum log(1+value) directly to avoid cancellation and the
    // inaccurate x87 log1pl implementation observed under Linux emulation.
    // |value|<0.01 makes the omitted terms after order 18 negligible compared
    // with the tail tolerances used by the inverse search.
    long double term=value,sum=value;
    for(int n=2;n<=18;++n) {term*=-value;sum+=term/n;}
    return sum;
}
// ELMAG Eq. (9). x and y are the lower and upper recoil
// energy fractions E_e/E_initial, not the emitted photon fractions. The upper
// photon cut therefore becomes a recoil upper bound 1-rMin.
// The small-interval expansion is valid when q=(y-x)/x is small. Using only
// the absolute width y-x, as the former implementation did, fails at high S
// where x=m_e^2/S is tiny, and can make the normalized CDF nonmonotone.
long double ComptonIntegral(long double x, long double y) {
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

// Pair-production kinematics: r_min=(1-sqrt(1-4/s))/2, with s=S/m_e^2.
// Rationalizing this expression avoids subtracting nearly equal numbers at
// high invariant energy. PairIntegral integrates ELMAG Eq. (11) from r_min;
// differences of logarithms use log1p so a very narrow tail remains resolvable.
long double PairMinimum(long double s) {
    return 2/(s*(1+std::sqrt(1-4/s)));
}
long double PairIntegral(long double r,long double s) {
    const long double r0=PairMinimum(s), delta=r-r0;
    if(delta<=0)return 0;
    const long double a=4/s, u=1-r, u0=1-r0;
    const long double logarithm=SecondaryLog1p(delta/r0)-SecondaryLog1p(-delta/u0);
    return ((1+a-a*a/2)*logarithm-2*delta
        -a*a/4*delta*(1/(r*r0)+1/(u*u0)))/(1+2*a*(1-a));
}

// This inverts the conditional secondary-energy CDF at a GIVEN invariant S.
// It does not select S or replace the fast background rejection sampler.
// Search log(energy fraction) to resolve tiny recoil/soft-lepton fractions.
// The probability tolerance scales with min(u,1-u), keeping both probability
// tails resolved instead of stopping when an absolute CDF residual looks small.
// A nonfinite/out-of-range CDF or exhaustion is a non-retryable calculation
// error: retrying only failed histories would bias the retained cascade sample.
long double InvertSecondaryCDF(const SecondaryCDF& cdf,long double lo,
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

}} // namespace Interactions::detail
