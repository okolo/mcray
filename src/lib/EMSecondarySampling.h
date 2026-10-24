#ifndef MCRAY_EM_SECONDARY_SAMPLING_H
#define MCRAY_EM_SECONDARY_SAMPLING_H
#include <functional>

namespace Interactions { namespace detail {
// Stable conditional secondary-energy integrals. Dimensionless s=S/m_e^2;
// Compton arguments are recoil fractions; pair arguments are lepton fractions.
long double ComptonIntegral(long double x,long double y);
long double PairMinimum(long double s);
long double PairIntegral(long double r,long double s);

using SecondaryCDF = std::function<long double(long double)>;
// lo/hi are logarithmic energy coordinates; uniform is the CDF probability.
// The implementation is independent of background invariant-S sampling.
long double InvertSecondaryCDF(const SecondaryCDF& cdf,long double lo,
                              long double hi,long double uniform);
}}
#endif
