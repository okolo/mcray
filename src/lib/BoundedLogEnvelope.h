#ifndef MCRAY_BOUNDED_LOG_ENVELOPE_H
#define MCRAY_BOUNDED_LOG_ENVELOPE_H
#include "MathUtils.h"
#include <queue>

namespace Utils {
// Convert an analytically derived upper bound to the double precision used by
// the rejection test. A rounded-down bound could accept too many proposals near
// a density peak and thereby bias the distribution, even if its integral is
// accurate. Inflate a positive bound by 1e-10, then advance one double towards
// +infinity so that the final conversion cannot round the inflated value down.
// Keep exact zero: an empty interval must not acquire artificial probability.
//
// The 1e-10 factor is an engineering allowance for evaluation roundoff, not a
// proven error bound for all intermediate operations, cross-section evaluations
// or implementations of exp/log. Long-double intermediates reduce roundoff but
// do not provide directed rounding for those operations. Independent envelope
// audits and the fatal g > M check in Accept remain necessary safeguards.
//
// This allowance has a different purpose from the constructor's 1% threshold:
// that threshold limits extra rejection work, not error in the output spectrum.
inline double SamplingUpperRound(long double value) {
    if(value<0 || !std::isfinite(value))throw new Exception("Invalid analytic sampling bound",false);
    if(value==0)return 0;
    return std::nextafter(static_cast<double>(value*(1+1e-10L)),std::numeric_limits<double>::infinity());
}

// Rejection sampling with a constant envelope for the density in log(x).
// If density(x)=f(x), then t=log(x) has density g(t)=x*f(x), including the
// Jacobian x. For invariant sampling x is S; the density need not be normalized.
// Uniform t proposals are accepted with probability g(t)/M. Any fixed M above
// g everywhere gives the intended distribution; excess M only costs attempts.
//
// Construction partitions the support and asks the density for an analytic
// upper bound on EACH WHOLE interval. Evaluated points give a lower estimate of
// max(g), used to decide when refinement is sufficient; they never establish M.
// The partition is an envelope-construction aid, not a binned approximation to
// the sampled spectrum. Draw still evaluates the original density at proposals.
// See docs/validation/bounded_rejection_sampling.md for the bound derivations.
class BoundedLogEnvelope {
public:
    struct Point {double t,x,y;}; // log(x), physical x, and x*f(x), respectively.
    struct Cell {
        Point lo,hi;double upper;
        bool operator<(const Cell& other) const {return upper<other.upper;}
    };
    BoundedLogEnvelope(const Function& density,double xmin,double xmax,double observedMaximum=0);
    bool Supported() const {return supported;}
    double Upper() const {return cells.empty()?0:cells.top().upper;}
    double ObservedMaximum() const {return maximum;}
    size_t Evaluations() const {return evaluations;}
    size_t Cells() const {return cells.size();}
    const std::priority_queue<Cell>& Partition() const {return cells;}
    Point Proposal(double uniform) const;
    // Keep M fixed once proposals start. Enlarging an underestimated M after
    // encountering g > M cannot undo the bias of earlier accepted proposals.
    // A violation must stop the calculation, rather than discard one history.
    bool Accept(const Point& p,double uniform) const;
    double Draw(mcray::Randomizer& rng,size_t* attempts=nullptr) const;
private:
    Point point(double t) const;
    // Endpoint y values are g=x*f(x), not f(x). The last argument is a hint:
    // a bound already below that value needs no further expensive tightening.
    // The density must still return a valid whole-interval bound in that case.
    Cell cell(Point lo,Point hi) const;
    const Function& density;
    double xmin,xmax,logmin,logmax;
    mutable double maximum=0;
    mutable size_t evaluations=0;
    bool supported=false;
    std::priority_queue<Cell> cells;
};
}
#endif
