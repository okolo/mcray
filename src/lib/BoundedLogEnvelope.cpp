#include "BoundedLogEnvelope.h"

namespace Utils {

BoundedLogEnvelope::BoundedLogEnvelope(const Function& density,double xmin,double xmax,double observedMaximum):
        density(density),xmin(xmin),xmax(xmax),logmin(std::log(xmin)),logmax(std::log(xmax)),maximum(observedMaximum) {
    if(!density.HasIntervalUpperBound())return;
    std::vector<Cell> storage;storage.reserve(32);
    cells=std::priority_queue<Cell>(std::less<Cell>(),std::move(storage));
    if(!(xmin>0 && xmax>xmin && std::isfinite(logmin) && std::isfinite(logmax)))
        throw new Exception("Invalid bounded sampler support",false);
    if(!(logmax>logmin))return;
    cells.push(cell(point(logmin),point(logmax)));
    // Refine the interval with the largest upper bound until the global M
    // is within 1% of an observed value. This is an efficiency target: it
    // does NOT permit g to exceed M by 1%, nor a 1% error in probabilities.
    // The budget limits construction time; failure selects CDF in the caller.
    for(int iteration=0;iteration<4096;++iteration) {
        if(cells.top().upper<=1.01*maximum){supported=true;return;}
        const auto parent=cells.top();cells.pop();
        const double t=(parent.lo.t+parent.hi.t)/2;
        if(t==parent.lo.t || t==parent.hi.t)return;
        const auto middle=point(t);
        cells.push(cell(parent.lo,middle));cells.push(cell(middle,parent.hi));
    }
    // No random numbers have been consumed. The caller may select CDF
    // before starting an event if an enclosure could not be tightened.
}

BoundedLogEnvelope::Point BoundedLogEnvelope::Proposal(double uniform) const {
    if(!(uniform>=0 && uniform<=1))throw new Exception("Invalid rejection proposal variate",false);
    return point(logmin+(logmax-logmin)*uniform);
}

bool BoundedLogEnvelope::Accept(const Point& p,double uniform) const {
    if(!(uniform>=0 && uniform<=1))throw new Exception("Invalid rejection acceptance variate",false);
    if(!supported || p.y>Upper())throw new Exception("Analytic rejection envelope violated",false);
    return p.y>0 && uniform*Upper()<=p.y;
}

double BoundedLogEnvelope::Draw(mcray::Randomizer& rng,size_t* attempts) const {
    for(int i=0;i<1000000;++i) {
        const auto p=Proposal(rng.Rand());const double u=rng.Rand();
        if(Accept(p,u)){if(attempts)*attempts=i+1;return p.x;}
    }
    // An unusually low acceptance rate is an error, not permission to
    // narrow the physical support or return the location of a sampled peak.
    throw new Exception("Bounded rejection attempt limit exceeded",false);
}

BoundedLogEnvelope::Point BoundedLogEnvelope::point(double t) const {
    const double x=t==logmin?xmin:t==logmax?xmax:std::clamp(std::exp(t),xmin,xmax);
    const double y=x*density(x);
    if(!std::isfinite(y)||y<0)throw new Exception("Invalid bounded sampling density",false);
    ++evaluations;maximum=std::max(maximum,y);return {t,x,y};
}

BoundedLogEnvelope::Cell BoundedLogEnvelope::cell(Point lo,Point hi) const {
    const double upper=density.LogDensityUpperBound(lo.x,hi.x,lo.y,hi.y,1.01*maximum);
    if(!std::isfinite(upper)||upper<0||upper<std::max(lo.y,hi.y))
        throw new Exception("Invalid analytic log-density enclosure",false);
    return {lo,hi,upper};
}

} // namespace Utils
