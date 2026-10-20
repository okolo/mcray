#ifndef MCRAY_BOUNDED_LOG_ENVELOPE_H
#define MCRAY_BOUNDED_LOG_ENVELOPE_H
#include "MathUtils.h"
#include <queue>

namespace Utils {
// Guard arithmetic and elementary-function roundoff separately from the 1%
// envelope construction budget. Analytic inequalities, not a sampled maximum,
// supply the enclosure. This is not an interval-arithmetic proof of libm.
inline double SamplingUpperRound(long double value) {
    if(value<0 || !std::isfinite(value))throw new Exception("Invalid analytic sampling bound",false);
    if(value==0)return 0;
    return std::nextafter(static_cast<double>(value*(1+1e-10L)),std::numeric_limits<double>::infinity());
}

class BoundedLogEnvelope {
public:
    struct Point {double t,x,y;};
    struct Cell {
        Point lo,hi;double upper;
        bool operator<(const Cell& other) const {return upper<other.upper;}
    };
    BoundedLogEnvelope(const Function& density,double xmin,double xmax,double observedMaximum=0):
        density(density),xmin(xmin),xmax(xmax),logmin(std::log(xmin)),logmax(std::log(xmax)),maximum(observedMaximum) {
        if(!density.HasIntervalUpperBound())return;
        std::vector<Cell> storage;storage.reserve(32);
        cells=std::priority_queue<Cell>(std::less<Cell>(),std::move(storage));
        if(!(xmin>0 && xmax>xmin && std::isfinite(logmin) && std::isfinite(logmax)))
            throw new Exception("Invalid bounded sampler support",false);
        cells.push(cell(point(logmin),point(logmax)));
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
    bool Supported() const {return supported;}
    double Upper() const {return cells.empty()?0:cells.top().upper;}
    double ObservedMaximum() const {return maximum;}
    size_t Evaluations() const {return evaluations;}
    size_t Cells() const {return cells.size();}
    const std::priority_queue<Cell>& Partition() const {return cells;}
    Point Proposal(double uniform) const {
        if(!(uniform>=0 && uniform<=1))throw new Exception("Invalid rejection proposal variate",false);
        return point(logmin+(logmax-logmin)*uniform);
    }
    bool Accept(const Point& p,double uniform) const {
        if(!(uniform>=0 && uniform<=1))throw new Exception("Invalid rejection acceptance variate",false);
        if(!supported || p.y>Upper())throw new Exception("Analytic rejection envelope violated",false);
        return p.y>0 && uniform*Upper()<=p.y;
    }
    double Draw(mcray::Randomizer& rng,size_t* attempts=nullptr) const {
        for(int i=0;i<1000000;++i) {
            const auto p=Proposal(rng.Rand());const double u=rng.Rand();
            if(Accept(p,u)){if(attempts)*attempts=i+1;return p.x;}
        }
        throw new Exception("Bounded rejection attempt limit exceeded",false);
    }
private:
    Point point(double t) const {
        const double x=std::clamp(std::exp(t),xmin,xmax),y=x*density(x);
        if(!std::isfinite(y)||y<0)throw new Exception("Invalid bounded sampling density",false);
        ++evaluations;maximum=std::max(maximum,y);return {t,x,y};
    }
    Cell cell(Point lo,Point hi) const {
        const double upper=density.LogDensityUpperBound(lo.x,hi.x,lo.y,hi.y,1.01*maximum);
        if(!std::isfinite(upper)||upper<0||upper<std::max(lo.y,hi.y))
            throw new Exception("Invalid analytic log-density enclosure",false);
        return {lo,hi,upper};
    }
    const Function& density;
    double xmin,xmax,logmin,logmax;
    mutable double maximum=0;
    mutable size_t evaluations=0;
    bool supported=false;
    std::priority_queue<Cell> cells;
};
}
#endif
