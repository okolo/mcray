#ifndef MCRAY_SECONDARY_ORACLE_H
#define MCRAY_SECONDARY_ORACLE_H
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <gsl/gsl_integration.h>
namespace secondary {
using Real = long double;
struct Audit {
    unsigned long long proposals=0, envelope=0, narrowing=0, argmax=0;
    double largestEnvelopeRatio=1.;
};
extern Audit audit;
// Direct differential laws: arXiv:1106.5508v2, equations (11), (12).
// Overall factors cancel in a distribution conditional on s.
inline Real pairDensity(Real r, Real s) {
    const Real u=1-r, a=4/s;
    return (r*r/u+u+a/u-a*a/(4*r*u*u))/r/(1+2*(1-a)*a);
}
inline Real comptonDensity(Real y, Real s) {
    const Real a=1/s, b=1-a;
    return ((1+y*y)/2-2*a*(y-a)*(1-y)/(y*b*b))/y;
}
inline Real pairMinimum(Real s) {return 2/(s*(1+std::sqrt(1-4/s)));}
struct Region {
    std::string id,view,kind;
    Real lo,hi;
    int multiplicity(Real first, Real second, bool firstElectron) const {
        auto inside=[&](Real x){return x>=lo && (x<hi || (hi==1 && x==1));};
        if(view=="gamma" || view=="low")return inside(first);
        if(view=="recoil" || view=="high")return inside(second);
        if(view=="electron")return inside(firstElectron?first:second);
        if(view=="positron")return inside(firstElectron?second:first);
        if(view=="inclusive")return inside(first)+inside(second);
        throw std::runtime_error("Unknown secondary observable");
    }
};
inline std::vector<Region> regions(bool pair,const std::vector<Real>& edges) {
    std::vector<Region> result;
    const std::vector<std::string> views=pair?std::vector<std::string>{"low","high","electron","positron","inclusive"}:std::vector<std::string>{"gamma","recoil"};
    for(const auto& view:views) {
        for(size_t i=1;i<edges.size();++i)result.push_back({view+".bin."+std::to_string(i-1),view,"bin",edges[i-1],edges[i]});
        for(Real cut:{1e-2L,1e-3L,1e-4L,1e-6L,1e-8L}) {
            const auto name=std::to_string(double(std::log10(cut)));
            result.push_back({view+".lower."+name,view,"tail",0,cut});
            result.push_back({view+".upper."+name,view,"tail",1-cut,1});
        }
        if(pair && (view=="low" || view=="high"))for(Real d:{1e-2L,1e-4L,1e-6L}) {
            if(view=="low")result.push_back({view+".endpoint."+std::to_string(double(d)),view,"tail",.5L*(1-d),.5L});
            else result.push_back({view+".endpoint."+std::to_string(double(d)),view,"tail",.5L,.5L*(1+d)});
        }
    }
    return result;
}
struct Quadrature {
    gsl_integration_glfixed_table* table;
    explicit Quadrature(int order):table(gsl_integration_glfixed_table_alloc(order)) {if(!table)throw std::runtime_error("Cannot allocate quadrature");}
    ~Quadrature(){gsl_integration_glfixed_table_free(table);}
    template<class F> Real logarithmic(F density,Real lo,Real hi) const {
        if(hi<=lo)return 0;
        if(!(lo>0))throw std::runtime_error("Nonpositive quadrature bound");
        const Real a=std::log(lo),b=std::log(hi);
        const int pieces=std::max(1,int(std::ceil((b-a)/.25L)));
        Real value=0;
        for(int j=0;j<pieces;++j)for(size_t i=0;i<table->n;++i) {
            double t,w;gsl_integration_glfixed_point(double(a+(b-a)*j/pieces),double(a+(b-a)*(j+1)/pieces),i,&t,&w,table);
            const Real x=std::exp(Real(t)),f=density(x);
            if(!std::isfinite(f)||f<0)throw std::runtime_error("Invalid independent differential density");
            value+=w*x*f;
        }
        return value;
    }
};
// Directly integrate disjoint event-pattern intervals. No difference of CDFs
// near one, no production primitive, root solver or secondary sampler.
inline std::vector<std::array<Real,3>> conditional(bool pair,Real s,Real threshold,
        const std::vector<Region>& obs,const Quadrature& q, bool supportOnly=false) {
    const Real lo=pair?pairMinimum(s):1/s;
    const Real hi=pair?.5L:1-threshold;
    std::vector<std::array<Real,3>> mass(obs.size());
    if(hi<=lo)return mass;
    std::vector<Real> breaks{lo,hi};
    for(const auto& o:obs)for(Real x:{o.lo,o.hi,1-o.lo,1-o.hi})if(x>lo && x<hi)breaks.push_back(x);
    std::sort(breaks.begin(),breaks.end());breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
    Real total=0;
    for(size_t k=1;k<breaks.size();++k) {
        const Real left=breaks[k-1],right=breaks[k],middle=(left+right)/2;
        const Real weight=supportOnly?1:q.logarithmic([&](Real x){return pair?pairDensity(x,s):comptonDensity(x,s);},left,right);
        total+=weight;
        const Real first=pair?middle:1-middle,second=pair?1-middle:middle;
        for(size_t j=0;j<obs.size();++j) {
            if(pair)for(bool electron:{false,true})mass[j][obs[j].multiplicity(first,second,electron)]+=weight/2;
            else mass[j][obs[j].multiplicity(first,second,true)]+=weight;
        }
    }
    if(!(total>0))throw std::runtime_error("Empty independent conditional density");
    for(auto& row:mass)for(auto& x:row)x/=total;
    return mass;
}
}
#endif
