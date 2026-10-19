#ifndef MCRAY_INVARIANT_CDF_DIAGNOSTIC_H
#define MCRAY_INVARIANT_CDF_DIAGNOSTIC_H
#include "NumericalLogCDF.h"
#include <chrono>
#include <fstream>
#include <iomanip>

bool LegacyRejectionSampleLogDistribution(const Utils::Function&,mcray::Randomizer&,double&,double&,double,double,double,int=3);
bool HistoricalSampleLogDistribution(const Utils::Function&,double,double&,double&,double,double,double);

namespace secondary {
// Fixed composite Gauss-Legendre oracle, with separately accumulated prefix
// and suffix masses. It does not use the production adaptive integrator.
struct InvariantReference {
    const Utils::Function& density;
    double lo,step;
    std::vector<long double> prefix,suffix;
    Quadrature quadrature;
    InvariantReference(const Utils::Function& f,double xmin,double xmax,int cells,int order):
        density(f),lo(std::log(xmin)),step((std::log(xmax)-lo)/cells),prefix(cells+1),suffix(cells+1),quadrature(order) {
        for (int i=0;i<cells;++i) prefix[i+1]=prefix[i]+mass(lo+i*step,lo+(i+1)*step);
        for (int i=cells-1;i>=0;--i) suffix[i]=suffix[i+1]+mass(lo+i*step,lo+(i+1)*step);
    }
    long double mass(double a,double b) const {
        long double sum=0;
        for (size_t i=0;i<quadrature.table->n;++i) {
            double t,w;gsl_integration_glfixed_point(a,b,i,&t,&w,quadrature.table);
            const double x=std::exp(t);sum+=static_cast<long double>(w)*x*density(x);
        }
        return sum;
    }
    long double probability(double x,bool right) const {
        const double t=std::log(x);
        const int i=std::clamp(int((t-lo)/step),0,int(prefix.size())-2);
        return (right?suffix[i+1]+mass(t,lo+(i+1)*step):prefix[i]+mass(lo+i*step,t))/prefix.back();
    }
};

inline void invariantDiagnostic(const Utils::Function& f,double lo,double hi,double tol,const std::string& output) {
    InvariantReference coarse(f,lo,hi,8192,16),fine(f,lo,hi,16384,24);
    Utils::NumericalLogCDF cdf(f,lo,hi,tol);
    std::ofstream out(output);
    out<<std::setprecision(21)<<"# total="<<cdf.IntegralValue()<<" total_error_estimate="<<cdf.IntegralError()
       <<" reference="<<fine.prefix.back()<<" reference_coarse="<<coarse.prefix.back()<<"\n";
    out<<"u\tx\tsmall_tail_reference\trelative_probability_error\treference_change_relative\tevaluations\titerations\n";
    for (double u:{1e-12,1e-10,1e-8,1e-6,1e-4,.001,.01,.1,.25,.5,.75,.9,.99,.999,1-1e-4,1-1e-6,1-1e-8,1-1e-10,1-1e-12}) {
        const auto count=cdf.Evaluations(),iterations=cdf.Iterations();
        const double x=cdf.Quantile(u),p=u>.5?1-u:u;
        const long double actual=fine.probability(x,u>.5),check=coarse.probability(x,u>.5);
        out<<u<<'\t'<<x<<'\t'<<actual<<'\t'<<actual/p-1<<'\t'<<(actual-check)/p<<'\t'<<cdf.Evaluations()-count<<'\t'<<cdf.Iterations()-iterations<<'\n';
    }
    // Full construction and inversion for every draw: no cached CDF shortcut.
    // Preserve repeated block timings for a fair comparison with the old code.
    out<<"# benchmark: sampler repeat N seconds checksum\n";
    for (int repeat=0;repeat<3;++repeat) for (int offset=0;offset<3;++offset) {
        const int mode=(repeat+offset)%3;
        mcray::Randomizer rng(20261020+repeat);Utils::MathUtils math;
        double sum=0,integral=0,x=0;
        const auto start=std::chrono::steady_clock::now();
        for (int i=0;i<128;++i) {
            const bool ok=mode==0?math.SampleLogDistributionCDF(f,rng,x,integral,lo,hi,tol):
                mode==1?LegacyRejectionSampleLogDistribution(f,rng,x,integral,lo,hi,tol):
                HistoricalSampleLogDistribution(f,rng.Rand(),x,integral,lo,hi,tol);
            if (!ok) throw std::runtime_error("Benchmark sampling failed");
            sum+=std::log(x);
        }
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        out<<"# benchmark "<<(mode==0?"cdf":mode==1?"rejection":"nr")<<' '<<repeat<<" 128 "<<seconds<<' '<<sum<<'\n';
    }
    if (!out) throw std::runtime_error("Cannot write invariant CDF diagnostics");
}
}
#endif
