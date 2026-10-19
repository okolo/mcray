#ifndef MCRAY_NUMERICAL_LOG_CDF_H
#define MCRAY_NUMERICAL_LOG_CDF_H

#include "MathUtils.h"
#include "Utils.h"
#include <gsl/gsl_integration.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

namespace Utils {

// Numerical inversion of x*f(x) in log(x). Quadrature errors are estimates,
// not certified bounds for an arbitrary black-box Function. No probability
// mass is deliberately truncated and no interpolation is used as the answer.
class NumericalLogCDF {
    struct Panel {
        double lo, hi, mass, error;
        bool operator<(const Panel& other) const { return error < other.error; }
    };
    struct Integral {
        double mass=0, error=0;
        std::vector<Panel> panels;
    };
    const Function& function;
    double xmin, xmax, span, tolerance;
    long double logmin, logmax;
    Integral normalization;
    std::vector<double> breakpoints;
    size_t evaluations=0, iterations=0;
    bool reverse=false;

    [[noreturn]] static void fail(const std::string& message) {
        throw new Exception("NumericalLogCDF: "+message, false);
    }
    double density(double t) {
        ++evaluations;
        // The PDF interface itself is binary64. Endpoint-relative double exp
        // avoids an expensive long-double transcendental on every quadrature
        // node while retaining small distances from either endpoint. Use an
        // absolute log for extremely broad intervals to avoid exp(t) overflow.
        const double coordinate=t<500?(reverse?xmax*std::exp(-t):xmin*std::exp(t)):
            std::exp(double(reverse?logmax-t:logmin+t));
        const double x=std::clamp(coordinate,xmin,xmax);
        const double f=function(x);
        if (!std::isfinite(f) || f<0) fail("negative or nonfinite density");
        const double value=x*f;
        if (!std::isfinite(value)) fail("log-density overflow");
        return value;
    }
    static double callback(double t, void* context) {
        return static_cast<NumericalLogCDF*>(context)->density(t);
    }
    Panel panel(double lo, double hi) {
        gsl_function f{callback,this};
        double result,error,absolute,asc;
        gsl_integration_qk15(&f,lo,hi,&result,&error,&absolute,&asc);
        if (!std::isfinite(result) || !std::isfinite(error) || result<0 || error<0)
            fail("invalid quadrature result");
        return {lo,hi,result,error};
    }
    Integral integrate(double end, bool retain=false, double absoluteTolerance=0) {
        Integral out;
        if (end==0) return out;
        // Avoid relying on a single rule over a broad log-energy interval.
        // Adaptive refinement below controls the integral, not a PDF maximum.
        std::priority_queue<Panel> panels;
        long double mass=0,error=0;
        auto add=[&](Panel p) {panels.push(p);mass+=p.mass;error+=p.error;};
        if (normalization.panels.empty()) {
            const size_t initial=std::max<size_t>(8,size_t(std::ceil(end/2)));
            if (initial>4096) fail("log interval is too wide");
            std::vector<double> cuts=breakpoints;
            for (size_t i=0;i<=initial;++i) cuts.push_back(end*i/initial);
            std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
            if (cuts.size()>4097) fail("too many integration breakpoints");
            for (size_t i=1;i<cuts.size();++i) add(panel(cuts[i-1],cuts[i]));
        } else {
            // Reuse the normalization partition, including zero-density cells.
            // A small nonzero endpoint slice must be integrated on its own;
            // a broad rule spanning a long zero plateau could miss that slice.
            auto include=[&](Panel p) {
                if (p.lo>=end) return;
                add(p.hi<=end?p:panel(p.lo,end));
            };
            if (reverse) {
                for (auto i=normalization.panels.rbegin();i!=normalization.panels.rend();++i)
                    include({span-i->hi,span-i->lo,i->mass,i->error});
            } else for (const auto& p:normalization.panels) include(p);
        }
        const double eps=tolerance/32;
        while (error>std::max(static_cast<long double>(absoluteTolerance),eps*mass)) {
            if (panels.size()>=4096) fail("quadrature did not converge");
            const Panel p=panels.top();panels.pop();
            const double mid=p.lo+(p.hi-p.lo)/2;
            if (!(p.lo<mid && mid<p.hi)) fail("quadrature reached floating-point resolution");
            const Panel left=panel(p.lo,mid),right=panel(mid,p.hi);
            mass+=static_cast<long double>(left.mass)+right.mass-p.mass;
            error+=static_cast<long double>(left.error)+right.error-p.error;
            panels.push(left);panels.push(right);
        }
        out.mass=double(mass);out.error=double(std::max(0.L,error));
        if (!std::isfinite(out.mass) || !std::isfinite(out.error) || out.mass<0)
            fail("integral overflow or cancellation");
        if (retain) {
            while (!panels.empty()) {out.panels.push_back(panels.top());panels.pop();}
            std::sort(out.panels.begin(),out.panels.end(),[](const Panel& a,const Panel& b){return a.lo<b.lo;});
        }
        return out;
    }
    double initialGuess(double probability) const {
        long double remaining=static_cast<long double>(probability)*normalization.mass;
        if (reverse) {
            for (auto i=normalization.panels.rbegin();i!=normalization.panels.rend();++i) {
                if (i->mass>0 && remaining<=i->mass)
                    return span-i->hi+(i->hi-i->lo)*double(remaining/i->mass);
                remaining-=i->mass;
            }
        } else {
            for (const auto& p:normalization.panels) {
                if (p.mass>0 && remaining<=p.mass)
                    return p.lo+(p.hi-p.lo)*double(remaining/p.mass);
                remaining-=p.mass;
            }
        }
        return span/2;
    }
public:
    NumericalLogCDF(const Function& f,double lo,double hi,double relativeTolerance):
        function(f),xmin(std::max(lo,f.Xmin())),xmax(std::min(hi,f.Xmax())) {
        if (std::isnan(f.Xmin()) || std::isnan(f.Xmax()) ||
            !std::isfinite(xmin) || !std::isfinite(xmax) || xmin<=0 || xmax<=xmin ||
            !std::isfinite(relativeTolerance) || relativeTolerance<=0 || relativeTolerance>.1)
            fail("invalid support or tolerance");
        // Separate probability accuracy from a possibly loose rate tolerance.
        // Quadrature receives 1/32 of this small-tail relative error budget.
        tolerance=std::min(relativeTolerance,1e-5);
        if (tolerance<4096*std::numeric_limits<double>::epsilon())
            fail("requested accuracy is below quadrature resolution");
        logmin=std::log(static_cast<long double>(xmin));
        logmax=std::log(static_cast<long double>(xmax));
        const long double ratio=static_cast<long double>(xmax)/xmin;
        span=double(std::isfinite(ratio)?std::log(ratio):logmax-logmin);
        if (!(span>0) || !std::isfinite(span)) fail("unresolved log support");
        for (double x:f.IntegrationBreakpoints()) {
            if (!std::isfinite(x)) fail("nonfinite integration breakpoint");
            if (x>xmin && x<xmax) {
                const long double r=static_cast<long double>(x)/xmin;
                breakpoints.push_back(double(std::isfinite(r)?std::log(r):std::log(static_cast<long double>(x))-logmin));
            }
        }
        normalization=integrate(span,true);
    }
    double IntegralValue() const {return normalization.mass;}
    double IntegralError() const {return normalization.error;}
    size_t Evaluations() const {return evaluations;}
    size_t Iterations() const {return iterations;}
    double Quantile(double u) {
        if (!std::isfinite(u) || u<0 || u>1 || !(normalization.mass>0))
            fail("invalid probability or empty distribution");
        if (u==0) return xmin;
        if (u==1) return xmax;
        reverse=u>.5;
        const double probability=reverse?1-u:u;
        const double target=probability*normalization.mass;
        if (!(target>0)) fail("tail probability underflow");
        const double budget=tolerance*target;
        double lo=0,hi=span,t=initialGuess(probability);
        for (size_t step=0;step<128;++step) {
            ++iterations;
            if (!(lo<t && t<hi)) t=lo+(hi-lo)/2;
            if (!(lo<t && t<hi)) fail("inverse CDF reached floating-point resolution");
            const Integral value=integrate(t,false,budget/32);
            const double residual=value.mass-target;
            const double uncertainty=value.error+probability*normalization.error;
            if (std::abs(residual)+uncertainty<=budget) {
                const double x=double(std::exp(reverse?logmax-t:logmin+t));
                if (!(xmin<=x && x<=xmax) || !std::isfinite(x)) fail("quantile outside support");
                return x;
            }
            if (residual<0) lo=t; else hi=t;
            const double slope=density(t);
            const double next=slope>0?t-residual/slope:std::numeric_limits<double>::quiet_NaN();
            // Safeguarded Newton; retain the complete physical bracket.
            t=std::isfinite(next) && next>lo && next<hi && std::abs(next-t)<.75*(hi-lo)
                ?next:lo+(hi-lo)/2;
        }
        fail("inverse CDF did not converge");
    }
};
}
#endif
