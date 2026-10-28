#include "weighted_common.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace weighted {
Guide::Guide(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read immutable importance guide");
    for (std::string line;std::getline(input,line);) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream row(line);char process;unsigned a,z;std::string code,extra;double score;
        if (!(row>>process>>a>>z>>code>>score) || row>>extra || (process!='P' && process!='D') ||
            !std::isfinite(score) || score<0. || a<2 || z>a ||
            !entries.emplace(std::make_tuple(process,a,z,code),score).second)
            throw std::runtime_error("Invalid or duplicate importance guide row");
    }
    if (entries.empty()) throw std::runtime_error("Empty importance guide");
}
double Guide::score(char process,unsigned a,unsigned z,const std::string& code) const {
    const auto found=entries.find({process,a,z,code});
    if (found==entries.end()) throw std::runtime_error("Missing physical channel in importance guide");
    return found->second;
}
namespace {
Selection choose(const std::vector<double>& p,const std::vector<double>& g,double fraction,double u,bool density) {
    if (p.empty() || p.size()!=g.size() || !std::isfinite(fraction) || fraction<0. || fraction>=1. ||
        !std::isfinite(u) || u<0. || u>=1.) throw std::runtime_error("Invalid weighted channel selection");
    double total=0.,guided=0.;
    for (size_t i=0;i<p.size();++i) {
        if (!std::isfinite(p[i]) || p[i]<0. || !std::isfinite(g[i]) || g[i]<0.)
            throw std::runtime_error("Invalid branch probability or guide score");
        total+=p[i];guided+=p[i]>0. ? (density?g[i]:p[i]*g[i]) : 0.;
    }
    if (std::abs(total-1.)>1e-12 || !std::isfinite(guided)) throw std::runtime_error("Unnormalized branch law");
    std::vector<double> q(p.size());double cumulative=0.;size_t chosen=p.size();
    for (size_t i=0;i<p.size();++i) {
        const double component=p[i]>0. ? (density?g[i]:p[i]*g[i]) : 0.;
        q[i]=guided>0. ? (1.-fraction)*p[i]+fraction*component/guided : p[i];
        if (p[i]>0. && !(q[i]>0.)) throw std::runtime_error("Importance support underflow");
        cumulative+=q[i];
        if (chosen==p.size() && u<cumulative) chosen=i;
    }
    if (chosen==p.size()) for (size_t i=p.size();i--;) if (q[i]>0.) {chosen=i;break;}
    if (chosen==p.size() || p[chosen]<=0.) throw std::runtime_error("Empty importance support");
    return {chosen,p[chosen]/q[chosen]};
}
}
Selection select(const std::vector<double>& p,const std::vector<double>& g,double fraction,double u) {
    return choose(p,g,fraction,u,false);
}
Selection selectDensity(const std::vector<double>& p,const std::vector<double>& g,double fraction,double u) {
    return choose(p,g,fraction,u,true);
}
double reweight(double w,double ratio,unsigned copies) {
    if (!copies) throw std::runtime_error("Zero splitting replication count");
    const double result=w*ratio/copies;
    if (!copies || !std::isfinite(w) || !(w>0.) || !std::isfinite(ratio) || !(ratio>0.) ||
        !std::isfinite(result) || !(result>0.)) throw std::runtime_error("Weight overflow or underflow; no clipping allowed");
    return result;
}
std::string channelCode(unsigned value,unsigned digits) {
    std::ostringstream out;out<<std::setw(digits)<<std::setfill('0')<<value;return out.str();
}
std::vector<double> officialProbabilities(const std::vector<double>& raw) {
    if (raw.empty()) throw std::runtime_error("Missing official PD branches");
    std::vector<double> p(raw.size());double remaining=1.;
    for (size_t i=0;i<raw.size();++i) {
        if (!std::isfinite(raw[i]) || raw[i]<0.) throw std::runtime_error("Invalid official branching fraction");
        p[i]=i+1==raw.size() ? remaining : std::min(remaining,raw[i]);remaining-=p[i];
    }
    return p;
}
}
