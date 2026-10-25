#pragma once
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace weighted {
class Guide {
public:
    explicit Guide(const std::string& path);
    double score(char process,unsigned a,unsigned z,const std::string& channel) const;
private:
    std::map<std::tuple<char,unsigned,unsigned,std::string>,double> entries;
};
struct Selection { size_t index; double likelihood; };
Selection select(const std::vector<double>& physical,const std::vector<double>& guide,
                 double fraction,double uniform);
// Defensive mixture with a topology proposal over physically possible channels.
Selection selectDensity(const std::vector<double>& physical,const std::vector<double>& density,
                        double fraction,double uniform);
double reweight(double weight,double likelihood,unsigned copies=1);
std::string channelCode(unsigned value,unsigned digits);
// Effective law of the pinned CRPropa nearest-node cumulative branch sampler.
std::vector<double> officialProbabilities(const std::vector<double>& branches);
}
