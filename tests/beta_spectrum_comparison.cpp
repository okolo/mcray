#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
struct Moments { double sum=0.,sum2=0.;unsigned n=0;void Add(double e){sum+=e;sum2+=e*e;++n;}
    double Mean()const{return sum/n;} double SE()const{return std::sqrt(std::max(0.,(sum2-sum*sum/n)/(n-1)/n));}};
std::array<Moments,2> mc,ref;
std::array<double,2> continuum{},discrete{};
unsigned expected=0;
void Read(const char* file,std::array<Moments,2>& moments,double gamma) {
    std::ifstream input(file);if(!input) throw std::runtime_error("Cannot open beta spectrum file");
    std::map<std::pair<unsigned,unsigned>,std::array<unsigned,2>> counts;
    for(std::string line;std::getline(input,line);) {
        if(line.empty() || line[0]=='#') continue;
        unsigned batch,event;std::string species;double weight,energy;
        std::istringstream row(line);if(!(row>>batch>>event>>species>>weight>>energy) || weight!=1. || !std::isfinite(energy) || energy<0.)
            throw std::runtime_error("Invalid beta product");
        auto& count=counts[{batch,event}];
        if(species=="electron" || species=="positron") { moments[0].Add(energy/1e6/gamma);++count[0]; }
        if(species=="neutrinoE" || species=="neutrinoAE") { moments[1].Add(energy/1e6/gamma);++count[1]; }
    }
    if(counts.size()!=expected) throw std::runtime_error("Incomplete conditional ensemble");
    for(const auto& entry:counts) if(entry.second!=std::array<unsigned,2>{{1,1}})
        throw std::runtime_error("Conditional event needs exactly one charged lepton and one neutrino");
}
TEST(BetaSpectrumReference, ConditionalEventsAreComplete) {for(int i=0;i<2;++i){EXPECT_EQ(mc[i].n,expected);EXPECT_EQ(ref[i].n,expected);}}
TEST(BetaSpectrumReference, ProgramMeansAgreeWithinSamplingUncertainty) {
    for(int i=0;i<2;++i) EXPECT_NEAR(mc[i].Mean(),ref[i].Mean(),5.*std::hypot(mc[i].SE(),ref[i].SE()));
}
TEST(BetaSpectrumReference, BothSamplersMatchTheApprovedDiscreteSpectrum) {
    for(int i=0;i<2;++i){EXPECT_NEAR(mc[i].Mean(),discrete[i],5.*mc[i].SE());EXPECT_NEAR(ref[i].Mean(),discrete[i],5.*ref[i].SE());}
}
int main(int argc,char** argv) {
    if(argc!=6) return 2;
    try {
        const double q=std::stod(argv[3]),gamma=std::stod(argv[4]),me=0.511;
        expected=std::stoul(argv[5]);if(q<=0. || gamma<1. || expected<2) return 2;
        Read(argv[1],mc,gamma);Read(argv[2],ref,gamma);
        double norm=0.;std::array<double,2> integrals{};
        const int intervals=20000;
        for(int i=0;i<=intervals;++i) {
            const double e=me+q*i/intervals,nu=q+me-e;
            const double weight=(i==0 || i==intervals)?1.:(i%2?4.:2.);
            const double density=weight*e*std::sqrt(std::max(0.,e*e-me*me))*nu*nu;
            norm+=density;integrals[0]+=density*e;integrals[1]+=density*nu;
        }
        for(int i=0;i<2;++i) continuum[i]=integrals[i]/norm;
        norm=0.;integrals={};
        for(int i=1;i<=50;++i) {
            const double e=me+q*i/50.,nu=q+me-e;
            const double density=e*std::sqrt(std::max(0.,e*e-me*me))*nu*nu;
            const double intervalMean=e-q/100.;
            norm+=density;integrals[0]+=density*intervalMean;integrals[1]+=density*(q+me-intervalMean);
        }
        for(int i=0;i<2;++i) {
            discrete[i]=integrals[i]/norm;
            std::cout<<"# species="<<(i?"neutrino":"charged-lepton")<<" continuum_mean_MeV="<<continuum[i]
              <<" discrete_mean_MeV="<<discrete[i]<<" relative_discretization_bias="<<discrete[i]/continuum[i]-1
              <<" mcray_mean_MeV="<<mc[i].Mean()<<" CRPropa_mean_MeV="<<ref[i].Mean()<<'\n';
        }
    } catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 2;}
    int testArgc=1;::testing::InitGoogleTest(&testArgc,argv);return RUN_ALL_TESTS();
}
