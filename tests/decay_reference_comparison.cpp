#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <tuple>
using Key=std::tuple<unsigned,unsigned,double>;
using Values=std::array<double,3>; // mass/MeV, rate/s, rate/Mpc
std::map<Key,Values> nativeRates,referenceRates;
void Read(const char* file,std::map<Key,Values>& result) {
    std::ifstream input(file); if(!input) throw std::runtime_error("Cannot read rate report");
    for(std::string line;std::getline(input,line);) {
        if(line.empty() || line[0]=='#') continue;
        unsigned a,z; double gamma; Values values;
        std::istringstream row(line);
        if(!(row>>a>>z>>gamma>>values[0]>>values[1]>>values[2]) || !result.emplace(Key{a,z,gamma},values).second)
            throw std::runtime_error("Invalid or duplicate rate report row");
    }
}
TEST(NuclearDecayReference, CompleteIsotopeGammaInventory) {
    ASSERT_EQ(nativeRates.size(),632u); ASSERT_EQ(referenceRates.size(),632u);
    for(const auto& entry:nativeRates) EXPECT_EQ(referenceRates.count(entry.first),1u);
}
TEST(NuclearDecayReference, TabulatedMassesAgree) {
    for(const auto& entry:nativeRates) {
        const auto found=referenceRates.find(entry.first); ASSERT_NE(found,referenceRates.end());
        EXPECT_NEAR(entry.second[0]/found->second[0],1.,1e-11);
    }
}
TEST(NuclearDecayReference, RatesAgreeWithinLifetimeSerializationTolerance) {
    unsigned stable=0;
    for(const auto& entry:nativeRates) {
        const auto found=referenceRates.find(entry.first); ASSERT_NE(found,referenceRates.end());
        for(int i=1;i<=2;++i) {
            EXPECT_TRUE(std::isfinite(entry.second[i])); EXPECT_GE(entry.second[i],0.);
            if(found->second[i]==0.) EXPECT_DOUBLE_EQ(entry.second[i],0.);
            else EXPECT_NEAR(entry.second[i]/found->second[i],1.,1e-5);
        }
        stable+=found->second[1]==0.;
    }
    EXPECT_EQ(stable,264u);
}
int main(int argc,char** argv) {
    if(argc!=3) return 2;
    try { Read(argv[1],nativeRates);Read(argv[2],referenceRates); }
    catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 2; }
    int testArgc=1;::testing::InitGoogleTest(&testArgc,argv);return RUN_ALL_TESTS();
}
