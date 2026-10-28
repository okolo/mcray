#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>

namespace {
struct History { double a=0.,z=0.; bool survived=false; };
std::map<std::pair<unsigned,unsigned>,History> histories;
unsigned expected=0,primaryA=0,primaryZ=0;
double rate=-1.,distance=-1.;
bool alphaOnly=false;
std::string primaryName;
unsigned alphaCount=0;
TEST(NuclearDecayOutput, AllHistoriesReachOutput) { EXPECT_EQ(histories.size(),expected); }
TEST(NuclearDecayOutput, BaryonNumberAndChargeIncludeLeptons) {
    for (const auto& entry:histories) {
        EXPECT_NEAR(entry.second.a,primaryA,1e-10);
        EXPECT_NEAR(entry.second.z,primaryZ,1e-10);
    }
}
TEST(NuclearDecayOutput, PrimarySurvivalHasExpectedExponentialProbability) {
    ASSERT_GE(rate,0.); ASSERT_GE(distance,0.);
    const double probability=std::exp(-rate*distance);
    unsigned survived=0; for (const auto& entry:histories) survived+=entry.second.survived;
    // Diagnostic five-sigma binomial check, not rare-output acceptance precision.
    const double sigma=std::sqrt(expected*probability*(1.-probability));
    EXPECT_LE(std::abs(survived-expected*probability),5.*sigma+1.);
}
TEST(NuclearDecayOutput, Be8ProducesTwoAlphasWhenRequested) {
    if (!alphaOnly) GTEST_SKIP()<<"Not the isolated Be8 control";
    EXPECT_EQ(alphaCount,2*expected);
}
}
int main(int argc,char** argv) {
    std::string file;
    for (int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if (arg=="--alpha-only") { alphaOnly=true; continue; }
        if (i+1>=argc) return 2;
        const std::string value=argv[++i];
        if (arg=="--mc") file=value;
        else if (arg=="--primaries") expected=std::stoul(value);
        else if (arg=="--A") primaryA=std::stoul(value);
        else if (arg=="--Z") primaryZ=std::stoul(value);
        else return 2;
    }
    if (file.empty() || !expected || primaryA<2 || !primaryZ) return 2;
    std::ostringstream name; name<<'A'; if(primaryA<10) name<<'0'; name<<primaryA<<'Z'; if(primaryZ<10) name<<'0'; name<<primaryZ;
    primaryName=name.str();
    std::ifstream input(file); if(!input) return 2;
    const std::regex identity("A([0-9]+)Z([0-9]+)");
    for (std::string line;std::getline(input,line);) {
        if (line.empty()) continue;
        if (line[0]=='#') {
            for (const auto& key: {std::string("event_rate_per_Mpc="),std::string("L_mpc=")}) {
                const auto pos=line.find(key); if(pos!=std::string::npos)
                    (key[0]=='L'?distance:rate)=std::stod(line.substr(pos+key.size()));
            }
            continue;
        }
        unsigned batch,event; std::string species; double weight,energy;
        std::istringstream row(line);
        if(!(row>>batch>>event>>species>>weight>>energy) || !std::isfinite(weight) ||
           weight!=1. || !std::isfinite(energy) || energy<0.) return 2;
        double a=0.,z=0.; std::smatch match;
        if(std::regex_match(species,match,identity)) { a=std::stoul(match[1]); z=std::stoul(match[2]); }
        else if(species=="proton") { a=1.; z=1.; }
        else if(species=="neutron") a=1.;
        else if(species=="electron") z=-1.;
        else if(species=="positron") z=1.;
        else if(species.rfind("neutrino",0)!=0 && species!="photon") return 2;
        auto& history=histories[{batch,event}]; history.a+=a; history.z+=z;
        history.survived|=species==primaryName;
        if (species=="A04Z02") ++alphaCount;
        if (alphaOnly && species!="A04Z02") return 2;
    }
    int testArgc=1; ::testing::InitGoogleTest(&testArgc,argv); return RUN_ALL_TESTS();
}
