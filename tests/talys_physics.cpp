#include "TalysPhotoDisintegration.h"
#include "PhotoDisintegration.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <unistd.h>

namespace {
using namespace mcray;
using namespace Interactions;
class FixedK : public BackgroundIntegral {
public:
    double GetRateS(const Utils::Function&, const Particle&, double) override { return 0.; }
    double GetRateK(const Utils::Function& sigma, const Particle&, double) override { return sigma.f(10*units.MeV); }
    double GetRateAndSampleS(const Utils::Function&, const Particle&, Randomizer&, double&, double) override { return 0.; }
    double GetRateAndSampleK(const Utils::Function&, const Particle&, Randomizer&, double&, double) override { return 0.; }
    BackgroundIntegral* Clone() const override { return new FixedK; }
};
class TalysTest : public ::testing::Test {
protected:
    std::filesystem::path dir;
    void SetUp() override {
        dir = std::filesystem::temp_directory_path() / ("mcray-talys-fixture-" + std::to_string(getpid()));
        std::filesystem::remove_all(dir); std::filesystem::create_directory(dir);
        Write("provenance.txt", "mcray-talys-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\ngrid_policy=crpropa-padding\n");
        Write("manifest.json", "{}\n");
        std::string grid;
        for (int i=1;i<=513;++i) grid += std::to_string(i) + '\n';
        Write("grid_talys.tsv", grid); Write("grid_external.tsv",grid);
        Write("masses.tsv", "56 26 52192 NIST\n50 25 46600 NIST\n47 24 43804 NIST\n42 21 39144 NIST\n4 2 3728 NIST\n3 1 2796 NIST\n2 1 1864 NIST\n");
        Write("terminals.tsv", "50 25\n47 24\n42 21\n");
        Write("total.tsv", Row("56 26",50) + Row("4 2",1) + Row("3 1",1) + Row("2 1",1));
        Write("channels.tsv", Row("56 26 510000",2) + Row("56 26 720000",6) + Row("56 26 300102",4)
              + Row("4 2 011000",1) + Row("3 1 100000",1) + Row("2 1 100000",1));
    }
    void TearDown() override { std::filesystem::remove_all(dir); }
    void Write(const std::string& file, const std::string& contents) { std::ofstream(dir/file) << contents; }
    std::string Row(const std::string& prefix, double sigma) {
        std::string row = prefix;
        for (int i=0;i<513;++i) row += " " + std::to_string(i==512 ? 0. : sigma);
        return row + '\n';
    }
    void ExpectError(const std::function<void()>& action, const std::string& text) {
        try { action(); FAIL() << "Expected error containing " << text; }
        catch (Utils::Exception* error) {
            std::string message=error->Message(); delete error;
            EXPECT_NE(message.find(text),std::string::npos) << message;
        }
    }
};
TEST_F(TalysTest, NetworkMassesRatesAndInterpolationUseTheAgreedContracts) {
    auto data = std::make_shared<TalysData>(dir.string());
    EXPECT_EQ(data->ValidateNetwork(56,26).size(),9u);
    Particle p(Nucleus,0.); data->SetParticle(p,56,26); p.Energy=1e20*units.eV;
    EXPECT_DOUBLE_EQ(p.Mass(),52192*units.MeV);
    TalysPhotoDisintegration process(new FixedK,data);
    const auto rates=process.ChannelRates(p);
    EXPECT_NEAR(rates.total/(1e-3*units.barn),50.,1e-12);
    EXPECT_NEAR(rates.exclusiveSum/rates.total,12./50.,1e-14);
    std::unique_ptr<RandomInteraction> clone(process.Clone());
    EXPECT_DOUBLE_EQ(process.Rate(p),clone->Rate(p));
    const auto& sigma=*data->Get(56,26).total;
    EXPECT_DOUBLE_EQ(sigma.f(0.9*units.MeV),0.); EXPECT_DOUBLE_EQ(sigma.f(514*units.MeV),0.);
    EXPECT_NEAR(sigma.f(std::sqrt(512.*513.)*units.MeV)/sigma.f(512*units.MeV),0.5,1e-10);
    data->SetParticle(p,50,25); EXPECT_DOUBLE_EQ(process.Rate(p),0.);
    data->SetParticle(p,1,1); EXPECT_EQ(p.Type,Proton); EXPECT_DOUBLE_EQ(p.Mass(),Particle::Mass(Proton));
    EXPECT_DOUBLE_EQ(process.Rate(p),0.);
    p.SetPsbNucleus(56); ExpectError([&] { process.Rate(p); },"tabulated nuclear mass");
    PhotoDisintegration psb(new FixedK);
    data->SetParticle(p,56,26); ExpectError([&] { psb.Rate(p); },"PSB requires");
}
TEST_F(TalysTest, OneChannelPerEventConservesBaryonsChargeEnergyAndWeight) {
    auto data=std::make_shared<TalysData>(dir.string());
    TalysPhotoDisintegration process(new FixedK,data);
    Particle p(Nucleus,0.); data->SetParticle(p,56,26); p.Energy=1e20*units.eV; p.Weight=2.5;
    Randomizer random(20261008);
    std::map<unsigned,int> counts;
    for (int event=0;event<12000;++event) {
        std::vector<Particle> products;
        ASSERT_TRUE(process.GetSecondaries(p,products,random));
        ASSERT_FALSE(products.empty()); ++counts[products.front().MassNumber()];
        unsigned a=0,z=0; double energy=0.;
        for (const auto& q:products) {
            q.ValidateIdentity(); a+=q.MassNumber(); z+=q.AtomicNumber(); energy+=q.Energy;
            EXPECT_DOUBLE_EQ(q.Weight,p.Weight); EXPECT_EQ(q.SourceParticle,p.SourceParticle);
            EXPECT_DOUBLE_EQ(q.Energy,p.Energy*q.MassNumber()/56.);
            if (q.Type==Nucleus) EXPECT_EQ(q.MassModel(),Particle::NuclearMassModel::TALYS);
        }
        EXPECT_EQ(a,56u); EXPECT_EQ(z,26u); EXPECT_NEAR(energy/p.Energy,1.,1e-14);
    }
    for (const auto& expected : std::map<unsigned,double>{{50,2./12.},{47,6./12.},{42,4./12.}}) {
        const double mean=12000*expected.second;
        EXPECT_NEAR(counts[expected.first],mean,6*std::sqrt(mean*(1-expected.second)));
    }
    EXPECT_TRUE(process.KillsPrimary()); EXPECT_EQ(p.Name(),"A56Z26"); EXPECT_DOUBLE_EQ(p.Weight,2.5);
}
TEST_F(TalysTest, MissingReachableDataAndUnapprovedTerminalsAreErrors) {
    Write("masses.tsv","56 26 52192 NIST\n50 25 46600 NIST\n47 24 43804 NIST\n42 21 39144 NIST\n4 2 3728 NIST\n3 1 2796 NIST\n");
    ExpectError([&] { TalysData data(dir.string()); },"Missing TALYS mass");
    SetUp(); Write("terminals.tsv","56 26\n");
    ExpectError([&] { TalysData data(dir.string()); },"Unapproved");
    SetUp(); Write("channels.tsv",Row("56 26 510000",2));
    TalysData data(dir.string());
    ExpectError([&] { data.ValidateNetwork(4,2); },"Incomplete reachable");
}
TEST_F(TalysTest, InvalidChannelsAndPositiveTotalWithoutExclusiveRateAreErrors) {
    Write("channels.tsv",Row("56 26 999999",1));
    ExpectError([&] { TalysData data(dir.string()); },"Invalid TALYS residual");
    SetUp(); Write("channels.tsv",Row("56 26 510000",0));
    auto data=std::make_shared<TalysData>(dir.string());
    TalysPhotoDisintegration process(new FixedK,data);
    Particle p(Nucleus,0.); data->SetParticle(p,56,26);
    ExpectError([&] { process.ChannelRates(p); },"zero exclusive-channel sum");
    SetUp(); Write("total.tsv",Row("56 26",0));
    auto zeroTotal=std::make_shared<TalysData>(dir.string());
    TalysPhotoDisintegration inconsistent(new FixedK,zeroTotal);
    zeroTotal->SetParticle(p,56,26);
    ExpectError([&] { inconsistent.ChannelRates(p); },"zero total rate");
    SetUp(); Write("total.tsv",Row("56 26",-1));
    ExpectError([&] { TalysData invalid(dir.string()); },"cross-section value");
}
}
