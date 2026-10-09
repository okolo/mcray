#include "GZK.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <unistd.h>

namespace {
using namespace mcray;
using namespace Interactions;
class GzkPhoton : public BackgroundIntegral {
public:
    double GetRateS(const Utils::Function& sigma,const Particle& p,double) override {
        return sigma.f(p.Mass()*(p.Mass()+units.GeV))*p.Energy/p.Mass();
    }
    double GetRateK(const Utils::Function& sigma,const Particle&,double) override { return sigma.f(0.5*units.GeV); }
    double GetRateAndSampleS(const Utils::Function& sigma,const Particle& p,Randomizer&,double& s,double error) override {
        s=p.Mass()*(p.Mass()+units.GeV); return GetRateS(sigma,p,error);
    }
    double GetRateAndSampleK(const Utils::Function& sigma,const Particle& p,Randomizer&,double& k,double error) override {
        k=0.5*units.GeV; return GetRateK(sigma,p,error);
    }
    BackgroundIntegral* Clone() const override { return new GzkPhoton; }
};
class NuclearGzkTest : public ::testing::Test {
protected:
    std::filesystem::path dir;
    void Write(const std::string& file,const std::string& text) { std::ofstream(dir/file)<<text; }
    void SetUp() override {
        dir=std::filesystem::temp_directory_path()/("mcray-gzk-policy-"+std::to_string(getpid()));
        std::filesystem::remove_all(dir); std::filesystem::create_directory(dir);
        Write("provenance.txt","mcray-talys-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\ngrid_policy=crpropa-padding\ngzk_compatibility=true\n");
        Write("manifest.json","{}\n");
        Write("masses.tsv","4 2 3727.4 NIST\n3 2 2808.4 NIST\n3 1 2793.9 CRPropaApprox\n2 1 1875.6 NIST\n");
        Write("decay_provenance.txt","mcray-decay-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\n");
        Write("decay_manifest.json","{}\n"); Write("decays.tsv","");
        Write("decay_stable.tsv","4 2\n3 2\n2 1\n");
        Write("fallback_decays.tsv","3 1 00001 1e-99\n");
        Write("unsupported_decays.tsv","");
    }
    void TearDown() override { std::filesystem::remove_all(dir); }
    std::shared_ptr<const NuclearDecayData> Load(bool compatibility) {
        return std::make_shared<NuclearDecayData>(dir.string(),
            std::make_shared<TalysData>(dir.string(),false,compatibility));
    }
    void Error(const std::function<void()>& action,const std::string& message) {
        try { action(); FAIL()<<"Expected "<<message; }
        catch (Utils::Exception* e) { const auto text=e->Message(); delete e; EXPECT_NE(text.find(message),std::string::npos)<<text; }
    }
};
TEST(NuclearGzk, RateModifierBoundaries) {
    EXPECT_DOUBLE_EQ(GZK::NuclearModifier(1,1),1.);
    EXPECT_DOUBLE_EQ(GZK::NuclearModifier(56,0),0.);
    for (unsigned int a : {2u,8u})
        EXPECT_NEAR(GZK::NuclearModifier(a,2),0.85*std::pow(2.,2./3.),1e-14);
    EXPECT_DOUBLE_EQ(GZK::NuclearModifier(9,2),1.7);
    EXPECT_DOUBLE_EQ(GZK::NuclearModifier(56,26),0.85*26);
}
TEST_F(NuclearGzkTest, SameBundleSelectsStrictVersusApproximateMassesAndDripping) {
    auto strict=Load(false),compatible=Load(true);
    EXPECT_EQ(strict->Nuclei().Find(3,1),nullptr);
    EXPECT_TRUE(compatible->Nuclei().Get(3,1).approximateMass);
    Error([&]{strict->ValidateNetwork(4,2,false,true);},"Missing TALYS isotope");
    const auto reached=compatible->ValidateNetwork(4,2,false,true);
    EXPECT_TRUE(reached.count({1,0})); EXPECT_TRUE(reached.count({1,1}));
    EXPECT_TRUE(reached.count({3,1}));
    EXPECT_FALSE(strict->Nuclei().CRPropaFallback()); EXPECT_TRUE(compatible->Nuclei().CRPropaFallback());
    EXPECT_TRUE(std::isfinite(compatible->Get(3,1)[0].restRate));
    EXPECT_NEAR(compatible->Get(3,1)[0].restRate*units.second*1e-99,1.,1e-14);
    Particle p(Nucleus,0.); compatible->Nuclei().SetParticle(p,3,1); p.Energy=3e20*units.eV; p.Weight=2.5;
    NuclearDecay decay(compatible); Randomizer random(20261010); std::vector<Particle> out;
    ASSERT_TRUE(decay.GetSecondaries(p,out,random)); ASSERT_EQ(out.size(),2u);
    EXPECT_EQ(out[0].Type,Neutron); EXPECT_EQ(out[1].MassNumber(),2u);
    EXPECT_NEAR((out[0].Energy+out[1].Energy)/p.Energy,1.,1e-14);
    for (const auto& product:out) EXPECT_DOUBLE_EQ(product.Weight,p.Weight);
}
TEST_F(NuclearGzkTest, StrictAcceptsVerifiedLightClosureAndBreaksPureSpectators) {
    auto data=Load(false); EXPECT_FALSE(data->ValidateNetwork(3,2,false,true).empty());
    Particle p(Nucleus,0.); data->Nuclei().SetParticle(p,3,2); p.Energy=3e20*units.eV;p.Weight=3.;
    for (unsigned int z : {0u,2u}) {
        std::vector<Particle> out; data->AppendResidual(p,2,z,2e20*units.eV,out);
        ASSERT_EQ(out.size(),2u);
        for (const auto& q:out) {
            EXPECT_EQ(q.Type,z ? Proton : Neutron); EXPECT_DOUBLE_EQ(q.Energy,1e20*units.eV);
            EXPECT_DOUBLE_EQ(q.Weight,3.); EXPECT_EQ(q.SourceParticle,p.SourceParticle);
        }
    }
}
TEST_F(NuclearGzkTest, EqualGammaRatesTargetSelectionAndClone) {
    auto data=Load(true); GZK process(new GzkPhoton,20261010,data),free(new GzkPhoton);
    Particle p(Nucleus,0.); data->Nuclei().SetParticle(p,4,2);p.Energy=p.Mass()*1e11;
    const auto rates=process.NuclearTargetRates(p);
    for (size_t i=0;i<2;++i) {
        Particle nucleon(i ? Neutron : Proton,0.);nucleon.Energy=nucleon.Mass()*1e11;
        EXPECT_NEAR(rates[i]/free.Rate(nucleon),0.85*std::pow(2.,2./3.),1e-13);
        EXPECT_DOUBLE_EQ(process.Rate(nucleon),free.Rate(nucleon));
    }
    std::unique_ptr<RandomInteraction> clone(process.Clone());
    EXPECT_DOUBLE_EQ(clone->Rate(p),process.Rate(p));
    p.SetPsbNucleus(4); Error([&]{process.Rate(p);},"selected nuclear-data mass");
}
TEST_F(NuclearGzkTest, SophiaProductsConserveBaryonsChargeEnergyAndWeights) {
    auto data=Load(true); GZK process(new GzkPhoton,20261010,data);
    Particle p(Nucleus,0.);data->Nuclei().SetParticle(p,4,2);p.Energy=4e20*units.eV;p.Weight=2.5;
    Randomizer random(20261010);unsigned onProton=0;
    for (int i=0;i<200;++i) {
        std::vector<Particle> out;ASSERT_TRUE(process.GetSecondaries(p,out,random));
        unsigned baryons=0;double charge=0.,energy=0.;bool spectator=false;
        for (const auto& q:out) {
            baryons+=q.MassNumber();charge+=q.ElectricCharge();energy+=q.Energy;
            EXPECT_DOUBLE_EQ(q.Weight,p.Weight);EXPECT_EQ(q.SourceParticle,p.SourceParticle);
            if (q.Type==Nucleus) {
                spectator=true;EXPECT_EQ(q.MassNumber(),3u);
                EXPECT_NEAR(q.Energy/p.Energy,0.75,1e-14);
                onProton+=q.AtomicNumber()==1;
            }
        }
        EXPECT_TRUE(spectator);EXPECT_EQ(baryons,4u);EXPECT_DOUBLE_EQ(charge,2.);
        EXPECT_NEAR(energy/p.Energy,1.,1e-6);
        EXPECT_DOUBLE_EQ(p.Energy,4e20*units.eV);
    }
    const auto rates=process.NuclearTargetRates(p);const double probability=rates[0]/(rates[0]+rates[1]);
    EXPECT_NEAR(onProton,200*probability,5*std::sqrt(200*probability*(1-probability)));
}
TEST_F(NuclearGzkTest, UnsupportedChannelsAndMissingPdRemainErrorsInBothPolicies) {
    auto data=Load(true);Error([&]{data->ValidateNetwork(4,2,true,true);},"photodisintegration data");
    Write("unsupported_decays.tsv","3 2 UNSUPPORTED\n");
    for (bool mode : {false,true}) { auto selected=Load(mode); Error([&]{selected->ValidateNetwork(3,2,false,true);},"Unsupported source decay mode"); }
}
}
