#include "NuclearDecay.h"
#include "PropagationEngine.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <unistd.h>
namespace {
using namespace mcray;
using namespace Interactions;
class LineBackground : public BackgroundIntegral {
public:
    double GetRateK(const Utils::Function& sigma,const Particle&,double) override { return sigma.f(1.5*units.MeV); }
    double GetRateS(const Utils::Function&,const Particle&,double) override { return 0.; }
    double GetRateAndSampleK(const Utils::Function&,const Particle&,Randomizer&,double&,double) override { return 0.; }
    double GetRateAndSampleS(const Utils::Function&,const Particle&,Randomizer&,double&,double) override { return 0.; }
    BackgroundIntegral* Clone() const override { return new LineBackground; }
};
class NuclearPhotonsTest : public ::testing::Test {
protected:
    std::filesystem::path dir;
    void Write(const std::string& name,const std::string& value) { std::ofstream(dir/name)<<value; }
    std::string Sigma(double value) { std::string s;for(int i=0;i<513;++i)s+=' '+std::to_string(value);return s+'\n'; }
    void SetUp() override {
        dir=std::filesystem::temp_directory_path()/("mcray-photons-"+std::to_string(getpid()));
        std::filesystem::remove_all(dir);std::filesystem::create_directory(dir);
        Write("provenance.txt","mcray-talys-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\ngrid_policy=crpropa-padding\n");
        Write("manifest.json","{}\n");Write("masses.tsv","4 2 3727 NIST\n3 2 2808 NIST\n");
        std::string grid;for(int i=1;i<=513;++i)grid+=std::to_string(i)+'\n';
        Write("grid_talys.tsv",grid);Write("grid_external.tsv",grid);Write("terminals.tsv","");
        Write("total.tsv","4 2"+Sigma(2.));Write("channels.tsv","4 2 100000"+Sigma(2.));
        Write("photon_provenance.txt","mcray-nuclear-photons-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\n");
        Write("photon_manifest.json","{}\n");Write("pd_photon_exceptions.tsv","");Write("pd_photon_transitions.tsv","4 2 3 2"+Sigma(2.));
        Write("pd_photons.tsv","4 2 3 2 1"+Sigma(1.));
        Write("decay_provenance.txt","mcray-decay-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\n");
        Write("decay_manifest.json","{}\n");Write("decay_stable.tsv","3 2\n");Write("decays.tsv","4 2 00001 1\n");
        Write("decay_photon_channels.tsv","4 2 00001\n");Write("decay_photons.tsv","4 2 00001 1000 1\n4 2 00001 2000 0\n");
    }
    void TearDown() override { std::filesystem::remove_all(dir); }
    std::shared_ptr<const TalysData> Load(bool photons=true) { return std::make_shared<TalysData>(dir.string(),true,false,photons); }
    Particle Primary(const TalysData& data) { Particle p(Nucleus,0.);data.SetParticle(p,4,2);p.Energy=p.Mass()*1e8;p.Weight=3.;return p; }
    void Error(const std::function<void()>& action,const std::string& part) {
        try { action();FAIL()<<"Expected "<<part; }catch(Utils::Exception* e) {const auto s=e->Message();delete e;EXPECT_NE(s.find(part),std::string::npos)<<s;}
    }
};
TEST_F(NuclearPhotonsTest, NativeConvolutionConditionalLinesAndDefaultSaturation) {
    auto data=Load();auto p=Primary(*data);TalysPhotoDisintegration process(new LineBackground,data,true);
    auto yields=process.PhotonYields(p,{3,2});ASSERT_EQ(yields.size(),1u);EXPECT_DOUBLE_EQ(yields[0].probability,0.5);
    Write("pd_photons.tsv","4 2 3 2 1"+Sigma(4.));data=Load();
    TalysPhotoDisintegration saturated(new LineBackground,data,true),strict(new LineBackground,data,true,true);
    yields=saturated.PhotonYields(p,{3,2});EXPECT_DOUBLE_EQ(yields[0].probability,1.);EXPECT_DOUBLE_EQ(yields[0].rawYield,2.);
    Error([&]{strict.PhotonYields(p,{3,2});},"exceeds one");
    Write("pd_photon_transitions.tsv","4 2 3 2"+Sigma(0.));data=Load();
    TalysPhotoDisintegration zero(new LineBackground,data,true),zeroStrict(new LineBackground,data,true,true);
    EXPECT_TRUE(std::isinf(zero.PhotonYields(p,{3,2})[0].rawYield));EXPECT_DOUBLE_EQ(zero.PhotonYields(p,{3,2})[0].probability,1.);
    Error([&]{zeroStrict.PhotonYields(p,{3,2});},"zero transition rate");
}
TEST_F(NuclearPhotonsTest, PdBoostWeightsAncestryAndDocumentedEnergyExcess) {
    Write("pd_photons.tsv","4 2 3 2 1"+Sigma(2.));auto data=Load();auto p=Primary(*data);
    TalysPhotoDisintegration process(new LineBackground,data,true);std::unique_ptr<RandomInteraction> clone(process.Clone());
    Randomizer random(20261011);double photonMean=0.;
    for(int i=0;i<1000;++i) {
        std::vector<Particle> products;ASSERT_TRUE(clone->GetSecondaries(p,products,random));ASSERT_EQ(products.size(),3u);
        const auto& photon=products.back();EXPECT_EQ(photon.Type,Photon);EXPECT_DOUBLE_EQ(photon.Weight,p.Weight);
        EXPECT_EQ(photon.SourceParticle,p.SourceParticle);EXPECT_GE(photon.Energy,0.);
        const double scale=products[0].Energy/products[0].Mass()*units.MeV;
        EXPECT_LE(photon.Energy,2.*scale);photonMean+=photon.Energy/scale;
        EXPECT_NEAR((products[0].Energy+products[1].Energy)/p.Energy,1.,1e-14);
    }
    EXPECT_NEAR(photonMean/1000.,1.,5./std::sqrt(3.*1000.));
}
TEST_F(NuclearPhotonsTest, DecayEmitsOncePerCompoundChannelAndLeavesNuclearEnergy) {
    auto nuclei=Load(false);auto data=std::make_shared<NuclearDecayData>(dir.string(),nuclei,true);
    NuclearDecay process(data,true);std::unique_ptr<RandomInteraction> clone(process.Clone());auto p=Primary(*nuclei);
    Randomizer random(20261011);double mean=0.;
    for(int i=0;i<1000;++i) {
        std::vector<Particle> out;ASSERT_TRUE(clone->GetSecondaries(p,out,random));ASSERT_EQ(out.size(),3u);
        EXPECT_EQ(out[0].Type,Photon);EXPECT_DOUBLE_EQ(out[0].Weight,p.Weight);EXPECT_EQ(out[1].Type,Neutron);
        EXPECT_NEAR((out[1].Energy+out[2].Energy)/p.Energy,1.,1e-14);
        const double scale=p.Energy/p.Mass()*units.MeV;EXPECT_LE(out[0].Energy,2.*scale);mean+=out[0].Energy/scale;
    }
    EXPECT_NEAR(mean/1000.,1.,5./std::sqrt(3.*1000.));
}

TEST_F(NuclearPhotonsTest, BernoulliLinesHaveIndependentExpectedMultiplicities) {
    Write("pd_photons.tsv","4 2 3 2 1"+Sigma(0.5)+"4 2 3 2 2"+Sigma(1.5));
    auto data=Load();auto p=Primary(*data);TalysPhotoDisintegration process(new LineBackground,data,true);
    Randomizer random(20261011);unsigned photons=0,both=0;
    for(unsigned i=0;i<4000;++i) {
        std::vector<Particle> out;ASSERT_TRUE(process.GetSecondaries(p,out,random));
        photons+=out.size()-2;both+=out.size()==4;
    }
    EXPECT_NEAR(photons,4000.,5.*std::sqrt(4000.*0.375));
    EXPECT_NEAR(both,4000.*0.1875,5.*std::sqrt(4000.*0.1875*0.8125));
}
TEST_F(NuclearPhotonsTest, CompoundEmissionAttachesPhotonLinesOnlyOnce) {
    Write("masses.tsv","12 6 11178 NIST\n8 4 7454 NIST\n4 2 3727 NIST\n3 2 2808 NIST\n");
    Write("decays.tsv","12 6 00200 1\n");Write("decay_stable.tsv","8 4\n4 2\n3 2\n");
    Write("decay_photon_channels.tsv","12 6 00200\n");Write("decay_photons.tsv","12 6 00200 1000 1\n");
    auto nuclei=Load(false);auto data=std::make_shared<NuclearDecayData>(dir.string(),nuclei,true);
    NuclearDecay process(data,true);Particle p(Nucleus,0.);nuclei->SetParticle(p,12,6);p.Energy=p.Mass()*1e8;
    Randomizer random(20261011);std::vector<Particle> out;ASSERT_TRUE(process.GetSecondaries(p,out,random));ASSERT_EQ(out.size(),4u);
    EXPECT_EQ(out[0].Type,Photon);for(size_t i=1;i<4;++i)EXPECT_EQ(out[i].MassNumber(),4u);
}
TEST_F(NuclearPhotonsTest, MissingCoverageAndInvalidLinesAreErrors) {
    Write("pd_photon_transitions.tsv","");Error([&]{Load();},"no transition");
    Write("pd_photons.tsv","");auto data=Load();Error([&]{data->RequirePhotonCoverage(4,2);},"Missing PD photon transition");
    auto nuclei=Load(false);Write("decay_photon_channels.tsv","");
    Error([&]{NuclearDecayData d(dir.string(),nuclei,true);},"Missing decay photon channel coverage");
    Write("decay_photon_channels.tsv","4 2 00001\n");Write("decay_photons.tsv","4 2 00001 1000 1.1\n");
    Error([&]{NuclearDecayData d(dir.string(),nuclei,true);},"Invalid decay photon");
}

TEST_F(NuclearPhotonsTest, ApprovedLightExceptionIsExplicitAndOtherGapsFail) {
    Write("pd_photon_transitions.tsv","");Write("pd_photons.tsv","");
    Write("pd_photon_exceptions.tsv","4 2 3 2\n");auto data=Load();data->RequirePhotonCoverage(4,2);
    auto p=Primary(*data);TalysPhotoDisintegration process(new LineBackground,data,true);
    EXPECT_TRUE(process.PhotonYields(p,{3,2}).empty());
    Write("pd_photon_exceptions.tsv","4 2 2 2\n");Error([&]{Load();},"Unapproved PD photon");
}

TEST_F(NuclearPhotonsTest, StrictDataErrorsCannotBeBypassedByEngineResampling) {
    if (!cosmology.IsInitialized()) cosmology.Init(0.1);
    Write("pd_photons.tsv","4 2 3 2 1"+Sigma(4.));auto data=Load();
    class CountedStrict : public TalysPhotoDisintegration {
    public:
        mutable unsigned attempts=0;
        explicit CountedStrict(std::shared_ptr<const TalysData> data):TalysPhotoDisintegration(new LineBackground,data,true,true){}
        double Rate(const Particle&) const override {return 1e5/units.Mpc;}
        bool GetSecondaries(Particle& p,std::vector<Particle>& out,Randomizer& rng) const override {
            ++attempts;return TalysPhotoDisintegration::GetSecondaries(p,out,rng);
        }
    };
    for(double variation:{0.,0.01}) {
        ParticleStack stack;Result result;PropagationEngine engine(stack,result,20261011);
        auto* process=new CountedStrict(data);engine.AddInteraction(process);
        auto p=Primary(*data);p.Time.setZ(0.01);stack.AddPrimary(p);engine.Run(variation);
        EXPECT_EQ(process->attempts,1u);EXPECT_EQ(engine.FailureCount(),1u);
    }
}
TEST_F(NuclearPhotonsTest, OffPreservesExistingProductsAndRandomStream) {
    auto data=Load();auto p=Primary(*data);TalysPhotoDisintegration loaded(new LineBackground,data),original(new LineBackground,Load(false));
    Randomizer a(19),b(19);std::vector<Particle> x,y;ASSERT_TRUE(loaded.GetSecondaries(p,x,a));ASSERT_TRUE(original.GetSecondaries(p,y,b));
    ASSERT_EQ(x.size(),y.size());for(size_t i=0;i<x.size();++i){EXPECT_EQ(x[i].Type,y[i].Type);EXPECT_DOUBLE_EQ(x[i].Energy,y[i].Energy);}
    EXPECT_DOUBLE_EQ(a.Rand(),b.Rand());
    auto decay=std::make_shared<NuclearDecayData>(dir.string(),data,true);NuclearDecay off(decay);
    x.clear();ASSERT_TRUE(off.GetSecondaries(p,x,a));ASSERT_EQ(x.size(),2u);
}
}
