#include "NuclearDecay.h"
#include "PropagationEngine.h"
#include "Output.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <unistd.h>

namespace {
using namespace mcray;
using namespace Interactions;
class DecayTest : public ::testing::Test {
protected:
    std::filesystem::path dir;
    void Write(const std::string& name,const std::string& text) { std::ofstream(dir/name)<<text; }
    void SetUp() override {
        dir=std::filesystem::temp_directory_path()/("mcray-decay-"+std::to_string(getpid()));
        std::filesystem::remove_all(dir); std::filesystem::create_directory(dir);
        Write("provenance.txt","mcray-talys-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\ngrid_policy=crpropa-padding\n");
        Write("manifest.json","{}\n");
        std::string grid; for (int i=1;i<=513;++i) grid+=std::to_string(i)+'\n';
        Write("grid_talys.tsv",grid); Write("grid_external.tsv",grid);
        Write("masses.tsv","56 26 52194 NIST\n56 27 52192 NIST\n56 25 52192 NIST\n55 26 51260 NIST\n8 4 7456 NIST\n4 2 3728 NIST\n");
        Write("terminals.tsv",""); Write("total.tsv",""); Write("channels.tsv","");
        Write("decay_provenance.txt","mcray-decay-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\n");
        Write("decay_manifest.json","{}\n");
        Write("decay_stable.tsv","56 27\n56 25\n55 26\n4 2\n");
        Write("decays.tsv","56 26 10000 4\n56 26 00001 1.3333333333333333\n8 4 00100 1e-16\n");
    }
    void TearDown() override { std::filesystem::remove_all(dir); }
    std::shared_ptr<const NuclearDecayData> Load() {
        return std::make_shared<NuclearDecayData>(dir.string(),std::make_shared<TalysData>(dir.string()));
    }
    void ExpectError(const std::function<void()>& action,const std::string& message) {
        try { action(); FAIL()<<"Expected error: "<<message; }
        catch (Utils::Exception* error) { std::string text=error->Message(); delete error; EXPECT_NE(text.find(message),std::string::npos)<<text; }
    }
};
TEST_F(DecayTest, ProperLifetimeGammaScalingStableStatesAndClones) {
    auto data=Load(); NuclearDecay process(data);
    Particle p(Nucleus,0.); data->Nuclei().SetParticle(p,56,26); p.Energy=1e8*p.Mass();
    EXPECT_NEAR(process.Rate(p)*units.second*1e8,1.,1e-14);
    p.Time.setZ(0.5); EXPECT_NEAR(process.Rate(p)*units.second*1e8,1.,1e-14);
    p.Energy*=2.; EXPECT_NEAR(process.Rate(p)*units.second*2e8,1.,1e-14);
    std::unique_ptr<RandomInteraction> clone(process.Clone()); EXPECT_DOUBLE_EQ(process.Rate(p),clone->Rate(p));
    data->Nuclei().SetParticle(p,56,27); EXPECT_DOUBLE_EQ(process.Rate(p),0.);
    p.SetType(Neutron); EXPECT_DOUBLE_EQ(process.Rate(p),0.);
    p.SetPsbNucleus(56); ExpectError([&]{process.Rate(p);},"tabulated nuclear mass");
}
TEST_F(DecayTest, BranchSamplingChargeBaryonsWeightsAndBetaGamma) {
    auto data=Load(); NuclearDecay process(data);
    Particle p(Nucleus,0.); data->Nuclei().SetParticle(p,56,26); p.Energy=1e8*p.Mass(); p.Weight=2.5;
    Randomizer random(20261008); unsigned beta=0;
    for (int i=0;i<20000;++i) {
        std::vector<Particle> products; ASSERT_TRUE(process.GetSecondaries(p,products,random));
        unsigned baryons=0; double charge=0.; bool hasElectron=false,hasNeutrino=false;
        for (const auto& product:products) {
            baryons+=product.MassNumber(); charge+=product.ElectricCharge();
            EXPECT_DOUBLE_EQ(product.Weight,p.Weight); EXPECT_TRUE(std::isfinite(product.Energy)); EXPECT_GE(product.Energy,0.);
            hasElectron|=product.Type==Electron; hasNeutrino|=product.Type==NeutrinoAE;
        }
        EXPECT_EQ(baryons,56u); EXPECT_DOUBLE_EQ(charge,26.);
        if (hasElectron) {
            ++beta; EXPECT_TRUE(hasNeutrino); EXPECT_EQ(products.back().AtomicNumber(),27u);
            EXPECT_NEAR(products.back().Energy/products.back().Mass(),1e8,1e-6);
        } else {
            double sum=0.; for(const auto& product:products) sum+=product.Energy;
            EXPECT_NEAR(sum/p.Energy,1.,1e-14);
        }
    }
    EXPECT_NEAR(beta/20000.,0.25,0.015); EXPECT_TRUE(process.KillsPrimary());
}
TEST_F(DecayTest, BetaPlusProducesPositronAndElectronNeutrino) {
    Write("decays.tsv","56 26 01000 2\n8 4 00100 1e-16\n");
    auto data=Load(); NuclearDecay process(data);
    Particle p(Nucleus,0.); data->Nuclei().SetParticle(p,56,26); p.Energy=1e8*p.Mass();
    Randomizer random(17); std::vector<Particle> products; ASSERT_TRUE(process.GetSecondaries(p,products,random));
    ASSERT_EQ(products.size(),3u); EXPECT_EQ(products[0].Type,Positron); EXPECT_EQ(products[1].Type,NeutrinoE);
    EXPECT_EQ(products[2].AtomicNumber(),25u); EXPECT_DOUBLE_EQ(products[2].ElectricCharge()+products[0].ElectricCharge(),26.);
}
TEST_F(DecayTest, ConditionalKernelPreservesSeededProductsAndRates) {
    auto data=Load(); NuclearDecay process(data);
    Particle p(Nucleus,0.); data->Nuclei().SetParticle(p,56,26); p.Energy=1e8*p.Mass();p.Weight=2.5;
    const double before=process.Rate(p);
    Randomizer ordinary(20261024),conditional(20261024);
    const auto& channels=data->Get(56,26);
    double sum=0.;for(const auto& c:channels) sum+=c.restRate;
    for(int i=0;i<100;++i) {
        // Consume precisely the branch draw used by ordinary production.
        double u=conditional.Rand()*sum;
        const auto* selected=&channels.back();
        for(const auto& c:channels) {u-=c.restRate;if(u<0.) {selected=&c;break;}}
        std::vector<Particle> a,b;
        ASSERT_TRUE(process.GetSecondaries(p,a,ordinary));
        ASSERT_TRUE(process.GetChannelSecondaries(p,selected->code,b,conditional));
        ASSERT_EQ(a.size(),b.size());
        for(size_t j=0;j<a.size();++j) {
            EXPECT_EQ(a[j].Name(),b[j].Name());
            EXPECT_DOUBLE_EQ(a[j].Energy,b[j].Energy);
            EXPECT_DOUBLE_EQ(a[j].Weight,b[j].Weight);
        }
    }
    EXPECT_DOUBLE_EQ(before,process.Rate(p));
    std::vector<Particle> products;
    ExpectError([&]{process.GetChannelSecondaries(p,"99999",products,conditional);},"Unknown conditional");
    data->Nuclei().SetParticle(p,56,27);
    ExpectError([&]{process.GetChannelSecondaries(p,"10000",products,conditional);},"unstable nuclear state");
}
TEST_F(DecayTest, AlphaEmissionIncludesBothIdenticalDaughters) {
    auto data=Load(); NuclearDecay process(data);
    Particle p(Nucleus,0.); data->Nuclei().SetParticle(p,8,4); p.Energy=8e19*units.eV;
    Randomizer random(17); std::vector<Particle> products; ASSERT_TRUE(process.GetSecondaries(p,products,random));
    ASSERT_EQ(products.size(),2u);
    for (const auto& product:products) { EXPECT_EQ(product.Name(),"A04Z02"); EXPECT_DOUBLE_EQ(product.Energy,0.5*p.Energy); }
    EXPECT_EQ(data->ValidateNetwork(8,4,false).size(),2u);
}
TEST_F(DecayTest, TransportSurvivalMatchesDilatedExponentialLifetime) {
    if (!cosmology.IsInitialized()) cosmology.Init(0.1);
    auto data=Load();
    class Counts : public Utils::TSmartReferencedObj<IOutput> {
    public:
        unsigned survived=0; double baryons=0.,charge=0.;
        void AddParticle(const Particle& p) override {
            if (p.Type==Nucleus && p.MassNumber()==56 && p.AtomicNumber()==26) ++survived;
            baryons+=p.Weight*p.MassNumber(); charge+=p.Weight*p.ElectricCharge();
        }
        void Flush(bool) override {}
    };
    ParticleStack stack; Result result;
    Particle primary(Nucleus,0.01); data->Nuclei().SetParticle(primary,56,26); primary.Energy=1e8*primary.Mass();
    CosmoTime end; end.setT(primary.Time.t()+std::log(2.)*1e8*units.second); result.SetEndTime(end);
    auto* counts=new Counts; result.AddOutput(counts);
    PropagationEngine engine(stack,result,20261009); engine.AddInteraction(new NuclearDecay(data));
    for (int i=0;i<10000;++i) stack.AddPrimary(primary);
    engine.Run(0.);
    EXPECT_EQ(engine.FailureCount(),0u);
    EXPECT_NEAR(counts->survived,5000.,250.); // Five binomial standard deviations.
    EXPECT_DOUBLE_EQ(counts->baryons,560000.); EXPECT_DOUBLE_EQ(counts->charge,260000.);
}
TEST_F(DecayTest, CompoundBetaProtonChannelUsesSequentialKinematics) {
    Write("decays.tsv","56 26 10010 2\n8 4 00100 1e-16\n");
    auto data=Load(); NuclearDecay process(data);
    Particle p(Nucleus,0.); data->Nuclei().SetParticle(p,56,26); p.Energy=1e8*p.Mass(); p.Weight=3.;
    Randomizer random(31); std::vector<Particle> products; ASSERT_TRUE(process.GetSecondaries(p,products,random));
    ASSERT_EQ(products.size(),4u); EXPECT_EQ(products[0].Type,Electron); EXPECT_EQ(products[1].Type,NeutrinoAE);
    EXPECT_EQ(products[2].Type,Proton); EXPECT_EQ(products[3].Name(),"A55Z26");
    const double afterBeta=1e8*52192*units.MeV;
    EXPECT_NEAR(products[2].Energy/afterBeta,1./56.,1e-14);
    EXPECT_NEAR(products[3].Energy/afterBeta,55./56.,1e-14);
    double charge=0.; unsigned a=0;
    for(const auto& product:products) { charge+=product.ElectricCharge(); a+=product.MassNumber(); EXPECT_DOUBLE_EQ(product.Weight,3.); }
    EXPECT_EQ(a,56u); EXPECT_DOUBLE_EQ(charge,26.);
    EXPECT_EQ(p.Name(),"A56Z26");
}
TEST_F(DecayTest, PreflightRejectsMissingStableStatesAndInvalidEndpoints) {
    auto data=Load(); EXPECT_EQ(data->ValidateNetwork(56,26,false).size(),4u);
    ExpectError([&]{data->ValidateNetwork(56,26,true);},"photodisintegration data");
    Write("decay_stable.tsv","56 27\n56 25\n4 2\n"); data=Load();
    ExpectError([&]{data->ValidateNetwork(56,26,false);},"Missing nuclear decay");
    Write("decays.tsv","56 26 99999 1\n"); ExpectError([&]{Load();},"Invalid residual");
    Write("decays.tsv","56 27 01000 1\n"); Write("decay_stable.tsv","");
    ExpectError([&]{Load();},"Nonpositive nuclear beta");
    Write("decays.tsv","56 26 10000 nan\n"); ExpectError([&]{Load();},"partial lifetime");
}
}
