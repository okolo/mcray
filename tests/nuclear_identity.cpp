#include "Particle.h"
#include "ParticleStack.h"
#include "PropagationEngine.h"
#include "Output.h"
#include "Logger.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <functional>
#include <unistd.h>

namespace {
using namespace mcray;
void ExpectIdentityError(const std::function<void()>& action) {
    try { action(); FAIL() << "Expected invalid identity error"; }
    catch (Utils::Exception* error) { delete error; }
}
TEST(NuclearIdentity, CopyPreservesMassIdentityAndParticleState) {
    Particle p(56u, 0.01);
    p.SetNucleus(50, 25, 46500 * units.MeV, Particle::NuclearMassModel::TALYS);
    p.Weight = 2.5; p.Energy = 1e20 * units.eV;
    int custom = 42; p.interactionData[0] = &custom;
    Particle copy = p;
    Particle assigned(Photon, 0.); assigned = p;
    for (const auto* q : {&copy, &assigned}) {
        EXPECT_EQ(q->Type, Nucleus); EXPECT_EQ(q->MassNumber(), 50u); EXPECT_EQ(q->AtomicNumber(), 25u);
        EXPECT_EQ(q->Name(), "A50Z25"); EXPECT_EQ(q->ElectricCharge(), 25);
        EXPECT_DOUBLE_EQ(q->Mass(), p.Mass()); EXPECT_EQ(q->MassModel(), p.MassModel());
        EXPECT_EQ(q->id, p.id); EXPECT_EQ(q->SourceParticle, p.SourceParticle);
        EXPECT_EQ(q->interactionData[0], &custom); EXPECT_DOUBLE_EQ(q->Weight, 2.5);
    }
    assigned.SetType(Proton);
    EXPECT_EQ(assigned.MassNumber(), 1u); EXPECT_EQ(assigned.AtomicNumber(), 1u);
    EXPECT_EQ(assigned.MassModel(), Particle::NuclearMassModel::Unset);
    EXPECT_DOUBLE_EQ(assigned.Mass(), Particle::Mass(Proton));
    assigned.SetType(Neutron); EXPECT_EQ(assigned.AtomicNumber(), 0u);
    assigned.SetType(Photon); EXPECT_EQ(assigned.MassNumber(), 0u);
}
TEST(NuclearIdentity, InvalidUpdatesAreAtomicAndUnsetNucleiAreRejected) {
    Particle p(56u, 0.);
    const double mass = p.Mass();
    ExpectIdentityError([&] { p.SetNucleus(2, 3, mass, Particle::NuclearMassModel::TALYS); });
    ExpectIdentityError([&] { p.SetNucleus(2, 1, 0., Particle::NuclearMassModel::TALYS); });
    ExpectIdentityError([&] { p.SetPsbNucleus(1); });
    EXPECT_EQ(p.Name(), "A56Z26"); EXPECT_DOUBLE_EQ(p.Mass(), mass);
    Particle unset(Nucleus, 0.);
    ParticleStack stack(-1);
    ExpectIdentityError([&] { stack.AddPrimary(unset); });
    ExpectIdentityError([&] { stack.AddSecondaries({unset}); });
    ExpectIdentityError([&] { unset.Mass(); });
    ExpectIdentityError([&] { Particle::Mass(Nucleus); });
    p.Type = Proton;
    ExpectIdentityError([&] { p.ValidateIdentity(); });
}
TEST(NuclearOutput, SpectrumAndTrajectoryLogRetainIsotopeIdentity) {
    if (!cosmology.IsInitialized()) cosmology.Init();
    Particle a(56u,0.01), b(56u,0.01);
    a.SetNucleus(50,25,46500*units.MeV,Particle::NuclearMassModel::TALYS);
    b.SetNucleus(50,24,46501*units.MeV,Particle::NuclearMassModel::TALYS);
    a.Energy=b.Energy=1e20*units.eV;
    const auto file=std::filesystem::temp_directory_path()/("mcray-isotope-spectrum-"+std::to_string(getpid()));
    SpectrumOutput spectrum(file.string(),a.Energy,2.);
    spectrum.AddParticle(a); spectrum.AddParticle(b); spectrum.Flush(true);
    std::ifstream input(file); std::string header; std::getline(input,header);
    EXPECT_NE(header.find("A50Z25"),std::string::npos); EXPECT_NE(header.find("A50Z24"),std::string::npos);
    std::ostringstream log; Logger logger(log); logger.log(a); logger.log(b);
    EXPECT_NE(log.str().find("\tA\tZ\tmass/MeV"),std::string::npos);
    EXPECT_NE(log.str().find("\t50\t25\t46500"),std::string::npos);
    EXPECT_NE(log.str().find("\t50\t24\t46501"),std::string::npos);
    std::filesystem::remove(file);
}
TEST(NuclearOutput, RawFlushSeparatesIsotopesAndWritesOneHeaderPerFile) {
    if (!cosmology.IsInitialized()) cosmology.Init();
    const auto dir = std::filesystem::temp_directory_path() / ("mcray-isotope-output-" + std::to_string(getpid()));
    std::filesystem::remove_all(dir);
    {
        RawOutput output(dir.string(), false);
        Particle a(56u, 0.01), b(56u, 0.01);
        a.SetNucleus(50, 25, 46500*units.MeV, Particle::NuclearMassModel::TALYS);
        b.SetNucleus(50, 24, 46501*units.MeV, Particle::NuclearMassModel::TALYS);
        a.Energy = b.Energy = 1e20*units.eV;
        output.AddParticle(a); output.AddParticle(b); output.Flush(false);
        output.AddParticle(a); output.Flush(true);
    }
    for (const auto& name : {"A50Z25", "A50Z24"}) {
        std::ifstream input(dir/name); ASSERT_TRUE(input.good());
        int headers=0, rows=0;
        for (std::string line; std::getline(input,line);) if (!line.empty()) line[0]=='#' ? ++headers : ++rows;
        EXPECT_EQ(headers,1); EXPECT_EQ(rows, std::string(name)=="A50Z25" ? 2 : 1);
    }
    EXPECT_FALSE(std::filesystem::exists(dir/"nucleus"));
    std::filesystem::remove_all(dir);
}
}

namespace {
TEST(NuclearPropagation, FailedCascadeProductsCannotLeakIntoTheNextHistory) {
    using namespace mcray;
    class FailingProcess : public RandomInteraction {
    public:
        double Rate(const Particle& p) const override {
            if (p.Type != Nucleus) return 0.;
            if (p.MassNumber() == 54) Utils::Exception::Throw("Injected daughter-rate failure");
            return 1e6/units.Mpc;
        }
        bool GetSecondaries(Particle& p,std::vector<Particle>& products,Randomizer&) const override {
            for (unsigned a : {55u,54u}) { Particle q=p; q.SetPsbNucleus(a); products.push_back(q); }
            return true;
        }
        RandomInteraction* Clone() const override { return new FailingProcess; }
    };
    if (!cosmology.IsInitialized()) cosmology.Init();
    ParticleStack stack(-1); Result result;
    auto* output=new TotalEnergyOutput; result.AddOutput(output);
    Particle photon(Photon,0.01); photon.Energy=1e15*units.eV; stack.AddPrimary(photon);
    Particle nucleus(56u,0.01); nucleus.Energy=1e20*units.eV;
    stack.AddPrimary(nucleus); stack.AddPrimary(nucleus);
    PropagationEngine engine(stack,result,20261008); engine.AddInteraction(new FailingProcess);
    PropagationEngine clone(engine,20261009); clone.Run();
    EXPECT_EQ(engine.FailureCount(),2u); EXPECT_EQ(clone.FailureCount(),2u);
    EXPECT_DOUBLE_EQ(output->TotalE(),photon.Energy);
    EXPECT_DOUBLE_EQ(output->TotalE(Nucleus),0.);
}
}
