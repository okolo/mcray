#include "PD.h"
#include "GZK.h"
#include "PPP.h"
#include "ProtonPP.h"
#include "Randomizer.h"
#include "PropagationEngine.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <functional>
#include <cmath>

using namespace mcray;
using namespace Interactions;
using namespace Utils;

namespace {
void Check(bool ok, const char* message) {
    if(!ok) throw std::runtime_error(message);
}
void Near(double a, double b, const char* message) {
    Check(std::abs(a-b) <= 1e-9*std::max(std::abs(a), std::abs(b)), message);
}
void Fails(const std::function<void()>& action, const std::string& text) {
    try { action(); }
    catch(const std::runtime_error& error) {
        Check(std::string(error.what()).find(text) != std::string::npos, "wrong error diagnostic");
        return;
    }
    throw std::runtime_error("expected a model/data error");
}

struct Fixture {
    boost::filesystem::path path;
    std::string oldTables;
    Fixture(): path(boost::filesystem::temp_directory_path()/boost::filesystem::unique_path("mcray-pd-%%%%-%%%%")),
               oldTables(tables_dir) {
        boost::filesystem::create_directories(path/"PD");
        tables_dir = path.string()+"/";
        Write("eps.txt", "# MeV\n\n10\n20\n30\n");
    }
    ~Fixture() { tables_dir = oldTables; boost::filesystem::remove_all(path); }
    void Write(const char* name, const std::string& data) {
        std::ofstream file((path/"PD"/name).string());
        file << data;
    }
};

struct Background : BackgroundIntegral {
    double pRate = 2., nRate = 6.;
    ParticleType lastType = ParticleTypeEOF;
    double lastGamma = 0.;
    double GetRateS(const Function& sigma, const Particle& p, double = 0) override {
        lastType = p.Type;
        lastGamma = p.Energy/p.Mass();
        if(p.Type == Proton) return pRate;
        if(p.Type == Neutron) return nRate;
        return sigma.f(p.Mass()*(p.Mass()+40.*units.MeV));
    }
    double GetRateAndSampleS(const Function& sigma, const Particle& p, Randomizer&, double& s, double = 0) override {
        s = (p.Type == Proton || p.Type == Neutron) ? 2.*units.GeV*units.GeV :
            p.Mass()*(p.Mass()+40.*units.MeV);
        return GetRateS(sigma, p);
    }
    double GetRateK(const Function&, const Particle&, double = 0) override { return 0.; }
    double GetRateAndSampleK(const Function&, const Particle&, Randomizer&, double&, double = 0) override { return 0.; }
    BackgroundIntegral* Clone() const override {
        Background* copy = new Background;
        copy->pRate = pRate;
        copy->nRate = nRate;
        return copy;
    }
};

Particle Primary(ParticleType type) {
    Particle p(type, 0.);
    p.Energy = 1e10*p.Mass();
    return p;
}
void Conserved(const Particle& p, const std::vector<Particle>& products) {
    double Z = 0., A = 0., energy = 0.;
    for(const Particle& sec : products) {
        Z += sec.Weight*sec.ElectricCharge();
        A += sec.Weight*sec.AtomicMass();
        energy += sec.Weight*sec.Energy;
    }
    Near(Z, p.Weight*p.ElectricCharge(), "charge lost");
    Near(A, p.Weight*p.AtomicMass(), "nucleons lost");
    Near(energy, p.Weight*p.Energy, "energy lost");
}
void TestPD() {
    Fixture fixture;
    fixture.Write("xs_pd_sum.txt", "1 1 0 2 0\n3 3 0 2 0\n26 30 0 2 0\n");
    std::string channels = "# channels\n\n1 1 100000 0 1 0\n1 1 10000 0 1 0\n3 3 1 0 1 0\n3 3 1000 0 1 0\n";
    for(int i = 1; i <= 60; ++i) channels += "26 30 " + ToString(i) + " 0 0 0\n";
    channels += "26 30 100000 0 2 0\n";
    fixture.Write("xs_pd_thin.txt", channels);
    PD pd(new Background);
    Randomizer random(7183);
    Particle iron = Primary(Fe56);
    double s;
    Check(pd.SampleS(iron, s, random), "PD sampling failed");
    std::vector<Particle> products;
    pd.SampleSecondaries(iron, products, s, random);
    Check(products.size() == 2 && products[0].Type == Fe55 && products[1].Type == Neutron,
          "channel beyond 50 was not selected");
    Conserved(iron, products);
    Particle hydrogen = Primary(H2);
    pd.SampleS(hydrogen, s, random);
    bool sawProton = false, sawNeutron = false;
    for(int i = 0; i < 100; ++i) {
        products.clear();
        pd.SampleSecondaries(hydrogen, products, s, random);
        Check(products.size() == 2, "residual nucleon lost");
        sawProton = sawProton || products[0].Type == Proton;
        sawNeutron = sawNeutron || products[0].Type == Neutron;
        Conserved(hydrogen, products);
    }
    Check(sawProton && sawNeutron, "residual nucleon channel missing");
    Particle lithium = Primary(Li6);
    lithium.Weight = 3.;
    pd.SampleS(lithium, s, random);
    int alpha = 0;
    for(int i = 0; i < 10000; ++i) {
        products.clear();
        pd.SampleSecondaries(lithium, products, s, random);
        Conserved(lithium, products);
        if(products[0].Type == H2) ++alpha;
    }
    Check(alpha > 4750 && alpha < 5250, "PD channel probabilities are biased");
    // A threshold interval must remain nonzero, even next to a zero grid point.
    s = lithium.Mass()*(lithium.Mass()+30.*units.MeV);
    products.clear();
    pd.SampleSecondaries(lithium, products, s, random);
    Conserved(lithium, products);
    Fails([&] { pd.SampleSecondaries(lithium, products, 0., random); }, "no positive channel");
    Particle missing = Primary(C12);
    Fails([&] { pd.Rate(missing); }, "missing cross sections");
    Particle end(ParticleTypeEOF, 0.);
    Check(pd.Rate(end) == 0. && !pd.SampleS(end, s, random), "sentinel treated as isotope");
    std::unique_ptr<RandomInteraction> copy(pd.Clone());
    boost::filesystem::remove_all(fixture.path/"PD");
    products.clear();
    Check(copy->GetSecondaries(iron, products, random), "clone reread deleted tables");
    Conserved(iron, products);
}
void TestInvalidTables() {
    Fixture f;
    Fails([&] { PD pd(new Background); }, "cannot open");
    f.Write("xs_pd_sum.txt", "3 3 0 1 0\n");
    f.Write("xs_pd_thin.txt", "3 3 100000 0 1 0\n"); // Li5 is absent.
    Fails([&] { PD pd(new Background); }, "unsupported residual");
    f.Write("xs_pd_thin.txt", "3 3 1 0 0 0\n");
    Fails([&] { PD pd(new Background); }, "without channels");
    f.Write("xs_pd_sum.txt", "3 3 1 2\n");
    Fails([&] { PD pd(new Background); }, "missing column");
    f.Write("eps.txt", "10\n10\n");
    Fails([&] { PD pd(new Background); }, "increasing");
}
void TestGZK() {
    Background* background = new Background;
    GZK gzk(background);
    Particle p = Primary(He4);
    Near(gzk.Rate(p), 16., "GZK total rate");
    Randomizer random(8713);
    int protons = 0;
    double s;
    for(int i = 0; i < 20000; ++i) {
        Check(gzk.SampleS(p, s, random), "GZK failed to sample");
        if(s > 0.) ++protons;
    }
    Check(protons > 4700 && protons < 5300, "GZK ignored component rates");
    background->pRate = 0.;
    for(int i = 0; i < 100; ++i)
        Check(gzk.SampleS(p, s, random) && s < 0., "zero-rate proton selected");
    background->nRate = 0.;
    Check(!gzk.SampleS(p, s, random), "zero total rate sampled");
}
void TestPPP() {
    Background* background = new Background;
    PPP ppp(background);
    Particle nucleus = Primary(Fe56), proton = Primary(Proton);
    Near(ppp.Rate(nucleus), 26.*ppp.Rate(proton), "PPP nuclear rate scaling");
    ppp.Rate(nucleus);
    Check(background->lastType == Proton, "PPP rate uses nuclear mass");
    Near(background->lastGamma, nucleus.Energy/nucleus.Mass(), "PPP changed gamma");
    Randomizer r1(77), r2(77);
    double s;
    ppp.SampleS(nucleus, s, r1);
    std::vector<Particle> a, b;
    ppp.SampleSecondariesMean(nucleus, a, s, r1);
    ppp.SampleSecondariesMean(proton, b, s, r2);
    Near(a[0].Energy, b[0].Energy, "PPP electron energy scaled with nuclear energy");
    ProtonPPcel loss(new Background);
    Near(nucleus.Energy*loss.Rate(nucleus), 26.*proton.Energy*loss.Rate(proton),
         "nuclear pair-production energy loss is inconsistent with PPP");
}
struct AlwaysPropagate : IResult {
    bool AbortRequested() const override { return false; }
    bool ContinuePropagation(const Particle&) const override { return true; }
    void Add(std::vector<Particle>) override {}
    cosmo_time GetMinTravelTimeLeft(const Particle&) const override { return 1.; }
};
struct InvalidModel : RandomInteraction {
    double Rate(const Particle&) const override { throw std::runtime_error("invalid nuclear model"); }
    bool GetSecondaries(Particle&, std::vector<Particle>&, Randomizer&) const override { return false; }
    RandomInteraction* Clone() const override { return new InvalidModel; }
};
struct InvalidClone : InvalidModel {
    RandomInteraction* Clone() const override { throw std::runtime_error("invalid nuclear clone"); }
};
void TestFatalErrors() {
    for(int threads : {1, 2}) {
        ParticleStack stack;
        stack.AddPrimary(Primary(He4));
        AlwaysPropagate result;
        PropagationEngine engine(stack, result, 1);
        engine.AddInteraction(new InvalidModel);
        Fails([&] { engine.RunMultithread(0., threads); }, "invalid nuclear model");
    }
    ParticleStack stack;
    AlwaysPropagate result;
    PropagationEngine engine(stack, result, 1);
    engine.AddInteraction(new InvalidClone);
    Fails([&] { engine.RunMultithread(0., 2); }, "invalid nuclear clone");
}
}
int main(int argc, char** argv) {
    try {
        Check(argc == 2, "usage: nuclear_interactions_test TABLES_DIRECTORY");
        tables_dir = std::string(argv[1])+"/";
        Check(Particle::FindNucleus(1, 0) == Proton && Particle::FindNucleus(0, 1) == Neutron,
              "nucleon lookup");
        TestPD(); TestInvalidTables(); TestGZK(); TestPPP(); TestFatalErrors();
        std::cout << "Nuclear interaction regression checks passed\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr << error.what() << '\n';
    } catch(Exception* error) {
        std::cerr << error->Message() << '\n';
        delete error;
    }
    return 1;
}
