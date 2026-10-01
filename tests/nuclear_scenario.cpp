#include "PhotoDisintegration.h"
#include "GZK.h"
#include "PPP.h"
#include "ProtonPP.h"
#include "PropagationEngine.h"
#include "Output.h"
#include "frozen_inoue.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace mcray;
using namespace Interactions;

// Each primary owns a complete weighted history, including correlated products.
class HistoryOutput : public Utils::TSmartReferencedObj<IOutput> {
public:
    struct Product { ParticleType type; double weight, energy; };
    struct History {
        std::array<double, ParticleTypeEOF> yield{};
        std::array<double, ParticleTypeEOF> energy{};
        std::vector<Product> products;
    };
    explicit HistoryOutput(bool particles) : particles(particles) {}
    std::map<unsigned long long, History> histories;
    void AddParticle(const Particle& p) override {
        if (!p.SourceParticle) throw std::runtime_error("Missing primary history");
        auto& h = histories[p.SourceParticle->id];
        h.yield[p.Type] += p.Weight;
        h.energy[p.Type] += p.Weight * p.Energy / units.eV;
        if (particles) h.products.push_back({p.Type, p.Weight, p.Weight * p.Energy / units.eV});
    }
    void Flush(bool) override {}
private:
    bool particles;
};

struct Options {
    std::string process = "none", output, background = "cmb";
    int a = 56, batches = 1, primaries = 100;
    unsigned long seed = 20261001;
    double perNucleon = 1e19, distance = 0., accuracy = 1e-4;
    bool scan = false, particles = false;
};

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--scan") { o.scan = true; continue; }
        if (arg == "--particle-output") { o.particles = true; continue; }
        if (i + 1 == argc) throw std::runtime_error("Missing value for " + arg);
        const std::string value = argv[++i];
        if (arg == "--process") o.process = value;
        else if (arg == "--background") o.background = value;
        else if (arg == "--output") o.output = value;
        else if (arg == "--a") o.a = std::stoi(value);
        else if (arg == "--batches") o.batches = std::stoi(value);
        else if (arg == "--primaries") o.primaries = std::stoi(value);
        else if (arg == "--seed") o.seed = std::stoul(value);
        else if (arg == "--energy-per-nucleon") o.perNucleon = std::stod(value);
        else if (arg == "--distance-mpc") o.distance = std::stod(value);
        else if (arg == "--accuracy") o.accuracy = std::stod(value);
        else throw std::runtime_error("Unknown option: " + arg);
    }
    if (o.a < 0 || o.a > 56 || o.batches < 1 || o.primaries < 1 ||
        !std::isfinite(o.perNucleon) || o.perNucleon <= 0 ||
        !std::isfinite(o.distance) || o.distance < 0 ||
        !std::isfinite(o.accuracy) || o.accuracy <= 0 || o.accuracy > 0.01)
        throw std::runtime_error("Invalid scenario parameters");
    if (o.process != "none" && o.process != "psb" && o.process != "gzk" &&
        o.process != "pair-cel" && o.process != "pair")
        throw std::runtime_error("Unknown process");
    if ((o.process == "gzk" && o.a > 1) ||
        ((o.process == "pair-cel" || o.process == "pair") && o.a != 1))
        throw std::runtime_error("Nuclear GZK/pair support is not implemented");
    if (o.process == "psb" && o.a < 2)
        throw std::runtime_error("PSB requires a nucleus");
    if (o.output.empty()) throw std::runtime_error("--output is required");
    if (o.background != "cmb" && o.background != "cmb-inoue")
        throw std::runtime_error("Unknown background");
    return o;
}

int main(int argc, char** argv) {
    try {
        const auto o = parse(argc, argv);
        cosmology.Init(0.1, 0.692, 67.8);
        const double fullDistance = -cosmology.z2t0(0.01) / units.Mpc;
        const double distance = o.distance > 0 ? o.distance : fullDistance;
        if (distance > fullDistance * (1 + 1e-12))
            throw std::runtime_error("Path exceeds the z=0.01 reference interval");
        const ParticleType type = o.a == 0 ? Neutron : o.a == 1 ? Proton :
            static_cast<ParticleType>(A02 + o.a - 2);
        const double temperature = 2.7255 * units.K;
        // Fixed-temperature support already exists in the production background.
        CompoundBackground photons;
        photons.AddComponent(new PlankBackground(temperature, 1e-3 * temperature,
                                                 1e3 * temperature, 0., 0.02, true));
        if (o.background == "cmb-inoue") photons.AddComponent(new FrozenInoueBaseline);
        PBackgroundIntegral background = new ContinuousBackgroundIntegral(
            photons, 0.02, std::pow(10., 0.01), 0.02, o.accuracy);
        PRandomInteraction randomInteraction;
        PCELInteraction cel;
        if (o.process == "psb") randomInteraction = new PhotoDisintegration(background);
        else if (o.process == "gzk") randomInteraction = new GZK(background, static_cast<int>(o.seed));
        else if (o.process == "pair-cel" || o.process == "pair") {
            cel = new ProtonPPcel(background);
            if (o.process == "pair") randomInteraction = new PPP(background);
        }
        std::ofstream out(o.output);
        if (!out) throw std::runtime_error("Cannot create output: " + o.output);
        out << std::setprecision(17);
        out << "# source_z=0.01 L_mpc=" << distance << " T_K=2.7255 OmegaLambda=0.692 H0=67.8"
            << " process=" << o.process << " background=" << o.background
            << " primary=" << Particle::Name(type)
            << " seed=" << o.seed << " accuracy=" << o.accuracy << '\n';
        out << "# particle_rows=" << o.particles << '\n';
        Particle primary(type, 0.01);
        primary.Energy = o.perNucleon * (o.a ? o.a : 1) * units.eV;
        if (o.scan) {
            out << "# E_per_nucleon_eV event_rate_per_Mpc fractional_loss_per_Mpc\n";
            for (int i = 0; i <= 60; ++i) {
                primary.Energy = std::pow(10., 16. + i * 0.1) * (o.a ? o.a : 1) * units.eV;
                out << primary.Energy / units.eV / (o.a ? o.a : 1) << '\t'
                    << (randomInteraction ? randomInteraction->Rate(primary) * units.Mpc : 0.) << '\t'
                    << (cel ? cel->Rate(primary) * units.Mpc : 0.) << '\n';
            }
            return 0;
        }
        const double rate = randomInteraction ? randomInteraction->Rate(primary) * units.Mpc : 0.;
        out << "# E_per_nucleon_eV=" << o.perNucleon << " event_rate_per_Mpc=" << rate << '\n';
        out << "# batch history species yield energy_eV\n";
        for (int batch = 0; batch < o.batches; ++batch) {
            ParticleStack stack;
            Result result;
            CosmoTime end;
            end.setT(std::min<cosmo_time>(0., primary.Time.t() + distance * units.Mpc));
            result.SetEndTime(end);
            auto* output = new HistoryOutput(o.particles);
            result.AddOutput(output);
            PropagationEngine engine(stack, result, o.seed + batch);
            if (randomInteraction) engine.AddInteraction(randomInteraction);
            if (cel) engine.AddInteraction(cel);
            for (int i = 0; i < o.primaries; ++i) {
                Particle injected(type, 0.01);
                injected.Energy = primary.Energy;
                stack.AddPrimary(injected);
            }
            engine.Run(o.process == "psb" || o.process == "none" ? 0. : 1e-3);
            if (output->histories.size() != static_cast<size_t>(o.primaries))
                throw std::runtime_error("Propagation discarded one or more primary histories; comparison is invalid");
            int history = 0;
            for (const auto& item : output->histories) {
                if (o.particles) {
                    for (const auto& product : item.second.products)
                        out << batch << '\t' << history << '\t' << Particle::Name(product.type)
                            << '\t' << product.weight << '\t' << product.energy << '\n';
                } else {
                    for (int t = 0; t < ParticleTypeEOF; ++t)
                        if (item.second.yield[t] != 0.)
                            out << batch << '\t' << history << '\t' << Particle::Name(static_cast<ParticleType>(t))
                                << '\t' << item.second.yield[t] << '\t' << item.second.energy[t] << '\n';
                }
                ++history;
            }
            out << "# batch=" << batch << " primaries=" << o.primaries << '\n';
        }
        return 0;
    } catch (Utils::Exception* error) {
        std::cerr << error->Message() << '\n'; delete error;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; }
    return 1;
}
