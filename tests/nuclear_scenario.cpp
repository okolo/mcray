#include "PhotoDisintegration.h"
#include "TalysPhotoDisintegration.h"
#include "NuclearDecay.h"
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
    struct Product { std::string name; double weight, energy; };
    struct History {
        std::map<std::string, double> yield;
        std::map<std::string, double> energy;
        std::vector<Product> products;
    };
    explicit HistoryOutput(bool particles) : particles(particles) {}
    std::map<unsigned long long, History> histories;
    void AddParticle(const Particle& p) override {
        if (!p.SourceParticle) throw std::runtime_error("Missing primary history");
        auto& h = histories[p.SourceParticle->id];
        h.yield[p.Name()] += p.Weight;
        h.energy[p.Name()] += p.Weight * p.Energy / units.eV;
        if (particles) h.products.push_back({p.Name(), p.Weight, p.Weight * p.Energy / units.eV});
    }
    void Flush(bool) override {}
private:
    bool particles;
};

struct Options {
    std::string process = "none", output, background = "cmb", talysData;
    int z = 0, a = 56, batches = 1, primaries = 100;
    unsigned long seed = 20261001;
    double perNucleon = 1e19, distance = 0., accuracy = 1e-4;
    double backgroundStepDecades = 0.01;
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
        else if (arg == "--talys-data") o.talysData = value;
        else if (arg == "--background") o.background = value;
        else if (arg == "--output") o.output = value;
        else if (arg == "--a") o.a = std::stoi(value);
        else if (arg == "--charge") o.z = std::stoi(value);
        else if (arg == "--batches") o.batches = std::stoi(value);
        else if (arg == "--primaries") o.primaries = std::stoi(value);
        else if (arg == "--seed") o.seed = std::stoul(value);
        else if (arg == "--energy-per-nucleon") o.perNucleon = std::stod(value);
        else if (arg == "--distance-mpc") o.distance = std::stod(value);
        else if (arg == "--accuracy") o.accuracy = std::stod(value);
        else if (arg == "--background-step-decades") o.backgroundStepDecades = std::stod(value);
        else throw std::runtime_error("Unknown option: " + arg);
    }
    if (o.a < 0 || o.a > 56 || o.z < 0 || o.z > o.a || o.batches < 1 || o.primaries < 1 ||
        !std::isfinite(o.perNucleon) || o.perNucleon <= 0 ||
        !std::isfinite(o.distance) || o.distance < 0 ||
        !std::isfinite(o.accuracy) || o.accuracy <= 0 || o.accuracy > 0.01 ||
        !std::isfinite(o.backgroundStepDecades) || o.backgroundStepDecades <= 0 ||
        o.backgroundStepDecades > 0.1)
        throw std::runtime_error("Invalid scenario parameters");
    if (o.process != "none" && o.process != "psb" && o.process != "gzk" &&
        o.process != "pair-cel" && o.process != "pair" && o.process != "talys" &&
        o.process != "decay" && o.process != "talys-decay")
        throw std::runtime_error("Unknown process");
    if ((o.process == "gzk" && o.a > 1) ||
        ((o.process == "pair-cel" || o.process == "pair") && o.a != 1))
        throw std::runtime_error("Nuclear GZK/pair support is not implemented");
    if (o.process == "psb" && o.a < 2)
        throw std::runtime_error("PSB requires a nucleus");
    if ((o.process == "talys" || o.process == "talys-decay") && (o.a != 56 || (o.z != 0 && o.z != 26) || o.talysData.empty()))
        throw std::runtime_error("TALYS requires Fe56 and --talys-data");
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
            Nucleus;
        const double temperature = 2.7255 * units.K;
        // Fixed-temperature support already exists in the production background.
        CompoundBackground photons;
        photons.AddComponent(new PlankBackground(temperature, 1e-3 * temperature,
                                                 1e3 * temperature, 0., 0.02, true));
        if (o.background == "cmb-inoue") photons.AddComponent(new FrozenInoueBaseline);
        PBackgroundIntegral background = new ContinuousBackgroundIntegral(
            photons, 0.02, std::pow(10., o.backgroundStepDecades), 0.02, o.accuracy);
        std::shared_ptr<const TalysData> talys;
        std::shared_ptr<const NuclearDecayData> decay;
        if (o.process == "talys" || o.process == "decay" || o.process == "talys-decay") {
            talys = std::make_shared<TalysData>(o.talysData);
            if (o.a<2 || o.talysData.empty()) throw std::runtime_error("Nuclear decays require --a >= 2 and --talys-data");
            if (o.process != "talys") {
                decay=std::make_shared<NuclearDecayData>(o.talysData,talys);
                decay->ValidateNetwork(o.a,(o.z ? o.z : CNucleus::getZ(o.a)),o.process=="talys-decay");
            } else talys->ValidateNetwork(56,26);
        }
        auto configure = [&](Particle& p) {
            if (type != Nucleus) return;
            if (talys) talys->SetParticle(p,o.a,(o.z ? o.z : CNucleus::getZ(o.a)));
            else p.SetPsbNucleus(o.a);
        };
        PRandomInteraction randomInteraction;
        PCELInteraction cel;
        if (o.process == "psb") randomInteraction = new PhotoDisintegration(background);
        else if (o.process == "talys" || o.process == "talys-decay") randomInteraction = new TalysPhotoDisintegration(background, talys);
        else if (o.process == "gzk") randomInteraction = new GZK(background, static_cast<int>(o.seed));
        else if (o.process == "pair-cel" || o.process == "pair") {
            cel = new ProtonPPcel(background);
            if (o.process == "pair") randomInteraction = new PPP(background);
        }
        std::ofstream out(o.output);
        if (!out) throw std::runtime_error("Cannot create output: " + o.output);
        out << std::setprecision(17);
        Particle primary(type, 0.01); configure(primary);
        primary.Energy = o.perNucleon * (o.a ? o.a : 1) * units.eV;
        out << "# source_z=0.01 L_mpc=" << distance << " T_K=2.7255 OmegaLambda=0.692 H0=67.8"
            << " process=" << o.process << " background=" << o.background
            << " primary=" << primary.Name()
            << " seed=" << o.seed << " accuracy=" << o.accuracy
            << " background_step_decades=" << o.backgroundStepDecades << '\n';
        out << "# particle_rows=" << o.particles << '\n';
        if (talys) out << talys->Metadata();
        if (decay) out << decay->Metadata();
        if (o.process == "decay") randomInteraction=new NuclearDecay(decay);
        if (o.scan) {
            out << "# E_per_nucleon_eV event_rate_per_Mpc fractional_loss_per_Mpc";
            if (talys && o.process != "decay") out << " exclusive_sum_per_Mpc exclusive_sum_over_total";
            out << '\n';
            for (int i = 0; i <= 60; ++i) {
                primary.Energy = std::pow(10., 16. + i * 0.1) * (o.a ? o.a : 1) * units.eV;
                out << primary.Energy / units.eV / (o.a ? o.a : 1) << '\t'
                    << (randomInteraction ? randomInteraction->Rate(primary) * units.Mpc : 0.) << '\t'
                    << (cel ? cel->Rate(primary) * units.Mpc : 0.);
                if (talys && o.process != "decay") {
                    const auto rates=static_cast<TalysPhotoDisintegration*>(static_cast<RandomInteraction*>(randomInteraction))->ChannelRates(primary);
                    out << '\t' << rates.exclusiveSum*units.Mpc << '\t'
                        << (rates.total > 0. ? rates.exclusiveSum/rates.total : 0.);
                }
                out << '\n';
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
            if (o.process == "talys-decay") engine.AddInteraction(new NuclearDecay(decay));
            for (int i = 0; i < o.primaries; ++i) {
                Particle injected(type, 0.01);
                configure(injected);
                injected.Energy = primary.Energy;
                stack.AddPrimary(injected);
            }
            engine.Run(o.process == "psb" || o.process == "talys" || o.process == "decay" || o.process == "talys-decay" || o.process == "none" ? 0. : 1e-3);
            if (engine.FailureCount() || output->histories.size() != static_cast<size_t>(o.primaries))
                throw std::runtime_error("Propagation discarded one or more primary histories; comparison is invalid");
            int history = 0;
            for (const auto& item : output->histories) {
                if (o.particles) {
                    for (const auto& product : item.second.products)
                        out << batch << '\t' << history << '\t' << product.name
                            << '\t' << product.weight << '\t' << product.energy << '\n';
                } else {
                    for (const auto& entry : item.second.yield)
                        if (entry.second != 0.)
                            out << batch << '\t' << history << '\t' << entry.first
                                << '\t' << entry.second << '\t' << item.second.energy.at(entry.first) << '\n';
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
