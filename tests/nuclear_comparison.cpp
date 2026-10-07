#include <gtest/gtest.h>
#include "gzk_spectral_statistics.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Row { int batch, history; std::string species; double yield, energy; };
struct Options {
    std::string mc, initial, final, report, process = "none";
    std::string gzkReport, gzkGridInitial, gzkGridFinal, gzkStepInitial, gzkStepFinal;
    int primaries = 20, batches = 2, a = 56;
    double energyPerNucleon = 1e19;
} options;
std::vector<Row> rows;
double eventRate = -1., distance = -1.;
bool particleRows = false;

std::string canonical(const std::string& name) {
    static const std::map<std::string, std::string> names{
        {"proton", "p"}, {"neutron", "n"}, {"electron", "e"}, {"positron", "ae"},
        {"photon", "ph"}, {"neutrinoE", "en"}, {"neutrinoAE", "aen"},
        {"neutrinoM", "mn"}, {"neutrinoAM", "amn"}, {"neutrinoT", "tn"}, {"neutrinoAT", "atn"}};
    auto it = names.find(name);
    return it == names.end() ? name : it->second;
}

int mass(const std::string& s) {
    if (s == "proton" || s == "neutron" || s == "p" || s == "n") return 1;
    if (s.size() >= 2 && s[0] == 'A') return std::stoi(s.substr(1));
    return 0;
}
std::string group(const std::string& species) {
    if (species == "p" || species == "proton") return "p";
    if (species == "n" || species == "neutron") return "n";
    const int a = mass(species);
    if (a == 56) return "Fe56";
    if (a >= 41) return "A41_55";
    if (a >= 21) return "A21_40";
    if (a >= 2) return "A02_20";
    return "other";
}

// PrintSpectrum uses E*coefficient*bin_population. Dividing by E restores
// populations up to one shared factor, cancelled against initial injection.
struct ReferenceRow { std::string species; double energy, count; };
std::vector<ReferenceRow> referenceRows(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read reference spectrum " + path);
    std::string line; std::getline(in, line);
    std::regex columns("([0-9]+)=([A-Za-z0-9_]+)");
    std::map<int, std::string> names;
    for (auto it = std::sregex_iterator(line.begin(), line.end(), columns);
         it != std::sregex_iterator(); ++it)
        names[std::stoi((*it)[1]) - 1] = (*it)[2];
    std::vector<ReferenceRow> result;
    while (std::getline(in, line)) {
        std::istringstream row(line);
        std::vector<double> values;
        for (double x; row >> x;) values.push_back(x);
        if (values.empty()) continue;
        if (values[0] <= 0 || names.empty()) throw std::runtime_error("Invalid reference columns");
        for (const auto& entry : names) {
            if (entry.first >= static_cast<int>(values.size())) throw std::runtime_error("Truncated spectrum row");
            if (entry.second != "All") result.push_back({entry.second, values[0], values[entry.first]/values[0]});
        }
    }
    return result;
}
std::map<std::string, double> spectrum(const std::string& path, bool energy = false) {
    std::map<std::string, double> sum;
    for (const auto& row : referenceRows(path)) sum[row.species] += row.count*(energy ? row.energy : 1.);
    return sum;
}

TEST(NuclearEndToEnd, SecondaryMomentsAgreeWithTransportCR) {
    // Pair spectra intentionally use different physical approximations.
    if (options.process != "gzk") GTEST_SKIP();
    const auto initial = spectrum(options.initial), counts = spectrum(options.final), energies = spectrum(options.final, true);
    const double norm = initial.at(options.a == 0 ? "n" : "p");
    ASSERT_GT(norm, 0.);
    ASSERT_GE(options.batches, 20);
    for (const std::string species : {"e", "ae", "ph", "en", "aen", "mn", "amn"}) {
        for (bool energy : {false, true}) {
            std::vector<double> batches(options.batches);
            for (const auto& row : rows) if (canonical(row.species) == species)
                batches.at(row.batch) += (energy ? row.energy/options.energyPerNucleon : row.yield)/options.primaries;
            double mean = 0., variance = 0.;
            for (double value : batches) mean += value/options.batches;
            for (double value : batches) variance += (value-mean)*(value-mean)/(options.batches-1);
            const double se = std::sqrt(variance/options.batches);
            const auto& ref = energy ? energies : counts;
            const auto it = ref.find(species);
            const double expected = it == ref.end() ? 0. : it->second/norm/(energy ? options.energyPerNucleon : 1.);
            std::cout << "secondary " << species << (energy ? " energy" : " count")
                      << " MC=" << mean << " SE=" << se << " TransportCR=" << expected << '\n';
            EXPECT_LE(std::abs(mean-expected), std::max(.001, .01*expected)+4*se);
            if (expected > .01) { EXPECT_LE(se, .01*expected); }
        }
    }
}

TEST(NuclearEndToEnd, PairSecondaryEnergyMatchesPrimaryLoss) {
    if (options.process != "pair") GTEST_SKIP();
    ASSERT_GE(options.batches, 20);
    std::vector<double> loss(options.batches, 1.), emitted(options.batches);
    for (const auto& row : rows) {
        const double fraction = row.energy/options.primaries/options.energyPerNucleon;
        if (canonical(row.species) == "p") loss.at(row.batch) -= fraction;
        if (canonical(row.species) == "e" || canonical(row.species) == "ae")
            emitted.at(row.batch) += fraction;
    }
    double meanLoss = 0., meanEmission = 0., residual = 0., variance = 0.;
    double lossVariance = 0., emissionVariance = 0.;
    for (int i = 0; i < options.batches; ++i) {
        meanLoss += loss[i]/options.batches;
        meanEmission += emitted[i]/options.batches;
        residual += (emitted[i]-loss[i])/options.batches;
    }
    for (int i = 0; i < options.batches; ++i) {
        variance += std::pow(emitted[i]-loss[i]-residual, 2)/(options.batches-1);
        lossVariance += std::pow(loss[i]-meanLoss, 2)/(options.batches-1);
        emissionVariance += std::pow(emitted[i]-meanEmission, 2)/(options.batches-1);
    }
    const double se = std::sqrt(variance/options.batches);
    const double lossSe = std::sqrt(lossVariance/options.batches);
    const double emissionSe = std::sqrt(emissionVariance/options.batches);
    std::cout << "pair energy lost=" << meanLoss << " emitted=" << meanEmission << " SE=" << se << '\n';
    EXPECT_LE(std::abs(residual), .01*meanLoss+4*se);
    EXPECT_LE(se, .01*meanLoss);
    const auto initial = spectrum(options.initial), final = spectrum(options.final, true);
    const double norm = initial.at("p")*options.energyPerNucleon;
    ASSERT_GT(norm, 0.);
    const double refLoss = 1.-final.at("p")/norm;
    const double refEmission = (final.at("e")+final.at("ae"))/norm;
    std::cout << "TransportCR pair energy lost=" << refLoss << " emitted=" << refEmission << '\n';
    EXPECT_NEAR(refEmission, refLoss, .01*refLoss);
    // Compare the transferred energy, rather than the much larger energy
    // retained by the proton. Pair counts and shapes use different models.
    std::cout << "pair cross-code loss SE=" << lossSe << " emitted SE=" << emissionSe << '\n';
    EXPECT_NEAR(meanLoss, refLoss, .01*refLoss+4*lossSe);
    EXPECT_NEAR(meanEmission, refEmission, .01*refLoss+4*emissionSe);
    EXPECT_LE(lossSe, .01*refLoss);
    EXPECT_LE(emissionSe, .01*refLoss);
}

TEST(NuclearEndToEnd, WriteDiagnosticSpectra) {
    if (options.report.empty()) GTEST_SKIP();
    ASSERT_TRUE(particleRows) << "Individual particle rows are required for spectra";
    ASSERT_GE(options.batches, 2);
    const auto initial = spectrum(options.initial);
    const double norm = initial.at(options.a == 0 ? "n" : "p");
    ASSERT_GT(norm, 0.);
    // Fixed centered log bins: ten per decade, six decades below injection.
    // Bins are diagnostic only; sharp features need numerical convergence and
    // a separately justified joint spectral acceptance test.
    auto bin = [](double ratio) { return std::max(-61, std::min(1, static_cast<int>(std::floor(10*std::log10(ratio)+.5)))); };
    using Key = std::pair<std::string, int>;
    std::map<Key, std::vector<double>> mc;
    std::map<Key, double> reference;
    for (const auto& row : rows) {
        ASSERT_GT(row.yield, 0.);
        ASSERT_GT(row.energy, 0.);
        const Key key{canonical(row.species), bin(row.energy/row.yield/options.energyPerNucleon)};
        auto& batches = mc[key];
        batches.resize(options.batches);
        batches.at(row.batch) += row.yield/options.primaries;
    }
    for (const auto& row : referenceRows(options.final)) {
        if (row.count == 0.) continue;
        const Key key{row.species, bin(row.energy/options.energyPerNucleon)};
        reference[key] += row.count/norm;
        mc[key].resize(options.batches);
    }
    std::ofstream out(options.report);
    ASSERT_TRUE(out.good());
    out << std::setprecision(17);
    out << "# Diagnostic counts per injected primary; bin=-61 underflow, bin=1 overflow\n"
        << "# Center energy = injection_energy * 10^(log_bin/10)\n"
        << "# Zero-observation SE is unavailable; the bound concerns history occupancy, not mean multiplicity\n"
        << "# species log_bin MC_mean batch_SE TransportCR zero_history_probability_99pct_upper\n";
    for (const auto& item : mc) {
        double mean = 0., variance = 0.;
        for (double value : item.second) mean += value/options.batches;
        for (double value : item.second) variance += (value-mean)*(value-mean)/(options.batches-1);
        out << item.first.first << '\t' << item.first.second << '\t' << mean << '\t';
        if (mean == 0.) out << "NA";
        else out << std::sqrt(variance/options.batches);
        out << '\t' << reference[item.first] << '\t';
        if (mean == 0.)
            out << -std::expm1(std::log(.01)/(static_cast<double>(options.batches)*options.primaries));
        else out << "NA";
        out << '\n';
    }
}

TEST(NuclearEndToEnd, GzkSpectralPilotReport) {
    if (options.gzkReport.empty()) GTEST_SKIP();
    ASSERT_EQ(options.process, "gzk");
    ASSERT_EQ(options.a, 1);
    ASSERT_DOUBLE_EQ(options.energyPerNucleon, 1e20);
    ASSERT_TRUE(particleRows);
    ASSERT_GE(options.batches, 20);
    const std::set<std::string> species{"ph", "e", "ae", "en", "aen", "mn", "amn"};
    using Key = std::pair<std::string, int>;
    std::map<Key, std::vector<double>> batches;
    std::map<Key, std::set<std::pair<int, int>>> histories;
    for (const auto& s : species) for (int bin = -1; bin <= 4; ++bin)
        batches[{s, bin}].resize(options.batches);
    for (const auto& row : rows) {
        const auto s = canonical(row.species);
        if (!species.count(s)) continue;
        ASSERT_GE(row.batch, 0);
        ASSERT_LT(row.batch, options.batches);
        ASSERT_GE(row.history, 0);
        ASSERT_LT(row.history, options.primaries);
        ASSERT_TRUE(std::isfinite(row.yield) && row.yield > 0.);
        const Key key{s, GzkSpectral::bin(row.energy/row.yield)};
        batches.at(key).at(row.batch) += row.yield/options.primaries;
        histories[key].insert({row.batch, row.history});
    }
    auto reference = [&](const std::string& initialPath, const std::string& finalPath) {
        std::map<Key, double> result;
        const auto initial = spectrum(initialPath);
        const double norm = initial.at("p");
        if (!std::isfinite(norm) || norm <= 0.) throw std::runtime_error("Invalid reference injection");
        for (const auto& row : referenceRows(initialPath)) if (row.count > 0.) {
            if (row.species != "p" || std::abs(row.energy/1e20-1.) > 1e-10)
                throw std::runtime_error("Spectral reference injection differs");
        }
        const auto finalRows = referenceRows(finalPath);
        std::set<double> centers;
        for (const auto& row : finalRows) centers.insert(row.energy);
        const double width = GzkSpectral::gridLogWidth({centers.begin(), centers.end()});
        for (const auto& row : finalRows) {
            if (!species.count(row.species) || row.count == 0.) continue;
            if (!std::isfinite(row.count) || row.count < 0.) throw std::runtime_error("Invalid reference yield");
            const auto fractions = GzkSpectral::cellFractions(row.energy, width);
            for (int region = 0; region < 6; ++region)
                result[{row.species, region-1}] += fractions[region]*row.count/norm;
        }
        return result;
    };
    auto ref = reference(options.initial, options.final);
    const bool controls = !options.gzkGridInitial.empty() && !options.gzkGridFinal.empty() &&
                          !options.gzkStepInitial.empty() && !options.gzkStepFinal.empty();
    auto grid = controls ? reference(options.gzkGridInitial, options.gzkGridFinal) : std::map<Key, double>{};
    auto step = controls ? reference(options.gzkStepInitial, options.gzkStepFinal) : std::map<Key, double>{};
    const double critical = GzkSpectral::critical(options.batches);
    std::ofstream out(options.gzkReport);
    ASSERT_TRUE(out.good());
    out << std::setprecision(17)
        << "# GZK pilot only: no global spectral acceptance is asserted\n"
        << "# Seven species, four fixed bins; alpha=0.01, family=28, budget=0.05\n"
        << "# Bin -1 is unsupported underflow; bin 4 is overflow\n"
        << "# Reference cells use constant dN/dlog(E); counts are shared by log-energy overlap\n"
        << "# Numerical envelope is |reference-grid|+|reference-step|, not a rigorous error bound\n"
        << "# Conditional outcome assumes valid Student batch intervals; occupancy alone does not establish validity\n"
        << "# Student_critical=" << critical << '\n'
        << "# species bin MC_mean SE half_width reference grid step envelope occupied_histories occupied_batches numerical_ready conditional_outcome estimated_history_factor\n";
    for (const auto& item : batches) {
        const Key key = item.first;
        const auto m = GzkSpectral::moments(item.second);
        const double expected = ref[key];
        const double halfWidth = critical*m.se;
        const double envelope = controls ? std::abs(expected-grid[key])+std::abs(expected-step[key]) : 0.;
        int occupied = 0;
        for (double value : item.second) if (value > 0.) ++occupied;
        const bool inside = key.second >= 0 && key.second < 4;
        const bool numerical = inside && controls && expected > 0. && envelope <= .002*expected;
        const bool candidate = inside && occupied == options.batches && m.se > 0.;
        const auto outcome = GzkSpectral::classify(m.mean, expected, halfWidth, envelope, numerical, candidate);
        const double margin = .05*expected-std::abs(m.mean-expected)-envelope;
        out << key.first << '\t' << key.second << '\t' << m.mean << '\t' << m.se << '\t';
        if (m.mean > 0. && m.se > 0.) out << halfWidth; else out << "NA";
        out << '\t' << expected << '\t';
        if (controls) out << grid[key] << '\t' << step[key] << '\t' << envelope;
        else out << "NA\tNA\tNA";
        out << '\t' << histories[key].size() << '\t' << occupied << '\t' << numerical
            << '\t' << GzkSpectral::name(outcome) << '\t';
        if (inside && controls && margin > 0. && m.se > 0.)
            out << std::max(1., std::pow(halfWidth/margin, 2));
        else out << "NA";
        out << '\n';
        if (key.second == 4) {
            EXPECT_DOUBLE_EQ(m.mean, 0.) << "Unexpected MC secondary above injection";
            EXPECT_DOUBLE_EQ(expected, 0.) << "Unexpected reference secondary above injection";
        }
    }
    std::cout << "GZK spectral pilot report: " << options.gzkReport
              << "; conditional intervals only, no global acceptance decision\n";
}

TEST(NuclearEndToEnd, EveryPrimaryHasACompleteHistory) {
    ASSERT_FALSE(rows.empty());
    std::map<int, std::set<int>> histories;
    for (const auto& r : rows) {
        ASSERT_TRUE(std::isfinite(r.yield) && r.yield >= 0);
        ASSERT_TRUE(std::isfinite(r.energy) && r.energy >= 0);
        histories[r.batch].insert(r.history);
    }
    ASSERT_EQ(histories.size(), static_cast<size_t>(options.batches));
    for (const auto& item : histories) EXPECT_EQ(item.second.size(), static_cast<size_t>(options.primaries));
}

TEST(NuclearEndToEnd, PsbPreservesWeightedBaryonNumberPerHistory) {
    if (options.process != "none" && options.process != "psb") GTEST_SKIP();
    std::map<std::pair<int, int>, double> baryons;
    for (const auto& r : rows) baryons[{r.batch, r.history}] += mass(r.species) * r.yield;
    for (const auto& h : baryons) EXPECT_NEAR(h.second, std::max(1, options.a), 1e-9);
}

TEST(NuclearEndToEnd, NoInteractionsPreserveInjection) {
    if (options.process != "none") GTEST_SKIP();
    for (const auto& r : rows) {
        EXPECT_EQ(mass(r.species), std::max(1, options.a));
        EXPECT_DOUBLE_EQ(r.yield, 1.);
        EXPECT_NEAR(r.energy / (std::max(1, options.a) * options.energyPerNucleon), 1., 1e-14);
    }
}

TEST(NuclearEndToEnd, PsbSurvivalMatchesStaticAttenuation) {
    if (options.process != "psb") GTEST_SKIP();
    ASSERT_GE(eventRate, 0.);
    ASSERT_GT(distance, 0.);
    const double expected = std::exp(-eventRate * distance);
    double survivors = 0.;
    for (const auto& r : rows) if (r.species == "A56") survivors += r.yield;
    const double count = static_cast<double>(options.batches) * options.primaries;
    const double observed = survivors / count;
    const double se = std::sqrt(expected * (1.-expected) / count);
    EXPECT_LE(std::abs(observed-expected), 4*se + 1e-4);
}

TEST(NuclearEndToEnd, NucleonEnergyMomentsAgreeWithTransportCR) {
    if (options.process != "gzk" && options.process != "pair-cel" && options.process != "pair")
        GTEST_SKIP();
    ASSERT_FALSE(options.initial.empty());
    ASSERT_FALSE(options.final.empty());
    const auto initial = spectrum(options.initial);
    const auto finalCounts = spectrum(options.final);
    const auto finalEnergy = spectrum(options.final, true);
    const std::string primary = options.a == 0 ? "n" : "p";
    ASSERT_TRUE(initial.count(primary));
    const double norm = initial.at(primary);
    ASSERT_GT(norm, 0.);
    auto canonical = [](const std::string& name) {
        if (name == "proton" || name == "p") return std::string("p");
        if (name == "neutron" || name == "n") return std::string("n");
        return std::string("other");
    };
    // Primary/converted nucleon yields and energy moments. Secondary spectra
    // require the separately planned full binned-output comparison.
    for (const std::string species : {"p", "n"}) {
        for (bool energy : {false, true}) {
            std::vector<double> batch(options.batches);
            for (const auto& r : rows)
                if (canonical(r.species) == species)
                    batch.at(r.batch) += (energy ? r.energy/options.energyPerNucleon : r.yield) / options.primaries;
            const auto& reference = energy ? finalEnergy : finalCounts;
            const auto it = reference.find(species);
            const double expected = it == reference.end() ? 0. : it->second/norm/(energy ? options.energyPerNucleon : 1.);
            double mean = 0., variance = 0.;
            for (double x : batch) mean += x/options.batches;
            if (options.batches > 1)
                for (double x : batch) variance += (x-mean)*(x-mean)/(options.batches-1);
            const double se = std::sqrt(variance/options.batches);
            const double tolerance = std::max(0.001, 0.01*expected);
            std::cout << "comparison " << species << (energy ? " energy" : " count")
                      << " MC=" << mean << " SE=" << se << " TransportCR=" << expected << '\n';
            EXPECT_LE(std::abs(mean-expected), tolerance+4*se);
            if (expected > 0.01) { EXPECT_LE(se, .01*expected); }
        }
    }
}

TEST(NuclearEndToEnd, PsbCompositionAgreesWithTransportCR) {
    if (options.process != "psb") GTEST_SKIP();
    if (options.initial.empty() || options.final.empty()) GTEST_SKIP() << "No TransportCR reference supplied";
    const auto initial = spectrum(options.initial), final = spectrum(options.final);
    ASSERT_EQ(options.a, 56);
    ASSERT_TRUE(initial.count("A56"));
    const double norm = initial.at("A56");
    ASSERT_GT(norm, 0.);
    std::map<std::string, double> expected;
    for (const auto& item : final) expected[group(item.first)] += item.second / norm;
    std::map<std::string, std::vector<double>> batches;
    for (const std::string name : {"Fe56", "A41_55", "A21_40", "A02_20", "p", "n"})
        batches[name].resize(options.batches);
    for (const auto& r : rows) {
        auto found = batches.find(group(r.species));
        if (found != batches.end()) found->second.at(r.batch) += r.yield / options.primaries;
    }
    ASSERT_GE(options.batches, 20) << "Reference comparison requires independent batches";
    // Six prespecified composition groups. 4 standard errors is a conservative
    // two-sided normal cutoff at family-wise alpha=0.01 for six comparisons.
    for (const auto& item : batches) {
        double mean = 0, variance = 0;
        for (double x : item.second) mean += x / options.batches;
        for (double x : item.second) variance += (x-mean)*(x-mean)/(options.batches-1);
        const double se = std::sqrt(variance/options.batches);
        const double ref = expected[item.first];
        std::cout << "comparison " << item.first << " MC=" << mean << " SE=" << se << " TransportCR=" << ref << '\n';
        const double numerical = std::max(0.001, 0.01 * std::abs(ref));
        EXPECT_LE(std::abs(mean-ref), numerical + 4*se) << item.first;
        if (ref > 0.01) { EXPECT_LE(se, 0.01*ref) << item.first << ": insufficient MC precision"; }
    }
}
} // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    try {
        for (int i=1; i<argc; ++i) {
            const std::string key=argv[i];
            if (i+1==argc) throw std::runtime_error("Missing option value");
            const std::string value=argv[++i];
            if (key=="--mc") options.mc=value;
            else if (key=="--initial") options.initial=value;
            else if (key=="--final") options.final=value;
            else if (key=="--spectral-report") options.report=value;
            else if (key=="--gzk-spectral-report") options.gzkReport=value;
            else if (key=="--gzk-grid-initial") options.gzkGridInitial=value;
            else if (key=="--gzk-grid-final") options.gzkGridFinal=value;
            else if (key=="--gzk-step-initial") options.gzkStepInitial=value;
            else if (key=="--gzk-step-final") options.gzkStepFinal=value;
            else if (key=="--process") options.process=value;
            else if (key=="--primaries") options.primaries=std::stoi(value);
            else if (key=="--batches") options.batches=std::stoi(value);
            else if (key=="--a") options.a=std::stoi(value);
            else if (key=="--energy-per-nucleon") options.energyPerNucleon=std::stod(value);
            else throw std::runtime_error("Unknown option "+key);
        }
        if (options.primaries<1 || options.batches<1) throw std::runtime_error("Invalid sample sizes");
        const int controls = !options.gzkGridInitial.empty() + !options.gzkGridFinal.empty() +
                             !options.gzkStepInitial.empty() + !options.gzkStepFinal.empty();
        if (controls != 0 && controls != 4)
            throw std::runtime_error("Provide all four GZK convergence paths or none");
        std::ifstream in(options.mc);
        if (!in) throw std::runtime_error("Cannot read MC output");
        for (std::string line; std::getline(in,line);) {
            if (line.empty()) continue;
            if (line[0]=='#') {
                if (line == "# particle_rows=1") particleRows = true;
                std::smatch match;
                if (std::regex_search(line, match, std::regex("event_rate_per_Mpc=([^ ]+)")))
                    eventRate=std::stod(match[1]);
                if (std::regex_search(line, match, std::regex("L_mpc=([^ ]+)")))
                    distance=std::stod(match[1]);
                continue;
            }
            Row row; std::istringstream fields(line);
            if (!(fields>>row.batch>>row.history>>row.species>>row.yield>>row.energy))
                throw std::runtime_error("Malformed MC row");
            rows.push_back(row);
        }
        return RUN_ALL_TESTS();
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
