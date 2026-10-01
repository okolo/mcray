#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
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
    std::string mc, initial, final, process = "none";
    int primaries = 20, batches = 2, a = 56;
    double energyPerNucleon = 1e19;
} options;
std::vector<Row> rows;
double eventRate = -1., distance = -1.;

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
std::map<std::string, double> spectrum(const std::string& path, bool energy = false) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read reference spectrum " + path);
    std::string line; std::getline(in, line);
    std::regex columns("([0-9]+)=([A-Za-z0-9_]+)");
    std::map<int, std::string> names;
    for (auto it = std::sregex_iterator(line.begin(), line.end(), columns);
         it != std::sregex_iterator(); ++it)
        names[std::stoi((*it)[1]) - 1] = (*it)[2];
    std::map<std::string, double> sum;
    while (std::getline(in, line)) {
        std::istringstream row(line);
        std::vector<double> values;
        for (double x; row >> x;) values.push_back(x);
        if (values.empty()) continue;
        if (values[0] <= 0 || names.empty()) throw std::runtime_error("Invalid reference columns");
        for (const auto& entry : names) {
            if (entry.first >= static_cast<int>(values.size())) throw std::runtime_error("Truncated spectrum row");
            if (entry.second != "All") sum[entry.second] += values[entry.first] / (energy ? 1. : values[0]);
        }
    }
    return sum;
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
    for (const auto& h : baryons) EXPECT_NEAR(h.second, options.a, 1e-9);
}

TEST(NuclearEndToEnd, NoInteractionsPreserveInjection) {
    if (options.process != "none") GTEST_SKIP();
    for (const auto& r : rows) {
        EXPECT_EQ(mass(r.species), options.a);
        EXPECT_DOUBLE_EQ(r.yield, 1.);
        EXPECT_NEAR(r.energy / (options.a * options.energyPerNucleon), 1., 1e-14);
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
            if (expected > 0.01) EXPECT_LE(se, .01*expected);
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
        if (ref > 0.01) EXPECT_LE(se, 0.01*ref) << item.first << ": insufficient MC precision";
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
            else if (key=="--process") options.process=value;
            else if (key=="--primaries") options.primaries=std::stoi(value);
            else if (key=="--batches") options.batches=std::stoi(value);
            else if (key=="--a") options.a=std::stoi(value);
            else if (key=="--energy-per-nucleon") options.energyPerNucleon=std::stod(value);
            else throw std::runtime_error("Unknown option "+key);
        }
        if (options.primaries<1 || options.batches<1) throw std::runtime_error("Invalid sample sizes");
        std::ifstream in(options.mc);
        if (!in) throw std::runtime_error("Cannot read MC output");
        for (std::string line; std::getline(in,line);) {
            if (line.empty()) continue;
            if (line[0]=='#') {
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
