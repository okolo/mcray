#include "TalysPhotoDisintegration.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>

namespace Interactions {
using namespace mcray;
using Utils::Exception;

const std::array<TalysData::Identity, 6> TalysData::Emitted = {{{1,0}, {1,1}, {2,1}, {3,1}, {3,2}, {4,2}}};
namespace {
const std::set<TalysData::Identity> approvedTerminals = {
    {5,2}, {5,3}, {6,2}, {8,4}, {9,5}, {12,5}, {16,9}, {20,11},
    {39,20}, {42,21}, {46,23}, {47,24}, {50,25}};
std::vector<std::string> ReadLines(const std::string& file) {
    std::ifstream input(file);
    if (!input) Exception::Throw("Cannot read TALYS data: " + file);
    std::vector<std::string> lines;
    for (std::string line; std::getline(input, line);)
        if (line.find_first_not_of(" \t\r") != std::string::npos && line[0] != '#') lines.push_back(line);
    return lines;
}
TalysData::Identity ReadIdentity(std::istringstream& row) {
    std::string tokens[2];
    unsigned int values[2]{};
    for (int i=0;i<2;++i) {
        if (!(row >> tokens[i]) || tokens[i].find_first_not_of("0123456789") != std::string::npos)
            Exception::Throw("Invalid nuclear identity in TALYS data");
        try {
            const auto value=std::stoull(tokens[i]);
            if (value>std::numeric_limits<unsigned int>::max())
                Exception::Throw("Invalid nuclear identity in TALYS data");
            values[i]=static_cast<unsigned int>(value);
        } catch (const std::exception&) { Exception::Throw("Invalid nuclear identity in TALYS data"); }
    }
    const auto a=values[0],z=values[1];
    if (a < 2 || !z || z > a) Exception::Throw("Invalid nuclear identity in TALYS data");
    return {a,z};
}
void RequireEnd(std::istringstream& row) {
    std::string extra;
    if (row >> extra) Exception::Throw("Unexpected column in TALYS data");
}
std::shared_ptr<const Utils::LinearFunc> ReadSigma(std::istringstream& row, const std::vector<double>& grid) {
    std::vector<double> values(grid.size());
    for (double& value : values) {
        if (!(row >> value) || !std::isfinite(value) || value < 0.)
            Exception::Throw("Invalid or missing TALYS cross-section value");
        value *= 1e-3 * units.barn; // mb -> mcray natural area units
    }
    RequireEnd(row);
    // LinearFunc copies scaled vectors and preserves its zero-handling flag.
    return std::make_shared<Utils::LinearFunc>(grid, values, new Utils::LogScale,
                                               new Utils::LogScale, 0., 0., true);
}
double CheckedRate(double value) {
    if (!std::isfinite(value) || value < 0.) Exception::Throw("Invalid TALYS convolved rate");
    return value;
}
}

TalysData::TalysData(const std::string& directory, bool loadCrossSections) {
    const auto metadata = ReadLines(directory + "/provenance.txt");
    if (metadata.size() < 3 || metadata[0] != "mcray-talys-v1" ||
        metadata[1] != "34b10e27375c0e9894d7c984ac15fb96e12975ef" ||
        metadata[2] != "grid_policy=crpropa-padding")
        Exception::Throw("Unsupported TALYS data format, revision or grid policy");
    for (const auto& line : metadata) fMetadata += "# TALYS " + line + '\n';
    for (const auto& line : ReadLines(directory + "/manifest.json")) fMetadata += "# TALYS manifest " + line + '\n';
    for (const auto& line : ReadLines(directory + "/masses.tsv")) {
        std::istringstream row(line); const auto key = ReadIdentity(row);
        double mass; std::string source;
        if (!(row >> mass >> source) || !std::isfinite(mass) || mass <= 0. ||
            (source != "NIST" && source != "NuDat3")) Exception::Throw("Invalid or unsupported TALYS mass");
        RequireEnd(row);
        if (!fIsotopes.emplace(key, Isotope{}).second) Exception::Throw("Duplicate TALYS mass");
        fIsotopes.at(key).mass = mass * units.MeV;
    }
    fMetadata += std::string("# TALYS load_cross_sections=") + (loadCrossSections ? "On\n" : "Off\n");
    if (!loadCrossSections) return;
    std::vector<double> external, talys;
    for (const auto& label : {"external", "talys"}) {
        auto& grid = std::string(label) == "external" ? external : talys;
        for (const auto& line : ReadLines(directory + "/grid_" + label + ".tsv")) {
            std::istringstream row(line); double k;
            if (!(row >> k) || !std::isfinite(k) || k <= 0. || (!grid.empty() && k*units.MeV <= grid.back()))
                Exception::Throw("Invalid TALYS photon-energy grid");
            RequireEnd(row); grid.push_back(k * units.MeV);
        }
        if (grid.size() != 513) Exception::Throw("TALYS grid must contain 513 padded nodes");
    }
    for (const auto& line : ReadLines(directory + "/terminals.tsv")) {
        std::istringstream row(line); const auto key = ReadIdentity(row); RequireEnd(row);
        if (!approvedTerminals.count(key) || !fIsotopes.count(key) || fIsotopes.at(key).terminal)
            Exception::Throw("Unapproved, massless or duplicate TALYS terminal state");
        fIsotopes.at(key).terminal = true;
        std::cerr << "TALYS: photodisintegration disabled for terminal A" << key.first << "Z" << key.second
                  << "; nuclear decays require a separately enabled process\n";
    }
    for (const auto& line : ReadLines(directory + "/total.tsv")) {
        std::istringstream row(line); const auto key = ReadIdentity(row);
        if (!fIsotopes.count(key)) Exception::Throw("Missing TALYS mass for total cross section");
        auto& isotope = fIsotopes.at(key);
        if (isotope.total || isotope.terminal) Exception::Throw("Duplicate or terminal TALYS total cross section");
        isotope.total = ReadSigma(row, key.first < 12 ? external : talys);
    }
    for (const auto& line : ReadLines(directory + "/channels.tsv")) {
        std::istringstream row(line); const auto key = ReadIdentity(row);
        if (!fIsotopes.count(key)) Exception::Throw("Missing TALYS mass for exclusive cross section");
        auto& isotope = fIsotopes.at(key);
        Channel channel;
        if (isotope.terminal || !(row >> channel.code) || channel.code.size() != 6 ||
            channel.code.find_first_not_of("0123456789") != std::string::npos)
            Exception::Throw("Invalid TALYS channel");
        unsigned int da = 0, dz = 0;
        for (size_t i = 0; i < Emitted.size(); ++i) {
            const auto count = channel.multiplicity[i] = channel.code[i] - '0';
            da += count * Emitted[i].first; dz += count * Emitted[i].second;
        }
        if (!da || da >= key.first || dz > key.second || key.second-dz > key.first-da ||
            (key.first-da > 1 && key.second-dz == 0)) Exception::Throw("Invalid TALYS residual A,Z");
        channel.residual = {key.first-da, key.second-dz};
        for (const auto& existing : isotope.channels)
            if (existing.code == channel.code) Exception::Throw("Duplicate TALYS channel");
        channel.sigma = ReadSigma(row, key.first < 12 ? external : talys);
        isotope.channels.push_back(std::move(channel));
    }
}

const TalysData::Isotope& TalysData::Get(unsigned int a, unsigned int z) const {
    const auto it = fIsotopes.find({a,z});
    if (it == fIsotopes.end()) Exception::Throw("Missing TALYS isotope A" + Utils::ToString(a) + "Z" + Utils::ToString(z));
    return it->second;
}
void TalysData::SetParticle(Particle& p, unsigned int a, unsigned int z) const {
    if (a == 1 && z <= 1) p.SetType(z ? Proton : Neutron);
    else p.SetNucleus(a, z, Get(a,z).mass, Particle::NuclearMassModel::TALYS);
}
std::set<TalysData::Identity> TalysData::ValidateNetwork(unsigned int a, unsigned int z) const {
    std::set<Identity> reached;
    std::vector<Identity> pending{{a,z}};
    while (!pending.empty()) {
        const auto key = pending.back(); pending.pop_back();
        if (!reached.insert(key).second || key.first == 1) continue;
        const auto& isotope = Get(key.first,key.second);
        if (isotope.terminal) continue;
        if (!isotope.total || isotope.channels.empty()) Exception::Throw("Incomplete reachable TALYS cross sections");
        for (const auto& channel : isotope.channels) {
            pending.push_back(channel.residual);
            for (size_t i = 0; i < Emitted.size(); ++i)
                if (channel.multiplicity[i]) pending.push_back(Emitted[i]);
        }
    }
    return reached;
}

TalysPhotoDisintegration::TalysPhotoDisintegration(BackgroundIntegral* background, std::shared_ptr<const TalysData> data)
    : fBackground(background), fData(std::move(data)) {
    if (!background || !fData) Exception::Throw("TALYS requires background integrals and immutable data");
}
RandomInteraction* TalysPhotoDisintegration::Clone() const {
    return new TalysPhotoDisintegration(fBackground->Clone(), fData);
}
const TalysData::Isotope& TalysPhotoDisintegration::Validate(const Particle& p) const {
    p.ValidateIdentity();
    const auto& isotope = fData->Get(p.MassNumber(),p.AtomicNumber());
    if (p.MassModel() != Particle::NuclearMassModel::TALYS || p.Mass() != isotope.mass)
        Exception::Throw("TALYS requires the pinned tabulated nuclear mass");
    if (!isotope.terminal && (!isotope.total || isotope.channels.empty()))
        Exception::Throw("Incomplete TALYS isotope encountered during propagation");
    return isotope;
}
double TalysPhotoDisintegration::Rate(const Particle& p) const {
    if (p.Type != Nucleus) return 0.;
    const auto& isotope = Validate(p);
    return isotope.terminal ? 0. : CheckedRate(fBackground->GetRateK(*isotope.total, p));
}
TalysPhotoDisintegration::Rates TalysPhotoDisintegration::ChannelRates(const Particle& p) const {
    Rates rates;
    if (p.Type != Nucleus) return rates;
    const auto& isotope = Validate(p);
    if (isotope.terminal) return rates;
    rates.total = CheckedRate(fBackground->GetRateK(*isotope.total, p));
    for (const auto& channel : isotope.channels) {
        const double rate = CheckedRate(fBackground->GetRateK(*channel.sigma, p));
        rates.exclusive.push_back(rate); rates.exclusiveSum += rate;
    }
    CheckedRate(rates.exclusiveSum);
    if (rates.total == 0. && rates.exclusiveSum > 0.) Exception::Throw("Positive TALYS exclusive sum with zero total rate");
    if (rates.total > 0. && rates.exclusiveSum <= 0.) Exception::Throw("Positive TALYS total rate with zero exclusive-channel sum");
    return rates;
}
bool TalysPhotoDisintegration::GetSecondaries(Particle& p, std::vector<Particle>& products, Randomizer& random) const {
    const Rates rates = ChannelRates(p);
    if (rates.total <= 0.) return false;
    const auto& isotope = Validate(p);
    double target = random.Rand() * rates.exclusiveSum;
    size_t selected = rates.exclusive.size();
    for (size_t i = 0; i < rates.exclusive.size(); ++i) {
        if (rates.exclusive[i] <= 0.) continue;
        selected = i;
        if (target < rates.exclusive[i]) break;
        target -= rates.exclusive[i];
    }
    if (selected == rates.exclusive.size()) Exception::Throw("Cannot select a positive TALYS channel");
    const auto& channel = isotope.channels[selected];
    // Prepare all products before modifying the caller's output vector.
    std::vector<Particle> event;
    auto emit = [&](const TalysData::Identity& identity) {
        Particle product = p;
        fData->SetParticle(product, identity.first, identity.second);
        product.Energy = p.Energy * identity.first / p.MassNumber();
        event.push_back(product);
    };
    emit(channel.residual);
    for (size_t i = 0; i < TalysData::Emitted.size(); ++i)
        for (unsigned int count = 0; count < channel.multiplicity[i]; ++count) emit(TalysData::Emitted[i]);
    products.insert(products.end(), event.begin(), event.end());
    return true;
}
}
