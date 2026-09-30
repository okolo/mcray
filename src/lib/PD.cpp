#include "PD.h"
#include "TableFunction.h"
#include "Randomizer.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>

using namespace mcray;
using namespace Utils;

namespace Interactions {
namespace {
// Decimal digits: n, p, d, t, He3, He4.
void Decode(int code, int (&counts)[6], int& Z, int& N) {
    const int charges[] = {0, 1, 1, 1, 2, 2};
    const int neutrons[] = {1, 0, 1, 2, 1, 2};
    Z = N = 0;
    for(int i = 5; i >= 0; --i) {
        counts[i] = code % 10;
        code /= 10;
        Z += charges[i]*counts[i];
        N += neutrons[i]*counts[i];
    }
}

bool HasData(std::string& line) {
    line = line.substr(0, line.find('#'));
    return line.find_first_not_of(" \t\r\n") != std::string::npos;
}

std::ifstream Open(const std::string& path) {
    std::ifstream file(path.c_str());
    if(!file) throw std::runtime_error("PD: cannot open " + path);
    return file;
}

std::shared_ptr<const Function> MakeSigma(const std::vector<double>& eps,
                                         const std::vector<double>& sigma,
                                         ParticleType type) {
    std::vector<double> s(eps.size());
    const double mass = Particle::Mass(type);
    for(size_t i = 0; i < s.size(); ++i)
        s[i] = mass*(mass + 2.*eps[i]*units.MeV);
    // Linear sigma on a logarithmic energy axis preserves exact zeroes and
    // additivity of partial cross sections, including threshold intervals.
    return std::shared_ptr<const Function>(
        new LinearFunc(s, sigma, new LogScale(), nullptr, 0., 0.));
}
}

PD::PD(BackgroundIntegral* background): fBackground(background),
    fIsotopes(EndRealNuclei-StartRealNuclei) {
    const std::string directory = tables_dir + "PD/";
    std::vector<double> eps;
    std::string line;
    std::ifstream energies = Open(directory + "eps.txt");
    while(std::getline(energies, line)) {
        if(!HasData(line)) continue;
        std::istringstream row(line);
        double value;
        std::string extra;
        if(!(row >> value) || !std::isfinite(value) || value <= 0. ||
           (!eps.empty() && value <= eps.back()) || (row >> extra))
            throw std::runtime_error("PD: eps.txt must contain increasing positive energies, one per row");
        eps.push_back(value);
    }
    if(eps.size() < 2) throw std::runtime_error("PD: eps.txt requires at least two energies");

    for(int thin = 0; thin != 2; ++thin) {
        const std::string path = directory + (thin ? "xs_pd_thin.txt" : "xs_pd_sum.txt");
        std::ifstream file = Open(path);
        while(std::getline(file, line)) {
            if(!HasData(line)) continue;
            std::istringstream row(line);
            int Z, N, code = 0;
            if(!(row >> Z >> N) || Z < 0 || N < 0 ||
               (thin && (!(row >> code) || code <= 0 || code > 999999)))
                throw std::runtime_error("PD: invalid isotope/channel in " + path);
            std::vector<double> sigma(eps.size());
            bool positive = false;
            for(size_t i = 0; i < sigma.size(); ++i) {
                if(!(row >> sigma[i]) || !std::isfinite(sigma[i]) || sigma[i] < 0.)
                    throw std::runtime_error("PD: invalid cross section or missing column in " + path);
                positive = positive || sigma[i] > 0.;
                sigma[i] *= units.barn*1e-3;
            }
            std::string extra;
            if(row >> extra) throw std::runtime_error("PD: extra column in " + path);
            const ParticleType type = Particle::FindNucleus(Z, N);
            // Tables may contain isotopes outside this particle catalogue.
            if(type < StartRealNuclei || type >= EndRealNuclei) continue;
            Isotope& data = fIsotopes[type-StartRealNuclei];
            if(thin) {
                int counts[6], dZ, dN;
                Decode(code, counts, dZ, dN);
                if(positive && (dZ > Z || dN > N ||
                   ((Z != dZ || N != dN) && Particle::FindNucleus(Z-dZ, N-dN) == ParticleTypeEOF)))
                    throw std::runtime_error("PD: unsupported residual for " +
                        std::string(Particle::Name(type)) + " channel " + ToString(code));
                for(const Channel& channel : data.channels)
                    if(channel.code == code) throw std::runtime_error("PD: duplicate channel in " + path);
                data.channels.push_back(Channel{code, MakeSigma(eps, sigma, type)});
            } else {
                if(data.total) throw std::runtime_error("PD: duplicate isotope in " + path);
                data.total = MakeSigma(eps, sigma, type);
            }
        }
    }
    // A thinned channel table may have a different normalization, but it must
    // cover every interval in which the total cross section is positive.
    for(size_t j = 0; j < fIsotopes.size(); ++j) {
        const Isotope& data = fIsotopes[j];
        if(!data.total) {
            if(!data.channels.empty()) throw std::runtime_error("PD: channel table has no total cross section");
            continue;
        }
        const double mass = Particle::Mass(ParticleType(j+StartRealNuclei));
        for(double energy : eps) {
            const double s = mass*(mass+2.*energy*units.MeV);
            double sum = 0.;
            for(const Channel& channel : data.channels) sum += channel.sigma->f(s);
            if(data.total->f(s) > 0. && sum <= 0.)
                throw std::runtime_error("PD: positive total cross section without channels for " +
                    std::string(Particle::Name(ParticleType(j+StartRealNuclei))));
        }
    }
}

PD::PD(BackgroundIntegral* background, const std::vector<Isotope>& isotopes):
    fBackground(background), fIsotopes(isotopes) {}

const PD::Isotope& PD::Data(ParticleType type) const {
    const Isotope& data = fIsotopes.at(type-StartRealNuclei);
    if(!data.total) throw std::runtime_error("PD: missing cross sections for " +
        std::string(Particle::Name(type)));
    return data;
}

double PD::Rate(const Particle& particle) const {
    if(particle.Type < StartRealNuclei || particle.Type >= EndRealNuclei) return 0.;
    return fBackground->GetRateS(*Data(particle.Type).total, particle);
}

bool PD::SampleS(const Particle& particle, double& s, Randomizer& random) const {
    if(particle.Type < StartRealNuclei || particle.Type >= EndRealNuclei) return false;
    return fBackground->GetRateAndSampleS(*Data(particle.Type).total, particle, random, s) != 0.;
}

RandomInteraction* PD::Clone() const {
    return new PD(fBackground->Clone(), fIsotopes);
}

int PD::Branching(const Isotope& data, double s, Randomizer& random) const {
    std::vector<double> weights;
    double total = 0.;
    int lastPositive = -1;
    for(size_t i = 0; i < data.channels.size(); ++i) {
        const double weight = data.channels[i].sigma->f(s);
        if(!std::isfinite(weight) || weight < 0.)
            throw std::runtime_error("PD: invalid channel cross section");
        weights.push_back(weight);
        total += weight;
        if(weight > 0.) lastPositive = int(i);
    }
    if(!std::isfinite(total) || total <= 0.)
        throw std::runtime_error("PD: no positive channel cross section at sampled energy");
    double target = random.Rand()*total;
    for(size_t i = 0; i < weights.size(); ++i) {
        if(weights[i] > 0. && target < weights[i]) return data.channels[i].code;
        target -= weights[i];
    }
    return data.channels[lastPositive].code; // rounding at the upper endpoint
}

void PD::SampleSecondaries(Particle& particle, std::vector<Particle>& secondaries,
                           double s, Randomizer& random) const {
    if(particle.Type < StartRealNuclei || particle.Type >= EndRealNuclei)
        throw std::runtime_error("PD: primary is not an isotope");
    const int code = Branching(Data(particle.Type), s, random);
    int counts[6], dZ, dN;
    Decode(code, counts, dZ, dN);
    const int Z = particle.ElectricCharge()-dZ;
    const int N = particle.AtomicMass()-particle.ElectricCharge()-dN;
    const ParticleType residual = Particle::FindNucleus(Z, N);
    if(Z < 0 || N < 0 || ((Z || N) && residual == ParticleTypeEOF))
        throw std::runtime_error("PD: unsupported residual nucleus");
    if(Z || N) {
        Particle sec = particle;
        sec.Type = residual;
        sec.Energy *= sec.Mass()/particle.Mass();
        sec.fCascadeProductionTime = particle.Time;
        secondaries.push_back(sec);
    }
    const ParticleType types[] = {Neutron, Proton, H2, H3, He3, He4};
    for(int i = 0; i < 6; ++i) {
        if(!counts[i]) continue;
        Particle sec = particle;
        sec.Type = types[i];
        sec.Weight *= counts[i];
        sec.Energy *= sec.Mass()/particle.Mass();
        sec.fCascadeProductionTime = particle.Time;
        secondaries.push_back(sec);
    }
}
}
