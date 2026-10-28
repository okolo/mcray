#pragma once
#include <gsl/gsl_cdf.h>
#include <cmath>
#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

namespace GzkSpectral {
constexpr int familySize = 28;
constexpr double alpha = .01;
constexpr double relativeBudget = .05;

inline int bin(double energy) {
    if (!std::isfinite(energy) || energy <= 0.)
        throw std::invalid_argument("Invalid particle energy");
    if (energy < 1e14) return -1;
    if (energy < 1e17) return 0;
    if (energy < 1e18) return 1;
    if (energy < 1e19) return 2;
    return energy <= 1e20 ? 3 : 4;
}

// Six regions: underflow, four approved intervals, and overflow.
inline std::array<double, 6> cellFractions(double energy, double logWidth) {
    if (!std::isfinite(energy) || energy <= 0. ||
        !std::isfinite(logWidth) || logWidth <= 0.)
        throw std::invalid_argument("Invalid reference cell");
    const double center = std::log10(energy);
    const double low = center-.5*logWidth, high = center+.5*logWidth;
    if (!(high > low)) throw std::invalid_argument("Unresolved reference cell width");
    const std::array<double, 7> edges{low, 14., 17., 18., 19., 20., high};
    std::array<double, 6> fractions{};
    // Constant dN/dlog(E) within each reference cell; include both tails.
    for (int i = 0; i < 6; ++i) {
        const double left = i == 0 ? low : std::max(low, edges[i]);
        const double right = i == 5 ? high : std::min(high, edges[i+1]);
        fractions[i] = std::max(0., right-left)/(high-low);
    }
    return fractions;
}

inline double gridLogWidth(const std::vector<double>& energies) {
    if (energies.size() < 2) throw std::invalid_argument("Reference grid is too short");
    for (double e : energies) if (!std::isfinite(e) || e <= 0.)
        throw std::invalid_argument("Invalid reference grid energy");
    const double span = std::log10(energies.back())-std::log10(energies.front());
    if (!(span > 0.)) throw std::invalid_argument("Reference grid is not increasing");
    const double density = std::round((energies.size()-1)/span);
    if (density < 1.) throw std::invalid_argument("Invalid reference grid density");
    const double width = 1./density;
    for (size_t i = 1; i < energies.size(); ++i) {
        if (energies[i] <= energies[i-1] ||
            std::abs(std::log10(energies[i])-std::log10(energies[i-1])-width) > 1e-5)
            throw std::invalid_argument("Reference grid is not uniform in log energy");
    }
    return width;
}

struct Moments { double mean, se; };
inline Moments moments(const std::vector<double>& batches) {
    if (batches.size() < 2) throw std::invalid_argument("At least two batches are required");
    long double mean = 0., sum = 0.;
    for (double value : batches) {
        if (!std::isfinite(value) || value < 0.)
            throw std::invalid_argument("Invalid batch yield");
        mean += value/batches.size();
    }
    for (double value : batches) sum += (value-mean)*(value-mean);
    return {static_cast<double>(mean),
            static_cast<double>(std::sqrt(sum/(batches.size()-1)/batches.size()))};
}

inline double critical(int batches, int family = familySize, double error = alpha) {
    if (batches < 2 || family < 1 || !std::isfinite(error) || error <= 0. || error >= 1.)
        throw std::invalid_argument("Invalid simultaneous interval settings");
    // Use the upper-tail quantile without subtracting a small probability from one.
    return gsl_cdf_tdist_Qinv(error/(2.*family), batches-1.);
}

enum class Outcome { Agreement, Disagreement, Inconclusive };
inline const char* name(Outcome outcome) {
    return outcome == Outcome::Agreement ? "agreement" :
           outcome == Outcome::Disagreement ? "disagreement" : "inconclusive";
}
inline Outcome classify(double mean, double reference, double halfWidth, double envelope,
                        bool numericalReady, bool statisticalReady) {
    if (!std::isfinite(mean) || !std::isfinite(reference) ||
        !std::isfinite(halfWidth) || !std::isfinite(envelope) ||
        mean < 0. || reference < 0. || halfWidth < 0. || envelope < 0.)
        throw std::invalid_argument("Invalid spectral interval input");
    if (!numericalReady || !statisticalReady || reference == 0.) return Outcome::Inconclusive;
    const double difference = std::abs(mean-reference);
    const double radius = halfWidth+envelope;
    const double budget = relativeBudget*reference;
    if (difference+radius <= budget) return Outcome::Agreement;
    if (difference-radius > budget) return Outcome::Disagreement;
    return Outcome::Inconclusive;
}
} // namespace GzkSpectral
