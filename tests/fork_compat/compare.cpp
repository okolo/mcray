// Differential kernel probe; reference implementations live in a separate TU.
#include "GridMagneticField.h"
#include "SpatialGeometry.h"
#include "fork_api.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace mcray;
namespace {
void Fields(const char *file, int mode, const char *queries) {
  ForkUnits(units.Mpc);
  void *reference = ForkField(file, mode);
  Interactions::SavedMF modern(
      file, static_cast<Interactions::SavedMF::FillingType>(mode));
  const CosmoTime time(0);
  std::ifstream input(queries);
  std::string tag;
  long double x[3];
  while (input >> tag >> x[0] >> x[1] >> x[2]) {
    for (auto &coordinate : x)
      coordinate *= units.Mpc;
    double expected[3];
    std::vector<double> actual(3);
    ForkValue(reference, x, expected);
    modern.GetValueGauss(x, time, actual);
    double absolute = 0, magnitude = 1;
    for (int i = 0; i < 3; ++i) {
      absolute = std::max(absolute, std::abs(actual[i] - expected[i]));
      magnitude =
          std::max({magnitude, std::abs(actual[i]), std::abs(expected[i])});
    }
    std::cout << "{\"category\":\"" << tag << "\",\"equal\":"
              << (std::memcmp(expected, actual.data(), sizeof(expected)) == 0
                      ? "true"
                      : "false")
              << ",\"absolute\":" << absolute
              << ",\"relative\":" << absolute / magnitude
              << ",\"scale_fork\":" << ForkScale(reference) / units.Mpc
              << ",\"scale_modern\":"
              << modern.MinVariabilityScale(time) / units.Mpc << "}\n";
  }
  ForkDeleteField(reference);
}
struct Moments {
  size_t n = 0;
  long double mean = 0, m2 = 0;
  void Add(long double value) {
    ++n;
    const auto delta = value - mean;
    mean += delta / n;
    m2 += delta * (value - mean);
  }
  long double Variance() const { return m2 / (n - 1); }
};
constexpr size_t metricCount = 26;
using Metrics = std::array<Moments, metricCount>;
void Observe(Metrics &moments, const long double *direction, double weight,
             const SpatialVector &axis, long double width) {
  // Squared chord avoids cancellation in 1-dot(axis,direction).
  long double chord = 0;
  for (size_t i = 0; i < 3; ++i)
    chord += (direction[i] - axis[i]) * (direction[i] - axis[i]);
  const auto u = chord / (2 * width);
  const auto bin = std::min(19, std::max(0, static_cast<int>(u * 20)));
  moments[0].Add(weight);
  moments[1].Add(weight * u);
  moments[2].Add(weight * u * u);
  for (size_t i = 0; i < 3; ++i)
    moments[3 + i].Add(weight * direction[i]);
  for (int i = 0; i < 20; ++i)
    moments[6 + i].Add(weight * (bin <= i));
}
void Cones(double degrees, double theta, double phi, size_t n) {
  const double angle = degrees * M_PI / 180;
  const SpatialVector axis = SourceGeometry::NormalizeAxis(
      {std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta),
       std::cos(theta)});
  SourceGeometry source({}, axis, SourceGeometry::AngularModel::Cone, degrees);
  Randomizer random(91827);
  void *reference = ForkCone(angle, theta, phi, 17311);
  Metrics old{}, modern{};
  const long double half = degrees * std::acos(-1.L) / 360;
  const long double width = 2 * std::sin(half) * std::sin(half);
  for (size_t i = 0; i < n; ++i) {
    auto sample = ForkNext(reference);
    Observe(old, sample.direction, sample.weight, axis, width);
    Particle particle(Photon, 0);
    particle.Weight = 1;
    source.Initialize(particle, random);
    Observe(modern, particle.Pdir, particle.Weight, axis, width);
  }
  ForkDeleteCone(reference);
  std::cout << "{\"n\":" << n << ",\"metrics\":[";
  for (size_t i = 0; i < metricCount; ++i) {
    if (i)
      std::cout << ',';
    const std::string name = i == 0   ? "weight"
                             : i == 1 ? "polar_mean"
                             : i == 2 ? "polar_second_moment"
                             : i < 6  ? "direction_" + std::to_string(i - 3)
                                      : "polar_cdf_" + std::to_string(i - 5);
    std::cout << "{\"name\":\"" << name << "\",\"mean_fork\":" << old[i].mean
              << ",\"mean_modern\":" << modern[i].mean
              << ",\"variance_fork\":" << old[i].Variance()
              << ",\"variance_modern\":" << modern[i].Variance() << '}';
  }
  std::cout << "]}\n";
}
} // namespace
int main(int argc, char **argv) {
  try {
    std::cout << std::setprecision(17);
    if (argc == 5 && std::string(argv[1]) == "field")
      Fields(argv[2], std::stoi(argv[3]), argv[4]);
    else if (argc == 6 && std::string(argv[1]) == "cone")
      Cones(std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]),
            std::stoul(argv[5]));
    else
      throw std::runtime_error("Invalid probe arguments");
  } catch (Utils::Exception *error) {
    std::cerr << error->Message() << '\n';
    delete error;
    return 1;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
