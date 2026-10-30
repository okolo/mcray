#include "SpatialGeometry.h"
#include <algorithm>
#include <array>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <functional>
#include <gtest/gtest.h>
#include <limits>
#include <memory>
#include <tuple>

namespace {
using namespace mcray;
using Model = SourceGeometry::AngularModel;
using High = boost::multiprecision::cpp_dec_float_100;

void Fatal(const std::function<void()>& action) {
  try { action(); FAIL() << "Expected fatal geometry error"; }
  catch (Utils::Exception* error) {
    std::unique_ptr<Utils::Exception> owner(error);
    EXPECT_FALSE(error->Retryable());
  }
}
coord_type Width(double degrees) {
  const coord_type half = degrees * std::acos(coord_type(-1)) / 360;
  return 2 * std::sin(half) * std::sin(half);
}
coord_type Offset(const Particle& particle, const SpatialVector& axis) {
  // The squared chord avoids cancellation in 1-dot for narrow profiles.
  coord_type squared = 0;
  for (size_t i = 0; i < 3; ++i) {
    const auto difference = particle.Pdir[i] - axis[i];
    squared += difference * difference;
  }
  return squared / 2;
}
TEST(VmfSource, ZeroConcentrationExactlyMatchesConeAndRandomStream) {
  for (double angle : {0., .001, 10., 180.}) {
    SourceGeometry cone({1, 2, 3}, {1, -2, 3}, Model::Cone, angle);
    SourceGeometry vmf({1, 2, 3}, {1, -2, 3}, Model::Vmf, angle, 0);
    Randomizer first(401715), second(401715);
    Particle a(56u, 0), b = a;
    a.Energy = b.Energy = 1e12 * units.eV;
    a.Weight = b.Weight = 2.5;
    a.fPrevId = b.fPrevId = 81;
    a.Ninteractions = b.Ninteractions = 9;
    const auto id = a.id;
    for (int draw = 0; draw < 100; ++draw) {
      cone.Initialize(a, first);
      vmf.Initialize(b, second);
      for (size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(a.Pdir[i], b.Pdir[i]);
        EXPECT_EQ(a.X[i], b.X[i]);
      }
      EXPECT_EQ(b.id, id);
      EXPECT_EQ(b.fPrevId, 81u);
      EXPECT_EQ(b.Ninteractions, 9);
      EXPECT_EQ(b.MassNumber(), 56u);
      EXPECT_EQ(b.AtomicNumber(), 26u);
      EXPECT_EQ(b.SourceParticle, a.SourceParticle);
      EXPECT_EQ(b.Weight, a.Weight);
      EXPECT_EQ(b.Energy, a.Energy);
    }
    EXPECT_EQ(first.Rand(), second.Rand());
  }
}
TEST(VmfSource, ZeroOpeningMatchesPencilWithoutRandomDraws) {
  for (double kappa : {0., 100., std::numeric_limits<double>::max()}) {
    SourceGeometry pencil({}, {1, 2, 3});
    SourceGeometry vmf({}, {1, 2, 3}, Model::Vmf, 0, kappa);
    Particle a(Photon, 0), b = a;
    Randomizer first(12), second(12);
    pencil.Initialize(a, first);
    vmf.Initialize(b, second);
    for (size_t i = 0; i < 3; ++i) EXPECT_EQ(a.Pdir[i], b.Pdir[i]);
    EXPECT_EQ(first.Rand(), second.Rand());
  }
}
TEST(VmfSource, ConstructorRejectsInvalidAndInapplicableParameters) {
  for (double kappa : {-1., std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()})
    Fatal([&] { SourceGeometry source({}, {0,0,1}, Model::Vmf, 10, kappa); });
  for (double angle : {-1., 181., std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()})
    Fatal([&] { SourceGeometry source({}, {0,0,1}, Model::Vmf, angle, 10); });
  Fatal([] { SourceGeometry source({}, {0,0,1}, Model::Cone, 10, 1); });
  Fatal([] { SourceGeometry source({}, {0,0,1}, Model::Pencil, 0, 1); });
  Fatal([] { SourceGeometry source({}, {0,0,0}, Model::Vmf, 10, 1); });
  const double tiny = std::numeric_limits<double>::denorm_min();
  if (!(Width(tiny) > 0))
    Fatal([&] { SourceGeometry source({}, {0,0,1}, Model::Vmf, tiny, 1); });
}
TEST(VmfSource, QuantilesMatchIndependentHundredDigitReference) {
  // Reference arithmetic is independent of the production expm1/log1p path.
  // 512 ulps relative allows trig/frame rounding in recovered chord offsets;
  // all selected offsets are normal in both double and extended long double.
  const std::array<std::pair<double,double>, 9> cases{{
      {0,180}, {1e-20,180}, {1e-12,10}, {1,10}, {100,10},
      {1e6,180}, {1e6,.001}, {1e300,180}, {100,1e-9}}};
  for (auto [kappa, angle] : cases) {
    SCOPED_TRACE(::testing::Message() << kappa << ", " << angle);
    SourceGeometry source({}, {0,0,1}, Model::Vmf, angle, kappa);
    Randomizer random(8307), reference(8307);
    Particle particle(Photon, 0);
    const High half = High(angle) * acos(High(-1)) / 360;
    const High width = 2 * sin(half) * sin(half);
    const High concentration(kappa), product = concentration * width;
    // exp(-1000)<1e-434, below the 100-digit reference precision.
    const High mass = product > 1000 ? High(1) : High(1 - exp(-product));
    for (int draw = 0; draw < 8; ++draw) {
      const High uniform(reference.Rand());
      const High azimuth = 2 * acos(High(-1)) * High(reference.Rand());
      const High expected = kappa == 0 ? High(uniform * width) :
          High(-log(High(1) - uniform * mass) / concentration);
      source.Initialize(particle, random);
      const auto value = expected.convert_to<coord_type>();
      const auto tolerance = 512 * std::numeric_limits<coord_type>::epsilon();
      EXPECT_LE(std::abs(Offset(particle, {0,0,1}) - value), tolerance * value);
      const High sine = sqrt(expected * (2 - expected));
      // The retained legacy frame for the z axis has transverse +y and -x.
      const std::array<High,3> direction{
          High(-sine*sin(azimuth)), High(sine*cos(azimuth)), High(1-expected)};
      for (size_t i = 0; i < 3; ++i)
        EXPECT_LE(std::abs(particle.Pdir[i] - direction[i].convert_to<coord_type>()),
                  tolerance * (i == 2 ? coord_type(1) : sine.convert_to<coord_type>()));
    }
    EXPECT_EQ(random.Rand(), reference.Rand());
  }
}
TEST(VmfSource, NonuniformSamplingPreservesNuclearPrimaryMetadata) {
  SourceGeometry source({1,2,3}, {1,-2,3}, Model::Vmf, 30, 100);
  Particle primary(56u, .01);
  primary.Energy = 1e20 * units.eV;
  primary.Weight = 7.25;
  primary.fPrevId = 81;
  primary.Ninteractions = 9;
  primary.interactionData[0] = &primary;
  const auto original = primary;
  Randomizer random(78);
  source.Initialize(primary, random);
  EXPECT_EQ(primary.id, original.id);
  EXPECT_EQ(primary.fPrevId, original.fPrevId);
  EXPECT_EQ(primary.SourceParticle, original.SourceParticle);
  EXPECT_EQ(primary.interactionData[0], original.interactionData[0]);
  EXPECT_EQ(primary.Ninteractions, original.Ninteractions);
  EXPECT_EQ(primary.Weight, original.Weight);
  EXPECT_EQ(primary.Energy, original.Energy);
  EXPECT_EQ(primary.MassNumber(), 56u);
  EXPECT_EQ(primary.AtomicNumber(), 26u);
  EXPECT_EQ(primary.Mass(), original.Mass());
  EXPECT_EQ(primary.Time.z(), original.Time.z());
  EXPECT_EQ(primary.fProductionTime.z(), original.fProductionTime.z());
}
TEST(VmfSource, ExtremeFiniteConcentrationsHaveFiniteUnitDirections) {
  for (double kappa : {1e-300, std::numeric_limits<double>::max()}) {
    for (double angle : {.001, 180.}) {
      SourceGeometry source({}, {0,0,1}, Model::Vmf, angle, kappa);
      Randomizer random(452), reference(452);
      Particle particle(Photon, 0);
      for (int draw = 0; draw < 20; ++draw) {
        source.Initialize(particle, random);
        coord_type norm = 0;
        for (auto component : particle.Pdir) {
          EXPECT_TRUE(std::isfinite(component));
          norm += component * component;
        }
        EXPECT_LE(std::abs(norm - 1), 64 * std::numeric_limits<coord_type>::epsilon());
        reference.Rand(); reference.Rand();
      }
      EXPECT_EQ(random.Rand(), reference.Rand());
    }
  }
}

using Parameters = std::tuple<double,double,int>; // kappa, degrees, axis index
class VmfDistribution : public ::testing::TestWithParam<Parameters> {};
TEST_P(VmfDistribution, FixedCdfEnsembleIsUniformAndWeightsArePreserved) {
  const auto [kappa, angle, axisIndex] = GetParam();
  const SpatialVector axis = SourceGeometry::NormalizeAxis(
      axisIndex == 0 ? SpatialVector{0,0,1} : SpatialVector{1,-2,3});
  SourceGeometry source({1,2,3}, axis, Model::Vmf, angle, kappa);
  Randomizer random(401715);
  Particle particle(Photon, 0);
  particle.Weight = 3;
  particle.Energy = 1e12 * units.eV;
  const auto id = particle.id;
  constexpr int count = 40000, binCount = 20;
  std::array<int,binCount> bins{};
  coord_type sum = 0;
  const coord_type width = Width(angle), eps = std::numeric_limits<coord_type>::epsilon();
  const coord_type mass = kappa == 0 ? 0 : -std::expm1(-coord_type(kappa)*width);
  for (int i = 0; i < count; ++i) {
    source.Initialize(particle, random);
    const auto offset = Offset(particle, axis);
    coord_type norm = 0;
    for (auto component : particle.Pdir) norm += component * component;
    EXPECT_LE(std::abs(norm - 1), 64 * eps);
    EXPECT_GE(offset, 0);
    // Chord recovery has absolute error O(eps*sqrt(U)) for a rotated axis.
    EXPECT_LE(offset, width + 64 * eps * (std::sqrt(width) + eps));
    const coord_type cdf = kappa == 0 ? offset / width :
        -std::expm1(-coord_type(kappa)*offset) / mass;
    ASSERT_GE(cdf, 0);
    ASSERT_LE(cdf, 1 + 64 * eps / std::sqrt(width));
    const int bin = std::min(binCount-1, int(cdf*binCount));
    ++bins[bin]; sum += cdf;
  }
  // Frozen checks: 40000 draws, seed 401715, 20 bins, six analytic sigmas.
  const double probability = 1. / binCount;
  for (int occupancy : bins)
    EXPECT_NEAR(occupancy, count * probability,
                6 * std::sqrt(count * probability * (1-probability)));
  EXPECT_NEAR(sum/count, .5, 6/std::sqrt(coord_type(12)*count));
  EXPECT_EQ(particle.Weight, 3);
  EXPECT_EQ(particle.Energy, 1e12 * units.eV);
  EXPECT_EQ(particle.id, id);
  for (size_t i = 0; i < 3; ++i) EXPECT_EQ(particle.X[i], coord_type(i+1));
}
std::vector<Parameters> Cases() {
  std::vector<Parameters> cases;
  for (double kappa : {0.,1.,100.,1000000.})
    for (double angle : {.001,10.,180.})
      for (int axis : {0,1}) cases.emplace_back(kappa,angle,axis);
  return cases;
}
INSTANTIATE_TEST_SUITE_P(FrozenGrid, VmfDistribution, ::testing::ValuesIn(Cases()));
} // namespace
