#include "ParticleStack.h"
#include "PropagationEngine.h"
#include "SpatialGeometry.h"
#include "Logger.h"
#include <cmath>
#include <functional>
#include <gtest/gtest.h>
#include <limits>
#include <memory>
#include <sstream>

namespace {
using namespace mcray;
class SpatialGeometryTest : public ::testing::Test {
  void SetUp() override {
    if (!cosmology.IsInitialized())
      cosmology.Init();
  }
};
void Fatal(const std::function<void()> &action) {
  try {
    action();
    FAIL() << "Expected invalid geometry";
  } catch (Utils::Exception *error) {
    std::unique_ptr<Utils::Exception> owner(error);
    EXPECT_FALSE(error->Retryable());
  }
}
class Capture : public Utils::TSmartReferencedObj<IOutput> {
public:
  std::vector<Particle> particles;
  void AddParticle(const Particle &p) override { particles.push_back(p); }
  void Flush(bool) override {}
};
class ProgressWatchdog : public ILogger {
public:
  unsigned calls=0;
  void log(const Particle&) override {
    // Make the regression fail promptly on old code instead of hanging CI.
    if (++calls>12)
      Utils::Exception::ThrowFatal("Test watchdog: variable-rate propagation stalled");
  }
};

TEST_F(SpatialGeometryTest, VariableRateRejectsHalfStepsThatCannotAdvanceTime) {
  Particle p(Photon,.01);
  p.Energy=1e12*units.eV;
  cosmo_time initial=p.Time.t();
  auto next=std::nextafter(initial,std::numeric_limits<cosmo_time>::infinity());
  const cosmo_time half=(next-initial)/2;
  // Adjacent normal binary values alternate their tie-to-even rounding.
  if (initial+half!=initial) {
    initial=std::nextafter(initial,-std::numeric_limits<cosmo_time>::infinity());
    next=std::nextafter(initial,std::numeric_limits<cosmo_time>::infinity());
  }
  ASSERT_EQ(initial+(next-initial)/2,initial);
  p.Time.setT(initial);
  CosmoTime end; end.setT(next);
  Result result; result.SetEndTime(end);
  Utils::SmartPtr<Capture> capture=new Capture;
  result.AddOutput(capture);
  ParticleStack stack; stack.AddPrimary(p);
  PropagationEngine engine(stack,result,17);
  ProgressWatchdog watchdog; engine.SetLogger(&watchdog);
  std::ostringstream errors;
  auto* previous=std::cerr.rdbuf(errors.rdbuf());
  engine.Run(.1);
  std::cerr.rdbuf(previous);
  result.Flush();
  EXPECT_EQ(engine.FailureCount(),1u);
  EXPECT_TRUE(capture->particles.empty());
  EXPECT_LT(watchdog.calls,12u);
  EXPECT_NE(errors.str().find("below available time resolution"),std::string::npos);
  EXPECT_EQ(errors.str().find("Test watchdog"),std::string::npos);
}
TEST_F(SpatialGeometryTest,
       PencilInitializesPositionAndPreservesRandomStreamAndIdentity) {
  SourceGeometry source({units.Mpc, -2 * units.Mpc, 3 * units.Mpc}, {0, 0, 5});
  Particle p(56u, 0.01);
  p.Energy = 1e20 * units.eV;
  p.Weight = 2.5;
  const auto birth = p.Time.t();
  const auto id = p.id;
  Randomizer random(17), control(17);
  source.Initialize(p, random);
  EXPECT_EQ(p.X[0], units.Mpc);
  EXPECT_EQ(p.X[1], -2 * units.Mpc);
  EXPECT_EQ(p.X[2], 3 * units.Mpc);
  EXPECT_EQ(p.Pdir[0], 0);
  EXPECT_EQ(p.Pdir[1], 0);
  EXPECT_EQ(p.Pdir[2], 1);
  EXPECT_EQ(p.MassNumber(), 56u);
  EXPECT_EQ(p.AtomicNumber(), 26u);
  EXPECT_EQ(p.id, id);
  EXPECT_EQ(p.Time.t(), birth);
  EXPECT_DOUBLE_EQ(p.Weight, 2.5);
  EXPECT_DOUBLE_EQ(p.Energy, 1e20 * units.eV);
  EXPECT_EQ(random.Rand(), control.Rand());
  const auto huge = SourceGeometry::NormalizeAxis({1e300L, 1e300L, -1e300L});
  const auto tiny = SourceGeometry::NormalizeAxis({1e-300L, 1e-300L, -1e-300L});
  for (size_t i = 0; i < 3; ++i)
    EXPECT_NEAR(huge[i], tiny[i], 1e-15);
  Fatal([] { SourceGeometry source({}, {0, 0, 0}); });
  Fatal([] {
    SourceGeometry source({},
                          {0, 0, std::numeric_limits<coord_type>::infinity()});
  });
  Fatal([] {
    SourceGeometry source({}, {0, 0, 1}, SourceGeometry::AngularModel::Pencil,
                          1);
  });
}
TEST_F(SpatialGeometryTest, ConeIsUniformInSolidAngleForRotatedAxis) {
  const SpatialVector axis = SourceGeometry::NormalizeAxis({1, -2, 3});
  for (double angle : {20., 180.}) {
    SourceGeometry source({}, axis, SourceGeometry::AngularModel::Cone, angle);
    Randomizer random(8721);
    Particle p(Photon, 0.01);
    p.Energy = 1e12 * units.eV;
    p.Weight = 3;
    const coord_type lower = std::cos(angle * std::acos(coord_type(-1)) / 180);
    coord_type sum = 0, sum2 = 0, perpendicular = 0;
    const SpatialVector transverse = SourceGeometry::NormalizeAxis({2, 1, 0});
    constexpr int count = 40000;
    for (int i = 0; i < count; ++i) {
      source.Initialize(p, random);
      coord_type dot = 0, cross = 0, norm = 0;
      for (size_t a = 0; a < 3; ++a) {
        dot += p.Pdir[a] * axis[a];
        cross += p.Pdir[a] * transverse[a];
        norm += p.Pdir[a] * p.Pdir[a];
      }
      EXPECT_GE(dot, lower - 2e-15);
      EXPECT_LE(dot, 1 + 2e-15);
      EXPECT_NEAR(norm, 1, 2e-15);
      sum += dot;
      sum2 += dot * dot;
      perpendicular += cross;
    }
    // Fixed seed; analytic uniform-cosine moments, with six standard errors.
    const coord_type width = 1 - lower, mean = (1 + lower) / 2;
    EXPECT_NEAR(sum / count, mean,
                6 * width / std::sqrt(coord_type(12) * count));
    EXPECT_NEAR(sum2 / count, (1 + lower + lower * lower) / 3,
                6 / std::sqrt(coord_type(count)));
    EXPECT_NEAR(perpendicular / count, 0, 6 / std::sqrt(coord_type(count)));
    EXPECT_DOUBLE_EQ(p.Weight, 3);
  }
}
TEST_F(SpatialGeometryTest, MembershipIsStrictAndIncludesBirth) {
  SphereObserver inside({1, 2, 3}, 2),
      outside({1, 2, 3}, 2, SphereObserver::StopWhen::Outside);
  Particle p(Photon, 0.01);
  p.X[0] = 1;
  p.X[1] = 2;
  p.X[2] = 3;
  EXPECT_TRUE(inside.Reached(p));
  EXPECT_FALSE(outside.Reached(p));
  p.X[0] = 3;
  EXPECT_FALSE(inside.Reached(p));
  EXPECT_FALSE(outside.Reached(p));
  p.X[0] = 3.1;
  EXPECT_FALSE(inside.Reached(p));
  EXPECT_TRUE(outside.Reached(p));
  Fatal([] { SphereObserver observer({}, 0); });
}
TEST_F(SpatialGeometryTest,
       StepNeverExceedsRemainingTimeAndResolutionFailuresAreFatal) {
  const auto radius = 0.01L * units.Mpc;
  Particle p(Photon, 0.01);
  for (auto mode :
       {SphereObserver::StopWhen::Inside, SphereObserver::StopWhen::Outside}) {
    SphereObserver observer({}, radius, mode);
    p.X[0] = radius;
    const cosmo_time remaining = 1e-8L * radius;
    EXPECT_EQ(observer.LimitStep(p, remaining), remaining);
    EXPECT_EQ(observer.LimitStep(p, 0), 0);
  }
  p.X[0] = 0;
  SphereObserver tooSmall({}, std::numeric_limits<double>::min());
  Fatal([&] { tooSmall.LimitStep(p, units.Mpc); });
}
TEST_F(SpatialGeometryTest,
       EngineStopsAtMembershipInBothRatePathsForPhotonAndNucleus) {
  for (double variable : {0., 0.1})
    for (bool nuclear : {false, true})
      for (auto mode : {SphereObserver::StopWhen::Inside,
                        SphereObserver::StopWhen::Outside}) {
        const coord_type radius = 0.01L * units.Mpc;
        const SpatialVector center = mode == SphereObserver::StopWhen::Inside
                                         ? SpatialVector{0, 0, 3 * radius}
                                         : SpatialVector{};
        auto observer = std::make_shared<SphereObserver>(center, radius, mode);
        Result result;
        result.SetObserver(observer);
        Particle p = nuclear ? Particle(56u, 0.01) : Particle(Photon, 0.01);
        p.Energy = 1e20 * units.eV;
        p.Weight = 2.5;
        CosmoTime end = p.Time;
        end += 5 * radius;
        result.SetEndTime(end);
        Utils::SmartPtr<Capture> capture = new Capture;
        result.AddOutput(capture);
        ParticleStack stack;
        stack.AddPrimary(p);
        PropagationEngine engine(stack, result, 17);
        engine.Run(variable);
        result.Flush();
        ASSERT_EQ(capture->particles.size(), 1u);
        const auto &final = capture->particles[0];
        EXPECT_TRUE(observer->Reached(final));
        EXPECT_LE(final.Time.t(), end.t());
        const coord_type expected =
            mode == SphereObserver::StopWhen::Inside ? 2 * radius : radius;
        const coord_type bound =
            (mode == SphereObserver::StopWhen::Inside ? 1.1e-3L : 1.1e-7L) *
            radius;
        EXPECT_NEAR(final.X[2], expected, bound);
        EXPECT_DOUBLE_EQ(final.Weight, p.Weight);
        EXPECT_DOUBLE_EQ(final.Energy, p.Energy);
        if (nuclear) {
          EXPECT_EQ(final.MassNumber(), 56u);
          EXPECT_EQ(final.AtomicNumber(), 26u);
        }
      }
}
TEST_F(SpatialGeometryTest,
       EngineStopsImmediatelyForInitialHitAndMissRespectsEndTime) {
  for (double variable : {0., 0.1}) {
    const auto radius = 0.01L * units.Mpc;
    for (bool hit : {true, false}) {
      auto observer = std::make_shared<SphereObserver>(
          hit ? SpatialVector{} : SpatialVector{3 * radius, 0, 0}, radius);
      Result result;
      result.SetObserver(observer);
      Particle p(Photon, 0.01);
      p.Energy = 1e12 * units.eV;
      CosmoTime end = p.Time;
      end += 5 * radius;
      result.SetEndTime(end);
      Utils::SmartPtr<Capture> capture = new Capture;
      result.AddOutput(capture);
      ParticleStack stack;
      stack.AddPrimary(p);
      PropagationEngine engine(stack, result, 17);
      engine.Run(variable);
      result.Flush();
      ASSERT_EQ(capture->particles.size(), 1u);
      EXPECT_EQ(observer->Reached(capture->particles[0]), hit);
      EXPECT_EQ(capture->particles[0].Time.t(), hit ? p.Time.t() : end.t());
    }
  }
}
} // namespace
