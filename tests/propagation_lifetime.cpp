#include "PropagationEngine.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <gtest/gtest.h>
#include <memory>
#include <omp.h>
#include <thread>

namespace {
using namespace mcray;
struct LifetimeState {
  std::atomic<int> clones{0}, destroyed{0}, activeDestructors{0};
  std::atomic<int> tokensDestroyed{0};
  std::atomic<bool> overlappingDestructors{false};
};

class SharedToken : public Utils::SmartReferencedObj {
public:
  explicit SharedToken(std::shared_ptr<LifetimeState> state)
      : state_(std::move(state)) {}
  ~SharedToken() override { ++state_->tokensDestroyed; }

private:
  std::shared_ptr<LifetimeState> state_;
};

// Each engine clone owns one reference to the same non-atomic intrusive token,
// just as Deflection3D clones share a field and engines share output filters.
class LifetimeProbe : public CELInteraction {
public:
  LifetimeProbe(std::shared_ptr<LifetimeState> state,
                const Utils::SmartPtr<SharedToken> &token, bool clone = false)
      : state_(std::move(state)), token_(token), clone_(clone) {}
  ~LifetimeProbe() override {
    if (!clone_)
      return;
    if (state_->activeDestructors.fetch_add(1) != 0)
      state_->overlappingDestructors = true;
    // Widen the overlap window without touching any physics/RNG state. The old
    // parallel cleanup reliably overlaps here; serialized cleanup cannot do so.
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    --state_->activeDestructors;
    ++state_->destroyed;
  }
  double Rate(const Particle &) const override { return 0; }
  CELInteraction *Clone() const override {
    ++state_->clones;
    return new LifetimeProbe(state_, token_, true);
  }

private:
  std::shared_ptr<LifetimeState> state_;
  Utils::SmartPtr<SharedToken> token_;
  bool clone_;
};

TEST(PropagationLifetime, WorkerCleanupSerializesSharedReferenceRelease) {
  if (omp_get_thread_limit() < 2)
    GTEST_SKIP() << "At least two OpenMP workers are required";
  omp_set_dynamic(0);
  for (const int requested : {2, 4}) {
    const int workers = std::min(requested, omp_get_thread_limit());
    auto state = std::make_shared<LifetimeState>();
    Utils::SmartPtr<SharedToken> token = new SharedToken(state);
    {
      ParticleStack stack;
      Result result;
      PropagationEngine engine(stack, result, 724106);
      engine.AddInteraction(new LifetimeProbe(state, token));
      ASSERT_EQ(token->RefCount(), 2);
      // Empty stacks isolate lifetime management. Every invocation still
      // creates and destroys the full worker team, including all interaction
      // clones.
      for (int batch = 1; batch <= 3; ++batch) {
        engine.RunMultithread(0, workers);
        EXPECT_EQ(state->clones, workers * batch);
        EXPECT_EQ(state->destroyed, workers * batch);
        EXPECT_EQ(state->activeDestructors, 0);
        EXPECT_FALSE(state->overlappingDestructors);
        EXPECT_EQ(token->RefCount(), 2);
      }
    }
    EXPECT_EQ(token->RefCount(), 1);
    EXPECT_EQ(state->tokensDestroyed, 0);
    token = nullptr;
    EXPECT_EQ(state->tokensDestroyed, 1);
  }
}
} // namespace
