// SPDX-License-Identifier: MIT
#include "SpatialGeometry.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace mcray {
namespace {
[[noreturn]] void GeometryError(const std::string &why) {
  Utils::Exception::ThrowFatal("Spatial geometry: " + why);
  std::terminate();
}
void ValidateVector(const SpatialVector &vector) {
  for (coord_type v : vector)
    if (!std::isfinite(v))
      GeometryError("coordinates/direction must be finite");
}
coord_type VmfLog1p(coord_type x) {
  if (std::abs(x) > coord_type(0.125L))
    return std::log1p(x);
  if (1 + x == 1)
    return x; // The quadratic correction is below coordinate precision.
  // log(1+x)=2*atanh(x/(2+x)). The series ratio is <=1/225 here,
  // so a few terms give coordinate precision without subtracting from one.
  // This avoids small-argument fyl2xp1 accuracy loss under x87 emulation.
  const coord_type y = x / (2 + x), square = y * y;
  coord_type power = y, sum = y;
  for (unsigned denominator = 3; ; denominator += 2) {
    power *= square;
    const coord_type next = sum + power / denominator;
    if (next == sum)
      return 2 * sum;
    sum = next;
  }
}
SpatialVector Cross(const SpatialVector &a, const SpatialVector &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}
} // namespace

SpatialVector SourceGeometry::NormalizeAxis(SpatialVector axis) {
  ValidateVector(axis);
  const coord_type scale =
      std::max({std::abs(axis[0]), std::abs(axis[1]), std::abs(axis[2])});
  if (!(scale > 0))
    GeometryError("source axis must not be zero");
  for (auto &v : axis)
    v /= scale;
  const coord_type norm = std::hypot(axis[0], axis[1], axis[2]);
  for (auto &v : axis)
    v /= norm;
  return axis;
}

SourceGeometry::SourceGeometry(SpatialVector position, SpatialVector axis,
                               AngularModel model, double halfAngleDegrees,
                               double kappa)
    : position_(position), axis_(NormalizeAxis(axis)), model_(model), kappa_(kappa) {
  ValidateVector(position_);
  if ((model != AngularModel::Pencil && model != AngularModel::Cone &&
       model != AngularModel::Vmf) ||
      !std::isfinite(halfAngleDegrees) || halfAngleDegrees < 0 ||
      halfAngleDegrees > 180 ||
      (model == AngularModel::Pencil && halfAngleDegrees != 0))
    GeometryError("angular model/half-angle is invalid");
  if (!std::isfinite(kappa) || kappa < 0 ||
      (model != AngularModel::Vmf && kappa != 0))
    GeometryError("kappa must be finite, nonnegative and used only with vmf");
  const coord_type half = halfAngleDegrees * std::acos(coord_type(-1)) / 360;
  // 1-cos(theta) evaluated without cancellation for a narrow cone.
  coneWidth_ = 2 * std::sin(half) * std::sin(half);
  if (model == AngularModel::Vmf) {
    if (halfAngleDegrees > 0 && !(coneWidth_ > 0))
      GeometryError("positive vmf cone width is not representable");
    if (kappa_ > 0 && coneWidth_ > 0) {
      // An overflowing positive product has exp(-kappa*U)=0 at this
      // precision. Avoid forming it, including on double-width long double.
      if (coneWidth_ > 1 && kappa_ >
          std::numeric_limits<coord_type>::max() / coneWidth_) {
        vmfMass_ = 1;
      } else {
        const coord_type product = kappa_ * coneWidth_;
        // The relative uniform-limit correction is at most product/2.
        // Below epsilon it is unresolvable; this also handles underflow.
        if (product > std::numeric_limits<coord_type>::epsilon())
          vmfMass_ = -std::expm1(-product);
      }
    }
  }
  SpatialVector reference{};
  size_t smallest = 0;
  for (size_t i = 1; i < 3; ++i)
    if (std::abs(axis_[i]) < std::abs(axis_[smallest]))
      smallest = i;
  reference[smallest] = 1;
  transverse_ = NormalizeAxis(Cross(axis_, reference));
  binormal_ = Cross(axis_, transverse_);
}

coord_type SourceGeometry::AngularOffset(coord_type uniform) const {
  if (model_ != AngularModel::Vmf || vmfMass_ == 0)
    return uniform * coneWidth_; // Includes the exact legacy cone/kappa=0 path.
  const coord_type offset = -VmfLog1p(-uniform * vmfMass_) / kappa_;
  // The inverse CDF is supported on [0,U]; bound possible edge rounding.
  return std::min(offset, coneWidth_);
}

void SourceGeometry::Initialize(Particle &particle, Randomizer &random) const {
  SpatialVector direction = axis_;
  // The legacy pencil/default and a zero-width cone consume no random values.
  if (model_ != AngularModel::Pencil && coneWidth_ > 0) {
    const coord_type offset = AngularOffset(random.Rand());
    const coord_type cosine = 1 - offset;
    const coord_type sine =
        std::sqrt(std::max(coord_type(0), offset * (2 - offset)));
    const coord_type azimuth = 2 * std::acos(coord_type(-1)) * random.Rand();
    for (size_t i = 0; i < 3; ++i)
      direction[i] =
          cosine * axis_[i] + sine * (std::cos(azimuth) * transverse_[i] +
                                      std::sin(azimuth) * binormal_[i]);
    direction = NormalizeAxis(direction);
  }
  for (size_t i = 0; i < 3; ++i) {
    particle.X[i] = position_[i];
    particle.Pdir[i] = direction[i];
  }
}

SphereObserver::SphereObserver(SpatialVector center, coord_type radius,
                               StopWhen stop)
    : center_(center), radius_(radius), stop_(stop) {
  ValidateVector(center_);
  if (!std::isfinite(radius_) || !(radius_ > 0) ||
      (stop_ != StopWhen::Inside && stop_ != StopWhen::Outside))
    GeometryError("sphere requires a finite positive radius and inside/outside "
                  "stop rule");
  // Preserve the fork's endpoint sampling scales, with a cap at remaining time.
  minimumStep_ = radius_ * (stop_ == StopWhen::Inside ? coord_type(1e-3L)
                                                      : coord_type(1e-7L));
  if (!(minimumStep_ > 0) || !std::isfinite(minimumStep_))
    GeometryError("sphere step is not representable");
}

SpatialVector SphereObserver::Offset(const Particle &particle) const {
  SpatialVector offset;
  for (size_t i = 0; i < 3; ++i)
    offset[i] = particle.X[i] - center_[i];
  ValidateVector(offset);
  return offset;
}

bool SphereObserver::Reached(const Particle &particle) const {
  const auto offset = Offset(particle);
  // Hypot avoids squaring overflow; equality remains outside both strict tests.
  const coord_type distance = std::hypot(offset[0], offset[1], offset[2]);
  return stop_ == StopWhen::Inside ? distance < radius_ : distance > radius_;
}

cosmo_time SphereObserver::LimitStep(const Particle &particle,
                                     cosmo_time remaining) const {
  if (!std::isfinite(remaining))
    GeometryError("remaining time must be finite");
  if (remaining <= 0)
    return remaining;
  const auto offset = Offset(particle);
  coord_type distance;
  if (stop_ == StopWhen::Outside)
    distance = radius_ - std::hypot(offset[0], offset[1], offset[2]);
  else
    distance = std::max({std::abs(offset[0]), std::abs(offset[1]),
                         std::abs(offset[2])}) -
               radius_;
  const cosmo_time step = std::min(remaining, std::max(minimumStep_, distance));
  if (!(step > 0) || !std::isfinite(step) ||
      particle.Time.t() + step <= particle.Time.t())
    GeometryError("observer step is below available time resolution");
  return step;
}
} // namespace mcray
