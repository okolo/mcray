// SPDX-License-Identifier: MIT
#ifndef MCRAY_SPATIAL_GEOMETRY_H
#define MCRAY_SPATIAL_GEOMETRY_H

#include "Particle.h"
#include "Randomizer.h"
#include <array>

namespace mcray {

using SpatialVector = std::array<coord_type, 3>;

// Pencil/cone preserve weights. Vmf uses the original uniform-cone proposal
// and multiplies the primary weight by the full-sphere vMF density.
class SourceGeometry {
public:
  enum class AngularModel : int { Pencil, Cone, Vmf };
  SourceGeometry(SpatialVector position, SpatialVector axis,
                 AngularModel model = AngularModel::Pencil,
                 double halfAngleDegrees = 0, double kappa = 0);
  static SpatialVector NormalizeAxis(SpatialVector axis);
  void Initialize(Particle &particle, Randomizer &random) const;

private:
  SpatialVector position_, axis_, transverse_, binormal_;
  AngularModel model_;
  coord_type coneWidth_;
  double kappa_;
  bool vmfActive_;
};

// One immutable observer is shared by stopping, step selection and output.
// Future crossing-based observers can implement this interface independently.
class SpatialObserver {
public:
  virtual ~SpatialObserver() = default;
  virtual bool Reached(const Particle &particle) const = 0;
  virtual cosmo_time LimitStep(const Particle &particle,
                               cosmo_time remaining) const = 0;
};

class SphereObserver final : public SpatialObserver {
public:
  enum class StopWhen : int { Inside, Outside };
  // Coordinates/radius use internal comoving length units.
  SphereObserver(SpatialVector center, coord_type radius,
                 StopWhen stop = StopWhen::Inside);
  bool Reached(const Particle &particle) const override;
  cosmo_time LimitStep(const Particle &particle,
                       cosmo_time remaining) const override;

private:
  SpatialVector Offset(const Particle &particle) const;
  SpatialVector center_;
  coord_type radius_, minimumStep_;
  StopWhen stop_;
};

} // namespace mcray
#endif
