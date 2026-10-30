#ifndef MCRAY_SPATIAL_CONFIGURATION_H
#define MCRAY_SPATIAL_CONFIGURATION_H

#include "SpatialGeometry.h"
#include <memory>
#include <yaml-cpp/yaml.h>

class SpatialConfiguration {
public:
  explicit SpatialConfiguration(const YAML::Node &effective);
  mcray::SourceGeometry Source() const;
  std::shared_ptr<const mcray::SpatialObserver> Observer() const;
  bool HasSpatialSettings() const;
  double EndRedshift() const { return endRedshift_; }
  bool PrimaryDirection() const { return primaryDirection_; }

private:
  mcray::SpatialVector position_{}, axis_{0, 0, 1}, center_{};
  mcray::SourceGeometry::AngularModel angular_ =
      mcray::SourceGeometry::AngularModel::Pencil;
  mcray::SphereObserver::StopWhen stop_ =
      mcray::SphereObserver::StopWhen::Inside;
  bool sphere_ = false, primaryDirection_ = false;
  double halfAngle_ = 0, kappa_ = 0, radius_ = 0, endRedshift_ = 0;
};
#endif
