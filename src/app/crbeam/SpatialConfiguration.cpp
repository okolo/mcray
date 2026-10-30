#include "SpatialConfiguration.h"
#include <cmath>
#include <stdexcept>

namespace {
mcray::SpatialVector Vector(const YAML::Node &node) {
  mcray::SpatialVector result;
  for (size_t i = 0; i < 3; ++i)
    result[i] = node[i].as<double>();
  return result;
}
mcray::SpatialVector Internal(mcray::SpatialVector vector) {
  for (auto &v : vector) {
    v *= mcray::units.Mpc;
    if (!std::isfinite(v))
      throw std::runtime_error(
          "Spatial coordinate is not representable in internal units");
  }
  return vector;
}
} // namespace
SpatialConfiguration::SpatialConfiguration(const YAML::Node &effective) {
  const auto source = effective["source"];
  if (source["position_Mpc"])
    position_ = Vector(source["position_Mpc"]);
  if (source["angular_distribution"]) {
    const auto angular = source["angular_distribution"];
    axis_ = Vector(angular["axis"]);
    if (angular["model"].as<std::string>() == "cone")
      angular_ = mcray::SourceGeometry::AngularModel::Cone;
    if (angular["model"].as<std::string>() == "vmf") {
      angular_ = mcray::SourceGeometry::AngularModel::Vmf;
      kappa_ = angular["kappa"].as<double>();
    }
    halfAngle_ = angular["half_angle_deg"].as<double>();
  }
  const auto propagation = effective["propagation"];
  if (propagation && propagation["end_redshift"])
    endRedshift_ = propagation["end_redshift"].as<double>();
  const auto observer = effective["observer"];
  if (observer && observer["model"].as<std::string>() == "sphere") {
    sphere_ = true;
    center_ = Vector(observer["center_Mpc"]);
    radius_ = observer["radius_Mpc"].as<double>();
    if (observer["stop_when"].as<std::string>() == "outside")
      stop_ = mcray::SphereObserver::StopWhen::Outside;
  }
  const auto output = effective["output"];
  if (output && output["primary_direction"])
    primaryDirection_ = output["primary_direction"].as<bool>();
  // Preflight unit conversions without creating output or sampling directions.
  Source();
  Observer();
}
mcray::SourceGeometry SpatialConfiguration::Source() const {
  return mcray::SourceGeometry(Internal(position_), axis_, angular_,
                               halfAngle_, kappa_);
}
std::shared_ptr<const mcray::SpatialObserver>
SpatialConfiguration::Observer() const {
  if (!sphere_)
    return {};
  const auto radius =
      static_cast<mcray::coord_type>(radius_) * mcray::units.Mpc;
  return std::make_shared<const mcray::SphereObserver>(Internal(center_),
                                                       radius, stop_);
}
bool SpatialConfiguration::HasSpatialSettings() const {
  return sphere_ || primaryDirection_ || endRedshift_ != 0 || halfAngle_ != 0 ||
         angular_ != mcray::SourceGeometry::AngularModel::Pencil ||
         position_ != mcray::SpatialVector{} ||
         axis_ != mcray::SpatialVector{0, 0, 1};
}
