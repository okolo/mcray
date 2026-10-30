#ifndef MCRAY_SPATIAL_OUTPUT_H
#define MCRAY_SPATIAL_OUTPUT_H

#include "Output.h"
#include "SpatialGeometry.h"
#include <memory>

namespace mcray {
// The fork-compatible layout has one immutable coordinate scale per file.
class SpatialOutput : public RawOutput {
public:
  SpatialOutput(std::string directory, bool overwrite, bool saveSourceDirection,
                bool saveId, std::shared_ptr<const SpatialObserver> observer);
  void WriteHeader(std::ostream &out, const Particle &first) override;
  void Write(std::ostream &out, const Particle &particle) override;

private:
  bool saveSourceDirection_, saveId_;
  std::shared_ptr<const SpatialObserver> observer_;
  std::map<std::string, coord_type> scales_;
};
} // namespace mcray
#endif
