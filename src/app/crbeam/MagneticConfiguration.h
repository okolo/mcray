#ifndef MCRAY_MAGNETIC_CONFIGURATION_H
#define MCRAY_MAGNETIC_CONFIGURATION_H

#include "CmdLine.h"
#include "Deflection3D.h"
#include <map>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

// Version 2 is normalized independently of the legacy scalar/flag bindings.
YAML::Node NormalizeMagneticField(const YAML::Node &requested,
                                  const YAML::Node &legacy,
                                  const std::string &base);
YAML::Node MagneticFieldFromCli(cors::cmdline::CmdLine &cli,
                                const YAML::Node &legacy);
void MagneticLegacyArguments(const YAML::Node &field,
                             std::vector<std::string> &arguments);

// Load immutable grids before output creation; retain their hashes from load
// time.
class ConfiguredMagneticField {
public:
  explicit ConfiguredMagneticField(const YAML::Node &field);
  Utils::SmartPtr<Interactions::MagneticField>
  Build(mcray::Randomizer &random) const;
  const YAML::Node &Provenance() const { return provenance_; }

private:
  void LoadGrids(const YAML::Node &field);
  Utils::SmartPtr<Interactions::MagneticField>
  BuildNode(const YAML::Node &field, mcray::Randomizer &random) const;
  YAML::Node field_, provenance_;
  std::map<std::string, Utils::SmartPtr<Interactions::MagneticField>> grids_;
};
#endif
