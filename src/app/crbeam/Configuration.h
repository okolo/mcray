#ifndef MCRAY_CRBEAM_CONFIGURATION_H
#define MCRAY_CRBEAM_CONFIGURATION_H

#include "CmdLine.h"
#include <cstdint>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

// Both input forms reach the existing process setup and validation path.
class Configuration {
public:
  Configuration(int argc, char **argv, cors::cmdline::CmdInfo *commands);
  int Argc() const { return static_cast<int>(fArguments.size()); }
  char **Argv() { return fArgv.data(); }
  bool IsYaml() const { return !fInputFile.empty(); }
  bool HasSeed() const { return fHasSeed; }
  unsigned long Seed() const { return fSeed; }
  YAML::Node &Effective() { return fEffective; }
  const YAML::Node &Effective() const { return fEffective; }
  void ResolveExecution(unsigned long seed, int threads, int batch,
                        const std::string &output, unsigned a, unsigned z);
  void SetAdditionalTable(const std::string &path);
  static std::string NuclearTables();
  static void VerifyBundle(const std::string &directory,
                           bool requireHashes = false);
  static std::string HashFile(const std::string &path);
  void WriteConfig(const std::string &directory) const;
  YAML::Node Provenance(const std::string &nuclearDirectory) const;

private:
  std::vector<std::string> fArguments, fOriginal;
  std::vector<char *> fArgv;
  YAML::Node fEffective;
  std::string fInputFile, fInputHash;
  bool fHasSeed = false, fAutomaticSeed = true;
  unsigned long fSeed = 0;
};

// A run is failed unless explicitly finalized, including pointer exceptions.
class RunRecord {
public:
  RunRecord(const std::string &directory, YAML::Node metadata);
  ~RunRecord();
  void Update(YAML::Node metadata);
  void Finish(bool interrupted, unsigned long failures = 0);

private:
  void Write();
  std::string fDirectory;
  YAML::Node fMetadata;
  bool fFinished = false;
};
#endif
