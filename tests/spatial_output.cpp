#include "SpatialOutput.h"
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <unistd.h>

namespace {
using namespace mcray;
namespace fs = std::filesystem;
class SpatialOutputTest : public ::testing::Test {
protected:
  fs::path root;
  void SetUp() override {
    if (!cosmology.IsInitialized())
      cosmology.Init();
    root = fs::temp_directory_path() /
           ("mcray-output-" + std::to_string(getpid()));
    fs::remove_all(root);
  }
  void TearDown() override { fs::remove_all(root); }
  static std::string Read(const fs::path &path) {
    std::ifstream file(path);
    return {std::istreambuf_iterator<char>(file), {}};
  }
  static coord_type Scale(const std::string &text) {
    return std::stold(text.substr(text.find("T/Mpc = ") + 8));
  }
  static std::vector<std::vector<coord_type>> Rows(const std::string &text) {
    std::vector<std::vector<coord_type>> result;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream fields(line);
      std::vector<coord_type> row;
      coord_type value;
      while (fields >> value)
        row.push_back(value);
      result.push_back(row);
    }
    return result;
  }
};
TEST_F(SpatialOutputTest, ScaleIsPerFileAndSurvivesLaterFlushAndIsotopeGroups) {
  Particle source(Photon, .01);
  source.Energy = 1e12 * units.eV;
  SpatialOutput output(root.string(), false, true, false, {});
  Particle first = source;
  first.SourceParticle = &source;
  first.Time += units.Mpc;
  first.X[0] = 2 * units.Mpc;
  output.AddParticle(first);
  output.Flush(false);
  Particle other = first;
  other.SetPsbNucleus(56);
  other.Time += units.Mpc;
  output.AddParticle(other);
  output.Flush(false);
  other.SetPsbNucleus(55);
  other.Time += units.Mpc;
  output.AddParticle(other);
  output.Flush(false);
  first.Time += 3 * units.Mpc;
  first.X[0] = 5 * units.Mpc;
  output.AddParticle(first);
  output.Flush(true);
  const auto photons = Read(root / "photon");
  const auto rows = Rows(photons);
  ASSERT_EQ(rows.size(), 2u);
  EXPECT_NEAR(Scale(photons), 1, 1e-12);
  EXPECT_NEAR(rows[0][2], 2, 1e-12);
  EXPECT_NEAR(rows[1][2], 5, 1e-12);
  EXPECT_NEAR(rows[1].back(), 4, 1e-12);
  EXPECT_NEAR(Scale(Read(root / Particle(56u, .01).Name())), 2, 1e-12);
  EXPECT_NEAR(Scale(Read(root / other.Name())), 3, 1e-12);
}
TEST_F(SpatialOutputTest, BirthHitUsesCoordinateFallbackAndAppendsRealFlight) {
  Particle source(Photon, .01);
  source.Energy = 1e12 * units.eV;
  source.X[0] = 3 * units.Mpc;
  source.X[2] = -2 * units.Mpc;
  auto observer = std::make_shared<SphereObserver>(
      SpatialVector{3 * units.Mpc, 0, -2 * units.Mpc}, units.Mpc);
  SpatialOutput output(root.string(), false, true, false, observer);
  Particle hit = source;
  hit.SourceParticle = &source;
  output.AddParticle(hit);
  output.Flush(false);
  hit.Time += units.Mpc;
  hit.X[0] += units.Mpc;
  output.AddParticle(hit);
  output.Flush(true);
  const auto text = Read(root / "photon");
  EXPECT_NE(text.substr(0, text.find('\n')).find("normalization = fallback"),
            std::string::npos);
  EXPECT_NE(text.find("\tin_target\tsource_x\tsource_y\tsource_z\tz_src\tE_src/"
                      "eV\tPathMpc\tFlightMpc\n"),
            std::string::npos);
  EXPECT_NEAR(Scale(text), 3, 1e-12);
  const auto rows = Rows(text);
  ASSERT_EQ(rows.size(), 2u);
  ASSERT_EQ(rows[0].size(), 17u);
  EXPECT_EQ(rows[0][9], 1);
  EXPECT_EQ(rows[0].back(), 0);
  EXPECT_NEAR(rows[1][2] * Scale(text), 4, 1e-12);
  EXPECT_NEAR(rows[1].back(), 1, 1e-12);
  for (const auto &row : rows)
    for (const auto value : row)
      EXPECT_TRUE(std::isfinite(value));
}
TEST_F(SpatialOutputTest, TinyPositiveTimeUsesLegacyChoiceAndDecimalHeader) {
  Particle source(Photon, 1e-8);
  source.Energy = units.eV;
  Particle p = source;
  p.SourceParticle = &source;
  p.Time += 1e-7L * units.Mpc;
  SpatialOutput output(root.string(), false, false, false, {});
  output.AddParticle(p);
  output.Flush(true);
  const auto text = Read(root / "photon");
  EXPECT_NE(text.find("normalization = first_particle_flight"),
            std::string::npos);
  const auto start = text.find("T/Mpc = ") + 8;
  const auto field = text.substr(start, text.find('\t', start) - start);
  EXPECT_EQ(field.find_first_of("eE"), std::string::npos);
  EXPECT_NEAR(Scale(text), (p.Time.t() - source.Time.t()) / units.Mpc, 1e-20);
}
} // namespace
