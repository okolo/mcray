#include "CRbeam.h"
#include "Sha256.h"
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <gtest/gtest.h>
#include <unistd.h>
namespace fs = std::filesystem;
namespace {
class ConfigurationTest : public ::testing::Test {
protected:
  fs::path root;
  std::string previousTables;
  std::vector<std::string> words;
  std::vector<char *> args;
  std::unique_ptr<CRbeam> program;
  void SetUp() override {
    root =
        fs::temp_directory_path() / ("mcray-yaml-" + std::to_string(getpid()));
    fs::remove_all(root);
    fs::create_directories(root / "configs");
    previousTables = Utils::tables_dir;
    omp_set_num_threads(1);
  }
  void TearDown() override {
    program.reset();
    Utils::tables_dir = previousTables;
    fs::remove_all(root);
  }
  fs::path Write(const std::string &text) {
    const auto path = root / "configs/input.yaml";
    std::ofstream(path) << text;
    return path;
  }
  void Parse(std::vector<std::string> options) {
    program.reset();
    words = {"CRbeam"};
    words.insert(words.end(), options.begin(), options.end());
    args.clear();
    for (auto &s : words)
      args.push_back(s.data());
    args.push_back(nullptr);
    program =
        std::make_unique<CRbeam>(static_cast<int>(words.size()), args.data());
  }
  void Error(const std::function<void()> &f, const std::string &message) {
    try {
      f();
      FAIL() << "Expected " << message;
    } catch (Utils::Exception *e) {
      const auto s = e->Message();
      delete e;
      EXPECT_NE(s.find(message), std::string::npos) << s;
    }
  }
  static std::string Read(const fs::path &file) {
    std::ifstream in(file);
    return {std::istreambuf_iterator<char>(in), {}};
  }
  std::string Common() {
    return "schema_version: 1\nsource:\n  redshift: 1.0e-8\n  particle: {type: "
           "proton}\n"
           "simulation: {primaries: 12, seed: 123456, threads: 1, batch_size: "
           "12}\n"
           "backgrounds: {ebl: {model: none}}\nmagnetic_field: {strength_G: "
           "0}\n"
           "electromagnetic: {cascade: false}\n";
  }
};
TEST_F(ConfigurationTest, RejectsDuplicatesUnknownTypesAliasesAndMixedInputs) {
  for (const auto &text :
       {"schema_version: 1\nschema_version: 1\n",
        "schema_version: 1\nsource: {redshift: .1, redshfit: .2}\n",
        "schema_version: 3\n",
        "schema_version: 1\nsource: {redshift: .1}\noutput: {directory: 123}\n",
        "schema_version: 1\nsource: {redshift: .1}\noutput: {name_suffix: "
        "true}\n",
        "schema_version: 1\nsource: {redshift: .1}\nelectromagnetic: {cascade: "
        "On}\n",
        "schema_version: 1\nsource: {redshift: .1}\nsimulation: {primaries: "
        "1.5}\n",
        "schema_version: 1\nsource: {redshift: .1}\nsimulation: {primaries: "
        "'10'}\n",
        "schema_version: 1\nsource: {redshift: .nan}\n",
        "schema_version: 1\nsource: &src {redshift: .1}\n",
        "schema_version: 1\nsource: !custom {redshift: .1}\n",
        "schema_version: 1\n---\nschema_version: 1\n",
        "schema_version: 1\nsource: {redshift: .1}\nnuclear: {data_directory: "
        "somewhere}\n",
        "schema_version: 1\nsource: {redshift: .1, particle: {type: nucleus, "
        "a: 56, z: 0}}\n",
        "schema_version: 1\nsource: {redshift: .1, spectrum: {type: power_law, "
        "index: 2}}\n"}) {
    SCOPED_TRACE(text);
    auto path = Write(text);
    Error([&] { Parse({path.string()}); }, "Configuration");
    EXPECT_FALSE(fs::exists(root / "output"));
  }
  const auto path = Write(Common());
  Error([&] { Parse({path.string(), "--seed", "1"}); }, "cannot mix");
  Error([&] { Parse({"--seed", "1", path.string()}); }, "cannot mix");
  // A CLI value ending in .yaml remains a CLI value, not positional input.
  Parse({"--output", (root / "results.yaml").string()});
}
TEST_F(ConfigurationTest, CompleteExampleAndAllBackgroundEnumsAreAccepted) {
  Parse({CRBEAM_FULL_EXAMPLE});
  EXPECT_TRUE(program->IsNuclearPrimary());
  EXPECT_EQ(program->PrimaryA(), 56u);
  EXPECT_EQ(program->NuclearPPPMode(), "Full");
  EXPECT_EQ(program->GZKResidualPolicy(), "CRPropa");
  for (const auto *model :
       {"none", "kneiske_elmag_best_fit", "kneiske_elmag_minimal",
        "inoue12_baseline", "inoue12_lower", "inoue12_upper", "franceschini08",
        "kneiske0309", "kneiske1001", "stecker05", "stecker16_lower",
        "stecker16_upper", "franceschini17"}) {
    const auto path = Write(std::string("schema_version: 1\nsource: {redshift: "
                                        ".01}\nbackgrounds: {ebl: {model: ") +
                            model + "}}\n");
    Parse({path.string()});
    EXPECT_EQ(program->EffectiveConfiguration()["backgrounds"]["ebl"]["model"]
                  .as<std::string>(),
              model);
  }
}
TEST_F(ConfigurationTest,
       CliYamlAndSavedSnapshotProduceIdenticalSeededHistories) {
  fs::path path;
  // Use a nontrivial sampled source spectrum to exercise RNG equivalence.
  std::string text = Common();
  const auto position = text.find("  particle: {type: proton}\n");
  text.insert(position + std::string("  particle: {type: proton}\n").size(),
              "  spectrum: {type: power_law, minimum_energy_eV: 1.0e+12, "
              "maximum_energy_eV: 1.0e+13, index: 2}\n");
  path = Write(text + "output: {directory: ../yaml-output}\n");
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  const auto yamlData = Read(root / "yaml-output/z0/proton");
  ASSERT_FALSE(yamlData.empty());
  auto config = YAML::LoadFile((root / "yaml-output/config.yaml").string());
  EXPECT_EQ(config["simulation"]["seed"].as<unsigned long>(), 123456ul);
  EXPECT_EQ(config["simulation"]["threads"].as<int>(), 1);
  EXPECT_EQ(config["simulation"]["batch_size"].as<int>(), 12);
  EXPECT_FALSE(fs::exists(root / "yaml-output/params"));
  const auto run = YAML::LoadFile((root / "yaml-output/run.yaml").string());
  EXPECT_EQ(run["status"].as<std::string>(), "completed");
  if (run["nuclear"].IsDefined()) {
    EXPECT_FALSE(run["nuclear"]["settings"].IsDefined());
    EXPECT_FALSE(run["nuclear"]["effective_metadata"].IsDefined());
  }
  EXPECT_FALSE(run["processes"].IsDefined());
  EXPECT_FALSE(run["input"]["argv"].IsDefined());
  EXPECT_FALSE(run["libraries"].IsDefined());
  EXPECT_EQ(
      run["config_sha256"].as<std::string>(),
      Configuration::HashFile((root / "yaml-output/config.yaml").string()));
  Parse({"--particle",
         "proton",
         "--z",
         "1e-8",
         "--emax",
         "1e13",
         "--emin_source",
         "1e12",
         "--power",
         "2",
         "--nparticles",
         "12",
         "--seed",
         "123456",
         "--batch",
         "12",
         "--background",
         "0",
         "--EGMF",
         "0",
         "--noEMcascade",
         "--output",
         (root / "cli-output").string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "cli-output/z0/proton"), yamlData);
  config["output"]["directory"] = (root / "replay").string();
  const auto replay = root / "configs/replay.yaml";
  std::ofstream(replay) << config;
  Parse({replay.string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "replay/z0/proton"), yamlData);
}
TEST_F(ConfigurationTest, LargeAndAutomaticSeedsAndTrajectoryBatchAreResolved) {
  const auto path = Write(
      Common() + "output: {directory: ../auto-output, trajectory_log: true}\n");
  auto config = YAML::LoadFile(path.string());
  config["simulation"]["seed"] = 5000000000ul;
  config["simulation"]["batch_size"] = "auto";
  std::ofstream(path) << config;
  Parse({path.string()});
  EXPECT_EQ(program->ConfiguredSeed(), 5000000000ul);
  ASSERT_EQ(program->run(), 0);
  auto effective = YAML::LoadFile((root / "auto-output/config.yaml").string());
  EXPECT_EQ(effective["simulation"]["seed"].as<unsigned long>(), 5000000000ul);
  EXPECT_EQ(effective["simulation"]["batch_size"].as<int>(), 12);
  config["simulation"]["seed"] = "auto";
  config["output"]["directory"] = (root / "random-seed").string();
  std::ofstream(path) << config;
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  effective = YAML::LoadFile((root / "random-seed/config.yaml").string());
  EXPECT_GT(effective["simulation"]["seed"].as<unsigned long>(), 0ul);
}
TEST_F(ConfigurationTest, NuclearMassesUseStandardTablesAndVerifyPayloads) {
  const auto nuclear = root / "tables/nuclear";
  fs::create_directories(nuclear);
  fs::create_directory_symlink(fs::absolute(previousTables) / "sophia2",
                               root / "tables/sophia2");
  Utils::tables_dir = (root / "tables").string() + "/";
  std::ofstream(nuclear / "provenance.txt")
      << "mcray-talys-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\ngrid_"
         "policy=crpropa-padding\n";
  std::ofstream(nuclear / "masses.tsv") << "26 13 24200 NIST\n";
  std::ofstream(nuclear / "manifest.json")
      << "{\"payload_sha256\":{\"masses.tsv\":\""
      << Configuration::HashFile((nuclear / "masses.tsv").string()) << "\"}}";
  const auto path =
      Write("schema_version: 1\nsource: {redshift: 1.0e-8, particle: {type: "
            "nucleus, a: 26, z: 13}}\n"
            "simulation: {primaries: 2, seed: 123, threads: 1}\nbackgrounds: "
            "{ebl: {model: none}}\n"
            "electromagnetic: {cascade: false}\nmagnetic_field: {strength_G: "
            "0}\noutput: {directory: ../nuclear-output}\n");
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_TRUE(fs::exists(root / "nuclear-output/z0/A26Z13"));
  const auto run = YAML::LoadFile((root / "nuclear-output/run.yaml").string());
  EXPECT_EQ(run["tables"]["nuclear_directory"].as<std::string>(),
            nuclear.string());
  EXPECT_EQ(run["tables"]["nuclear_manifests"]["manifest.json"]
               ["payload_sha256"]["masses.tsv"]
                   .as<std::string>(),
            Configuration::HashFile((nuclear / "masses.tsv").string()));
  EXPECT_FALSE(run["nuclear"]["settings"].IsDefined());
  EXPECT_FALSE(run["nuclear"]["effective_metadata"].IsDefined());
  EXPECT_FALSE(run["processes"].IsDefined());
  std::ofstream(nuclear / "masses.tsv", std::ios::app) << "27 13 25000 NIST\n";
  auto config = YAML::LoadFile(path.string());
  config["output"]["directory"] = (root / "corrupt").string();
  std::ofstream(path) << config;
  Parse({path.string()});
  Error([&] { program->run(); }, "checksum mismatch");
  EXPECT_FALSE(fs::exists(root / "corrupt"));
}
TEST_F(ConfigurationTest,
       BackgroundPathsAreRelativeToYamlAndFailuresRecordStatus) {
  auto path = Write(Common() + "output: {directory: ../failed}\n");
  auto config = YAML::LoadFile(path.string());
  config["backgrounds"]["additional"]["file"] = "missing.dat";
  std::ofstream(path) << config;
  Parse({path.string()});
  Error([&] { program->run(); }, "Cannot hash");
  EXPECT_EQ(YAML::LoadFile((root / "failed/run.yaml").string())["status"]
                .as<std::string>(),
            "failed");
  EXPECT_EQ(YAML::LoadFile((root / "failed/config.yaml")
                               .string())["backgrounds"]["additional"]["file"]
                .as<std::string>(),
            (root / "configs/missing.dat").string());
}
TEST_F(ConfigurationTest,
       MissingStandardTablesFailWithRunStatusInsteadOfCrashing) {
  fs::create_directories(root / "empty-tables");
  Utils::tables_dir = (root / "empty-tables").string() + "/";
  const auto path =
      Write(Common() + "output: {directory: ../missing-standard}\n");
  Parse({path.string()});
  Error([&] { program->run(); }, "Missing or incomplete SOPHIA");
  EXPECT_EQ(
      YAML::LoadFile((root / "missing-standard/run.yaml").string())["status"]
          .as<std::string>(),
      "failed");
}
TEST_F(ConfigurationTest, ValidAdditionalBackgroundUsesYamlRelativeFile) {
  const auto table = root / "configs/custom.dat";
  std::ofstream(table)
      << "0 0.01\n0.0001 0.001 0.001\n0.001 0.001 0.001\n0.01 0.001 0.001\n";
  const auto path = Write(Common() + "output: {directory: ../custom-output}\n");
  auto config = YAML::LoadFile(path.string());
  config["backgrounds"]["additional"]["file"] = "custom.dat";
  config["source"]["spectrum"]["type"] = "monoenergetic";
  config["source"]["spectrum"]["energy_eV"] = 1e12;
  std::ofstream(path) << config;
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  const auto run = YAML::LoadFile((root / "custom-output/run.yaml").string());
  EXPECT_EQ(run["tables"]["additional_sha256"].as<std::string>(),
            Configuration::HashFile(table.string()));
  const auto effective =
      YAML::LoadFile((root / "custom-output/config.yaml").string());
  EXPECT_EQ(effective["backgrounds"]["additional"]["file"].as<std::string>(),
            table.string());
}

std::string GridInput() {
  std::ostringstream out;
  for (int i=0;i<2;++i) for (int j=0;j<2;++j) for (int k=0;k<2;++k)
    out << i << ' ' << j << ' ' << k << " 0 1e-12 0\n";
  return out.str();
}
TEST_F(ConfigurationTest, GridCliYamlAndSnapshotAreEquivalentAndHashLoadedInput) {
  const auto grid=root/"configs/map.dat";
  std::ofstream(grid)<<GridInput();
  auto config=YAML::Load(Common()); config["schema_version"]=2;
  config["magnetic_field"]=YAML::Load("{model: grid, file: map.dat, boundary: mirror}");
  config["output"]["directory"]=(root/"grid-yaml").string();
  auto path=Write(YAML::Dump(config)); Parse({path.string()});
  ASSERT_EQ(program->run(),0);
  const auto data=Read(root/"grid-yaml/z0/proton"); ASSERT_FALSE(data.empty());
  const auto run=YAML::LoadFile((root/"grid-yaml/run.yaml").string());
  EXPECT_EQ(run["magnetic_field"]["files_sha256"][grid.string()].as<std::string>(),Configuration::HashFile(grid.string()));
  EXPECT_FALSE(run["magnetic_field"]["settings"].IsDefined());
  auto saved=YAML::LoadFile((root/"grid-yaml/config.yaml").string());
  EXPECT_EQ(saved["schema_version"].as<int>(),2);
  EXPECT_FALSE(saved["magnetic_field"]["strength_G"].IsDefined());
  EXPECT_EQ(saved["magnetic_field"]["file"].as<std::string>(),grid.string());
  Parse({"--particle","proton","--z","1e-8","--nparticles","12","--seed","123456","--batch","12",
    "--background","0","--noEMcascade","--magnetic-field","grid","--magnetic-field-file",grid.string(),
    "--magnetic-field-boundary","mirror","--output",(root/"grid-cli").string()});
  ASSERT_EQ(program->run(),0); EXPECT_EQ(Read(root/"grid-cli/z0/proton"),data);
  saved["output"]["directory"]=(root/"grid-replay").string();
  path=Write(YAML::Dump(saved)); Parse({path.string()}); ASSERT_EQ(program->run(),0);
  EXPECT_EQ(Read(root/"grid-replay/z0/proton"),data);
}
TEST_F(ConfigurationTest, VersionTwoRejectsInapplicableAndMalformedFieldSettings) {
  const auto base=YAML::Load(Common());
  for (const char* field : {
      "{model: grid}", "{model: grid, file: 42}", "{model: grid, file: map.dat, strength_G: 0}",
      "{model: grid, file: map.dat, boundary: unknown}", "{model: grid, file: map.dat, realization: per_primary}",
      "{model: grid, file: map.dat, file: duplicate}", "{model: none, modes: 1}",
      "{model: sum, components: []}", "{model: sum, components: [{model: sum, components: []}]}",
      "{model: sum, components: [{model: grid, file: map.dat, boundary: nope}]}",
      "{model: kolmogorov, modes: 1}", "{model: kolmogorov, modes: 1.5}",
      "{model: kolmogorov, minimum_scale_Mpc: 5, maximum_scale_Mpc: 1}",
      "{model: monochromatic, strength_G: .nan}", "{model: monochromatic, strength_G: '1'}",
      "{model: grid, file: map.dat, deflection_accuracy: 0}"}) {
    SCOPED_TRACE(field);
    auto path=Write(std::string("schema_version: 2\nsource: {redshift: .01}\nmagnetic_field: ")+field+"\n");
    Error([&]{Parse({path.string()});},"magnetic_field");
  }
  auto duplicate=Write("schema_version: 2\nsource: {redshift: .01}\nmagnetic_field: {model: none}\nmagnetic_field: {model: grid, file: map.dat}\n");
  Error([&]{Parse({duplicate.string()});},"duplicate");
  auto cfg=YAML::Clone(base); cfg["magnetic_field"]=YAML::Load("{model: grid, file: map.dat}");
  auto path=Write(YAML::Dump(cfg)); Error([&]{Parse({path.string()});},"Configuration");
}
TEST_F(ConfigurationTest, NewFieldCliRejectsConflictsAndBadNumbers) {
  for (auto options : std::vector<std::vector<std::string>>{
      {"--magnetic-field-file","map.dat"}, {"--magnetic-field","sum"},
      {"--magnetic-field","none","--magnetic-field","grid"},
      {"--magnetic-field","grid","--magnetic-field-file","map.dat","--magnetic-field-boundary",""},
      {"--magnetic-field","grid","--magnetic-field-file","map.dat","--EGMF","0"},
      {"--magnetic-field","grid","--magnetic-field-file","map.dat","--randomEGMF"},
      {"--magnetic-field","grid","--magnetic-field-file","map.dat","--accEGMF","10foo"},
      {"--magnetic-field","monochromatic","--turbulentEGMF"},
      {"--magnetic-field","monochromatic","--EGMF","oops"},
      {"--magnetic-field","kolmogorov","--nmodesEGMF","2.5"},
      {"--magnetic-field","none","--magnetic-field-boundary","zero"}}) {
    Error([&]{Parse(options);},"magnetic_field");
  }
  Parse({"--magnetic-field","none"});
  EXPECT_EQ(program->EffectiveConfiguration()["magnetic_field"]["model"].as<std::string>(),"none");
}
TEST_F(ConfigurationTest, GridFilesFailBeforeOutput) {
  auto cfg=YAML::Load(Common()); cfg["schema_version"]=2;
  cfg["magnetic_field"]=YAML::Load("{model: grid, file: bad.dat}");
  cfg["output"]["directory"]=(root/"invalid-grid").string();
  auto path=Write(YAML::Dump(cfg));
  Parse({path.string()}); Error([&]{program->run();},"Cannot hash");
  EXPECT_FALSE(fs::exists(root/"invalid-grid"));
  std::ofstream(root/"configs/bad.dat")<<"0 0 0 0 0 0\n";
  Parse({path.string()}); Error([&]{program->run();},"at least two");
  EXPECT_FALSE(fs::exists(root/"invalid-grid"));
}
TEST_F(ConfigurationTest, SumPreservesGridUnitsAndSupportsStochasticPerPrimaryReplay) {
  const auto grid=root/"configs/map.dat"; std::ofstream(grid)<<GridInput();
  auto cfg=YAML::Load(Common()); cfg["schema_version"]=2;
  cfg["magnetic_field"]=YAML::Load(
    "{model: sum, realization: per_primary, components: ["
    "{model: grid, file: map.dat, boundary: nearest},"
    "{model: kolmogorov, strength_G: 1e-12, minimum_scale_Mpc: 0.05, maximum_scale_Mpc: 5, modes: 8}]}");
  cfg["output"]["directory"]=(root/"sum").string();
  auto path=Write(YAML::Dump(cfg)); Parse({path.string()}); ASSERT_EQ(program->run(),0);
  const auto data=Read(root/"sum/z0/proton"); ASSERT_FALSE(data.empty());
  auto saved=YAML::LoadFile((root/"sum/config.yaml").string());
  EXPECT_EQ(saved["magnetic_field"]["components"].size(),2u);
  saved["output"]["directory"]=(root/"sum-replay").string();
  path=Write(YAML::Dump(saved)); Parse({path.string()}); ASSERT_EQ(program->run(),0);
  EXPECT_EQ(Read(root/"sum-replay/z0/proton"),data);
}
TEST_F(ConfigurationTest, ExplicitWaveModelMatchesLegacyCli) {
  const auto examples=fs::path(CRBEAM_FULL_EXAMPLE).parent_path();
  Parse({(examples/"grid-field.yaml").string()});
  EXPECT_EQ(program->EffectiveConfiguration()["magnetic_field"]["model"].as<std::string>(),"grid");
  Parse({(examples/"composed-field.yaml").string()});
  EXPECT_EQ(program->EffectiveConfiguration()["magnetic_field"]["model"].as<std::string>(),"sum");
  for (bool turbulent : {false,true}) {
    std::vector<std::string> common={"--z","1e-8","--nparticles","2","--seed","17","--background","0",
      "--noEMcascade","--EGMF","1e-12"};
    if (turbulent) {common.insert(common.end(),{"--turbulentEGMF","--nmodesEGMF","8"});}
    auto legacy=common; legacy.insert(legacy.end(),{"--output",(root/"legacy").string()});
    Parse(legacy); ASSERT_EQ(program->run(),0);
    common.insert(common.end(),{"--magnetic-field",turbulent?"kolmogorov":"monochromatic",
      "--output",(root/"explicit").string()});
    Parse(common); ASSERT_EQ(program->run(),0);
    EXPECT_EQ(Read(root/"legacy/z0/proton"),Read(root/"explicit/z0/proton"));
    fs::remove_all(root/"legacy");fs::remove_all(root/"explicit");
  }
}
TEST_F(ConfigurationTest, Sha256UsesStandardDigest) {
  EXPECT_EQ(Sha256Bytes("abc"),"ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  const auto path = Write("abc");
  EXPECT_EQ(Configuration::HashFile(path.string()),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
}
TEST_F(ConfigurationTest, Sha256HandlesPaddingAndMultipleBlocks) {
  EXPECT_EQ(Configuration::HashFile(Write(std::string(0, 'a')).string()),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  EXPECT_EQ(Configuration::HashFile(Write(std::string(55, 'a')).string()),
            "9f4390f8d30c2dd92ec9f095b65e2b9ae9b0a925a5258e241c9f1e910f734318");
  EXPECT_EQ(Configuration::HashFile(Write(std::string(56, 'a')).string()),
            "b35439a4ac6f0948b6d6f9e3c6af0f5f590ce20f1bde7090ef7970686ec6738a");
  EXPECT_EQ(Configuration::HashFile(Write(std::string(63, 'a')).string()),
            "7d3e74a05d7db15bce4ad9ec0658ea98e3f06eeecf16b4c6fff2da457ddc2f34");
  EXPECT_EQ(Configuration::HashFile(Write(std::string(64, 'a')).string()),
            "ffe054fe7ae0cb6dc65c3af9b61d5209f439851db43d0ba5997337df154668eb");
  EXPECT_EQ(Configuration::HashFile(Write(std::string(65, 'a')).string()),
            "635361c48bb9eab14198e76ea8ab7f1a41685d6ad62aa9146d301d4f17eb0ae0");
  EXPECT_EQ(Configuration::HashFile(Write(std::string(1000000, 'a')).string()),
            "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

} // namespace

TEST_F(ConfigurationTest, SpatialConeCliYamlAndSnapshotReplayAreIdentical) {
  auto config = YAML::Load(Common());
  config["schema_version"] = 2;
  config["magnetic_field"] = YAML::Load("{model: none}");
  config["source"]["position_Mpc"] = YAML::Load("[1, -2, 3]");
  config["source"]["angular_distribution"] =
      YAML::Load("{model: cone, axis: [1, 2, 3], half_angle_deg: 20}");
  config["output"] =
      YAML::Load("{directory: ../spatial-yaml, primary_direction: true}");
  const auto path = Write(YAML::Dump(config));
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  const auto data = Read(root / "spatial-yaml/z0/proton");
  EXPECT_NE(data.find("spatial_format = 1"), std::string::npos);
  EXPECT_NE(data.find("normalization = first_particle_flight"),
            std::string::npos);
  EXPECT_NE(data.find("\tPathMpc\tFlightMpc\n"), std::string::npos);
  Parse({"--particle",
         "proton",
         "--z",
         "1e-8",
         "--nparticles",
         "12",
         "--seed",
         "123456",
         "--batch",
         "12",
         "--background",
         "0",
         "--EGMF",
         "0",
         "--noEMcascade",
         "--source-position-Mpc",
         "1,-2,3",
         "--source-angular-model",
         "cone",
         "--source-axis",
         "1,2,3",
         "--source-half-angle-deg",
         "20",
         "--output-primary-direction",
         "--output",
         (root / "spatial-cli").string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "spatial-cli/z0/proton"), data);
  config = YAML::LoadFile((root / "spatial-yaml/config.yaml").string());
  config["output"]["directory"] = (root / "spatial-replay").string();
  const auto replay = Write(YAML::Dump(config));
  Parse({replay.string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "spatial-replay/z0/proton"), data);
}

TEST_F(ConfigurationTest, SpatialBirthHitAndOutsideStopProduceFiniteRows) {
  for (const std::string stop : {"inside", "outside"}) {
    auto config = YAML::Load(Common());
    config["schema_version"] = 2;
    config["magnetic_field"] = YAML::Load("{model: none}");
    config["source"]["position_Mpc"] = YAML::Load("[1, -2, 3]");
    config["observer"] = YAML::Load(
        "{model: sphere, center_Mpc: [1, -2, 3], radius_Mpc: 0.000001}");
    config["observer"]["stop_when"] = stop;
    config["output"]["directory"] = (root / stop).string();
    config["output"]["primary_direction"] = true;
    const auto path = Write(YAML::Dump(config));
    Parse({path.string()});
    ASSERT_EQ(program->run(), 0);
    const auto data = Read(root / stop / "z0/proton");
    EXPECT_NE(data.find(stop == "inside"
                            ? "normalization = fallback"
                            : "normalization = first_particle_flight"),
              std::string::npos);
    std::istringstream lines(data);
    std::string line;
    size_t count = 0;
    while (std::getline(lines, line)) {
      if (line.empty() || line[0] == '#')
        continue;
      std::istringstream fields(line);
      std::vector<long double> values;
      long double value;
      while (fields >> value)
        values.push_back(value);
      ASSERT_EQ(values.size(), 17u) << line;
      for (const auto v : values)
        EXPECT_TRUE(std::isfinite(v));
      EXPECT_EQ(values[9], 1);
      if (stop == "inside")
        EXPECT_EQ(values.back(), 0);
      else
        EXPECT_GT(values.back(), 0);
      ++count;
    }
    EXPECT_EQ(count, 12u);
  }
}

TEST_F(ConfigurationTest, SpatialEndRedshiftCapsPropagationWithoutObserver) {
  auto config = YAML::Load(Common());
  config["schema_version"] = 2;
  config["magnetic_field"] = YAML::Load("{model: none}");
  config["propagation"]["end_redshift"] = 5e-9;
  config["output"]["directory"] = (root / "end-time").string();
  const auto path = Write(YAML::Dump(config));
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  std::string data;
  for (const auto &entry : fs::directory_iterator(root / "end-time"))
    if (entry.is_directory())
      data = Read(entry.path() / "proton");
  ASSERT_FALSE(data.empty());
  const auto expected =
      (mcray::CosmoTime(5e-9).t() - mcray::CosmoTime(1e-8).t()) /
      mcray::units.Mpc;
  std::istringstream lines(data);
  std::string line;
  size_t count = 0;
  while (std::getline(lines, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream fields(line);
    long double value = 0, last = 0;
    while (fields >> value)
      last = value;
    EXPECT_NEAR(last, expected, std::abs(expected) * 1e-12);
    ++count;
  }
  EXPECT_EQ(count, 12u);
}

TEST_F(ConfigurationTest, SpatialInputValidatesVectorsRegionsAnglesAndEndTime) {
  for (const auto &options : std::vector<std::vector<std::string>>{
           {"--source-axis", "0,0,0"},
           {"--source-axis", "1,2"},
           {"--source-axis", "1,2,3,"},
           {"--source-axis", "nan,0,1"},
           {"--source-position-Mpc", "1,,3"},
           {"--source-half-angle-deg", "10"},
           {"--source-angular-model", "cone", "--source-half-angle-deg", "181"},
           {"--source-angular-model", "gaussian"},
           {"--observer-radius-Mpc", "1"},
           {"--observer", "sphere"},
           {"--observer", "sphere", "--observer-radius-Mpc", "0"},
           {"--observer", "none", "--observer-center-Mpc", "0,0,0"},
           {"--observer", "sphere", "--observer-radius-Mpc", "1",
            "--observer-stop-when", "enter"},
           {"--source-axis", "0,0,1", "--source-axis", "1,0,0"},
           {"--end-redshift", "-1"},
           {"--end-redshift", ".01"}}) {
    auto args = options;
    args.insert(args.end(), {"--z", ".01"});
    Error([&] { Parse(args); }, "");
  }
  for (const char *extra :
       {"source: {redshift: .01, position_Mpc: [1,2]}\n",
        "source: {redshift: .01, angular_distribution: {axis: [0,0,0]}}\n",
        "source: {redshift: .01, angular_distribution: {model: pencil, "
        "half_angle_deg: 10}}\n",
        "source: {redshift: .01, angular_distribution: {axis: ['1',0,0]}}\n",
        "source: {redshift: .01}\nobserver: {model: sphere, radius_Mpc: .1, "
        "radius_Mpc: .2}\n",
        "source: {redshift: .01}\nobserver: {model: none, radius_Mpc: .1}\n",
        "source: {redshift: .01}\nobserver: {model: sphere, radius_Mpc: .1, "
        "center_Mpc: [0,0]}\n",
        "source: {redshift: .01}\npropagation: {end_redshift: .01}\n",
        "source: {redshift: .01}\noutput: {primary_direction: 'true'}\n"}) {
    auto path = Write(std::string("schema_version: 2\n") + extra);
    Error([&] { Parse({path.string()}); }, "");
  }
  auto v1 = Write(
      "schema_version: 1\nsource: {redshift: .01, position_Mpc: [0,0,0]}\n");
  Error([&] { Parse({v1.string()}); }, "unknown");
}

TEST_F(ConfigurationTest,
       SpatialCliNormalizesAxisAndRecordsAllEffectiveSettings) {
  Parse({"--z",
         ".01",
         "--source-position-Mpc",
         "1, -2, 3",
         "--source-axis",
         "0,0,5",
         "--source-angular-model",
         "cone",
         "--source-half-angle-deg",
         "20",
         "--end-redshift",
         ".005",
         "--observer",
         "sphere",
         "--observer-center-Mpc",
         "0,0,1",
         "--observer-radius-Mpc",
         ".1",
         "--observer-stop-when",
         "outside",
         "--output-primary-direction"});
  const auto &cfg = program->EffectiveConfiguration();
  EXPECT_EQ(cfg["schema_version"].as<int>(), 2);
  EXPECT_EQ(cfg["source"]["position_Mpc"][1].as<double>(), -2);
  EXPECT_EQ(cfg["source"]["angular_distribution"]["axis"][2].as<double>(), 1);
  EXPECT_EQ(
      cfg["source"]["angular_distribution"]["half_angle_deg"].as<double>(), 20);
  EXPECT_EQ(cfg["observer"]["stop_when"].as<std::string>(), "outside");
  EXPECT_EQ(cfg["propagation"]["end_redshift"].as<double>(), .005);
  EXPECT_TRUE(cfg["output"]["primary_direction"].as<bool>());
}

namespace {
TEST_F(ConfigurationTest, PhotonThinningCliYamlAndSnapshotRetainExplicitPolicy) {
  auto cfg = YAML::Load(Common());
  cfg["schema_version"] = 2;
  cfg["sampling"]["thinning"]["model"] = "photon";
  cfg["sampling"]["thinning"]["alpha"] = 2.5;
  cfg["sampling"]["thinning"]["threshold_eV"] = 1e12;
  cfg["output"]["directory"] = (root / "photon-yaml").string();
  auto path = Write(YAML::Dump(cfg));
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  const auto data = Read(root / "photon-yaml/z0/proton");
  ASSERT_FALSE(data.empty());
  auto saved = YAML::LoadFile((root / "photon-yaml/config.yaml").string());
  const auto thinning = saved["sampling"]["thinning"];
  EXPECT_EQ(thinning["model"].as<std::string>(), "photon");
  EXPECT_DOUBLE_EQ(thinning["alpha"].as<double>(), 2.5);
  EXPECT_DOUBLE_EQ(thinning["threshold_eV"].as<double>(), 1e12);
  Parse({"--particle", "proton", "--z", "1e-8", "--nparticles", "12", "--seed", "123456",
         "--batch", "12", "--background", "0", "--EGMF", "0", "--noEMcascade",
         "--thinning-model", "photon", "--thinning", "2.5", "--thinning-threshold-eV", "1e12",
         "--output", (root / "photon-cli").string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "photon-cli/z0/proton"), data);
  saved["output"]["directory"] = (root / "photon-replay").string();
  path = Write(YAML::Dump(saved));
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "photon-replay/z0/proton"), data);
}
TEST_F(ConfigurationTest, PhotonThinningRejectsInvalidCliAndYamlBeforeOutput) {
  for (const auto& options : std::vector<std::vector<std::string>>{
        {"--thinning", "-1"}, {"--thinning", "2"}, {"-t", "2junk"}, {"-t", "nan"},
        {"--magnetic-field", "none", "--thinning", "0.5junk"},
        {"--thinning-model", "unknown"}, {"--thinning-model", "photon"},
        {"--thinning-threshold-eV", "1e12"},
        {"--thinning-model", "energy", "--thinning-threshold-eV", "1e12"},
        {"--thinning-model", "photon", "--thinning-threshold-eV", "-1"},
        {"--thinning-model", "photon", "--thinning-threshold-eV", "1junk"},
        {"--thinning-model", "photon", "--thinning-threshold-eV", "nan"},
        {"--thinning-model", "photon", "--thinning-threshold-eV", "1e12", "--thinning", "10.1"},
        {"--thinning-model", "photon", "--thinning-threshold-eV", "1e12", "-t", "2junk"},
        {"--thinning-model", "energy", "--thinning", "2"},
        {"--thinning-model", "energy", "--thinning-model", "photon"}}) {
    Error([&] { Parse(options); }, "Configuration");
  }
  for (const auto& settings : {"{model: photon}", "{model: unknown}",
        "{model: energy, threshold_eV: 1e12}", "{model: photon, threshold_eV: 0}",
        "{model: photon, threshold_eV: .inf}", "{model: photon, threshold_eV: 1e12, alpha: -1}",
        "{model: photon, threshold_eV: 1e12, alpha: 11}",
        "{model: photon, threshold_eV: '1e12'}"}) {
    SCOPED_TRACE(settings);
    const auto path = Write(std::string("schema_version: 2\nsampling: {thinning: ") + settings + "}\n");
    Error([&] { Parse({path.string()}); }, "Configuration");
  }
  Error([&] { Parse({Write("schema_version: 1\nsampling: {thinning: {model: photon, threshold_eV: 1e12}}\n").string()}); }, "unknown configuration key");
}
TEST_F(ConfigurationTest, ExplicitEnergyPolicyPreservesLegacySeededOutput) {
  std::vector<std::string> options{"--particle", "electron", "--z", "1e-3", "--emax", "1e13",
      "--emin", "1e10", "--nparticles", "20", "--seed", "723", "--batch", "20",
      "--background", "0", "--EGMF", "0", "--output", (root / "legacy-energy").string()};
  Parse(options);
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(program->EffectiveConfiguration()["schema_version"].as<int>(), 1);
  const auto legacy = Read(root / "legacy-energy/z0/photon");
  ASSERT_FALSE(legacy.empty());
  options.back() = (root / "explicit-energy").string();
  options.insert(options.end(), {"--thinning-model", "energy"});
  Parse(options);
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "explicit-energy/z0/photon"), legacy);
}
} // namespace

namespace {
TEST_F(ConfigurationTest, PhotonThinningIsActiveInRealCascadeAndReplaysAcrossInterfaces) {
  auto cfg = YAML::Load(Common());
  cfg["schema_version"] = 2;
  cfg["source"]["redshift"] = .001;
  cfg["source"]["particle"]["type"] = "electron";
  cfg["source"]["spectrum"]["type"] = "monoenergetic";
  cfg["source"]["spectrum"]["energy_eV"] = 1e13;
  // Four independent primaries suffice for this interface/replay check; the
  // separate Bernoulli moment test supplies the statistical ensemble.
  cfg["simulation"]["primaries"] = 4;
  cfg["simulation"]["batch_size"] = 4;
  cfg["simulation"]["minimum_energy_eV"] = 1e10;
  cfg["electromagnetic"]["cascade"] = true;
  cfg["sampling"]["thinning"]["model"] = "photon";
  cfg["sampling"]["thinning"]["threshold_eV"] = 1e12;
  cfg["sampling"]["thinning"]["alpha"] = 0;
  cfg["output"]["directory"] = (root / "no-photon-thinning").string();
  auto path = Write(YAML::Dump(cfg));
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  const auto unthinned = Read(root / "no-photon-thinning/z0/photon");
  ASSERT_FALSE(unthinned.empty());
  cfg["sampling"]["thinning"]["alpha"] = 1;
  cfg["output"]["directory"] = (root / "active-photon").string();
  path = Write(YAML::Dump(cfg));
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  const auto active = Read(root / "active-photon/z0/photon");
  ASSERT_FALSE(active.empty());
  EXPECT_NE(active, unthinned);
  std::istringstream rows(active);
  std::string line;
  bool corrected = false;
  while (std::getline(rows, line)) {
    if (line.empty() || line[0] == '#') continue;
    double energy, weight;
    ASSERT_TRUE(bool(std::istringstream(line) >> energy >> weight));
    corrected = corrected || weight > 1;
  }
  EXPECT_TRUE(corrected);
  Parse({"--particle", "electron", "--z", ".001", "--emax", "1e13", "--emin", "1e10",
         "--nparticles", "4", "--seed", "123456", "--batch", "4", "--background", "0", "--EGMF", "0",
         "--thinning-model", "photon", "--thinning", "1", "--thinning-threshold-eV", "1e12",
         "--output", (root / "active-photon-cli").string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "active-photon-cli/z0/photon"), active);
  auto saved = YAML::LoadFile((root / "active-photon/config.yaml").string());
  saved["output"]["directory"] = (root / "active-photon-replay").string();
  path = Write(YAML::Dump(saved));
  Parse({path.string()});
  ASSERT_EQ(program->run(), 0);
  EXPECT_EQ(Read(root / "active-photon-replay/z0/photon"), active);
  cfg["sampling"]["thinning"]["alpha"] = 10;
  cfg["sampling"]["thinning"]["threshold_eV"] = 1e300;
  cfg["output"]["directory"] = (root / "failed-photon").string();
  path = Write(YAML::Dump(cfg));
  Parse({path.string()});
  Error([&] { program->run(); }, "Propagation failed; partial output must not be used");
  const auto failed = YAML::LoadFile((root / "failed-photon/run.yaml").string());
  EXPECT_EQ(failed["status"].as<std::string>(), "failed");
}
} // namespace
