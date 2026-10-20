#include "CRbeam.h"
#include <gtest/gtest.h>
#include <functional>
#include <filesystem>
#include <fstream>
#include <unistd.h>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {
class CRbeamInputTest : public ::testing::Test {
protected:
    std::vector<std::string> words;
    std::vector<char*> argv;
    std::unique_ptr<CRbeam> program;

    void Parse(std::initializer_list<std::string> options) {
        program.reset();
        words = {"CRbeam"};
        words.insert(words.end(), options.begin(), options.end());
        argv.clear();
        for (auto& word : words) argv.push_back(word.data());
        const int argc = static_cast<int>(argv.size());
        argv.push_back(nullptr);
        program = std::make_unique<CRbeam>(argc, argv.data());
    }
    void ExpectError(const std::function<void()>& action, const std::string& text) {
        try {
            action();
            FAIL() << "Expected input error containing: " << text;
        } catch (Utils::Exception* error) {
            const std::string message = error->Message();
            delete error;
            EXPECT_NE(message.find(text), std::string::npos) << message;
        }
    }
};

TEST_F(CRbeamInputTest, MissingPDPolicyIsExplicitAndIndependentOfPhotonYields) {
    Parse({}); EXPECT_EQ(program->PDMissingDataPolicy(),"Strict");
    Parse({"--particle","nucleus","--A","56","--photodisintegration","TALYS",
           "--talys-data","unused","--pd-missing-data","CRPropa","--pd-photon-yields","Strict"});
    EXPECT_EQ(program->PDMissingDataPolicy(),"CRPropa");
    EXPECT_EQ(program->PDPhotonYieldPolicy(),"Strict");
    ExpectError([&] { Parse({"--pd-missing-data","unknown"}); }, "--pd-missing-data must be");
    ExpectError([&] { Parse({"--pd-missing-data","CRPropa"}); }, "requires a nuclear TALYS primary");
}
TEST_F(CRbeamInputTest, ExistingDefaultRemainsNonNuclear) {
    Parse({});
    EXPECT_FALSE(program->IsNuclearPrimary());
    EXPECT_EQ(program->PrimaryA(), 0u);
    EXPECT_EQ(program->PrimaryZ(), 0u);
}
TEST_F(CRbeamInputTest, AutomaticChargeUsesPsbRepresentative) {
    Parse({"--particle", "nucleus", "--A", "56"});
    EXPECT_TRUE(program->IsNuclearPrimary());
    EXPECT_EQ(program->PrimaryA(), 56u);
    EXPECT_EQ(program->PrimaryZ(), 26u);
}
TEST_F(CRbeamInputTest, ExplicitZeroChargeSelectsAutomatically) {
    Parse({"--particle", "nucleus", "--A", "56", "--Z", "0"});
    EXPECT_EQ(program->PrimaryZ(), 26u);
}
TEST_F(CRbeamInputTest, ExplicitChargeAndShortAliasesAreAccepted) {
    Parse({"-particle", "nucleus", "-A", "56", "-Z", "26"});
    EXPECT_EQ(program->PrimaryA(), 56u);
    EXPECT_EQ(program->PrimaryZ(), 26u);
}
TEST_F(CRbeamInputTest, NonNuclearInputsIgnoreNuclearValues) {
    const auto max = std::to_string(std::numeric_limits<unsigned int>::max());
    Parse({"--particle", "photon", "--A", max, "--Z", max});
    EXPECT_FALSE(program->IsNuclearPrimary());
    EXPECT_EQ(program->PrimaryA(), 0u);
    EXPECT_EQ(program->PrimaryZ(), 0u);
    Parse({"--primary", "10", "--A", "1", "--Z", "2"});
    EXPECT_FALSE(program->IsNuclearPrimary());
}
TEST_F(CRbeamInputTest, NucleusRequiresMassNumberAtLeastTwo) {
    for (const auto& a : {"0", "1"})
        ExpectError([&] { Parse({"--particle", "nucleus", "--A", a}); }, "--A >= 2");
    ExpectError([&] { Parse({"--particle", "nucleus"}); }, "--A >= 2");
}
TEST_F(CRbeamInputTest, ChargeCannotExceedMassNumber) {
    ExpectError([&] { Parse({"--particle", "nucleus", "--A", "56", "--Z", "57"}); }, "--Z must not exceed --A");
}
TEST_F(CRbeamInputTest, MissingPsbEntryRequiresExplicitCharge) {
    ExpectError([&] { Parse({"--particle", "nucleus", "--A", "57"}); }, "specify --Z explicitly");
}
TEST_F(CRbeamInputTest, DifferentPsbChainPrimariesAreAccepted) {
    for (const auto& isotope : std::vector<std::pair<unsigned,unsigned>>{{2,1},{4,2},{12,6},{16,8},{55,25}}) {
        Parse({"--particle", "nucleus", "--A", std::to_string(isotope.first)});
        EXPECT_EQ(program->PrimaryA(),isotope.first);
        EXPECT_EQ(program->PrimaryZ(),isotope.second);
        Parse({"--particle", "nucleus", "--A", std::to_string(isotope.first),
               "--photodisintegration", "PSB"});
        EXPECT_EQ(program->PrimaryA(),isotope.first);
    }
}
TEST_F(CRbeamInputTest, PsbRejectsUnsupportedIdentityWithoutProjection) {
    for (const auto& isotope : std::vector<std::pair<unsigned,unsigned>>{{56,25},{57,26}})
        ExpectError([&] { Parse({"--particle", "nucleus", "--A", std::to_string(isotope.first),
             "--Z", std::to_string(isotope.second), "--photodisintegration", "PSB"}); }, "PSB requires its representative isotope");
}
TEST_F(CRbeamInputTest, TalysInputAcceptsExplicitIsotopesAndChecksDataAtRunTime) {
    for (const auto& isotope : std::vector<std::pair<unsigned,unsigned>>{{12,6},{56,25},{57,26},{238,92}}) {
        Parse({"--particle", "nucleus", "--A", std::to_string(isotope.first),
               "--Z", std::to_string(isotope.second), "--photodisintegration", "TALYS", "--talys-data", "unused"});
        EXPECT_EQ(program->PrimaryA(),isotope.first);
        EXPECT_EQ(program->PrimaryZ(),isotope.second);
    }
}
TEST_F(CRbeamInputTest, OffDoesNotSilentlySubstituteAnUnsupportedIsotope) {
    for (const auto& isotope : std::vector<std::pair<unsigned,unsigned>>{{26,13},{56,25},{238,92}}) {
        Parse({"--particle", "nucleus", "--A", std::to_string(isotope.first),
               "--Z", std::to_string(isotope.second), "--photodisintegration", "Off"});
    }
}
TEST_F(CRbeamInputTest, UnsignedInputsRejectSignsFractionsJunkAndOverflow) {
    const auto overflow = std::to_string(static_cast<unsigned long long>(std::numeric_limits<unsigned int>::max()) + 1);
    for (const auto& value : std::vector<std::string>{"-1", "+1", "1.5", "1e2", "56junk", "", " 56", overflow}) {
        SCOPED_TRACE(value);
        // Syntax is checked even when the physical A,Z values are unused.
        ExpectError([&] { Parse({"--A", value}); }, "Invalid unsigned integer for --A");
        ExpectError([&] { Parse({"--Z", value}); }, "Invalid unsigned integer for --Z");
    }
}
TEST_F(CRbeamInputTest, MissingValuesDoNotUseDefaultsOrReadPastArgv) {
    ExpectError([&] { Parse({"--A"}); }, "Missing argument for --A");
    ExpectError([&] { Parse({"--Z", "--A", "56"}); }, "Missing argument for --Z");
    ExpectError([&] { Parse({"--particle"}); }, "Missing argument for --particle");
}
TEST_F(CRbeamInputTest, UnknownOrConflictingSelectorsAreErrors) {
    ExpectError([&] { Parse({"--particle", "unknown"}); }, "Unknown --particle value");
    ExpectError([&] { Parse({"--particle", "nucleus", "--primary", "10", "--A", "56"}); }, "Use either --particle or --primary");
}
TEST_F(CRbeamInputTest, EveryExistingNamedParticleIsAccepted) {
    for (int type = Electron; type <= Proton; ++type) {
        Parse({"--particle", Particle::Name(static_cast<ParticleType>(type))});
        EXPECT_FALSE(program->IsNuclearPrimary());
    }
}
TEST_F(CRbeamInputTest, PhotodisintegrationDefaultsToOff) {
    Parse({"--particle", "nucleus", "--A", "56"});
    EXPECT_EQ(program->PhotodisintegrationMode(), "Off");
}
TEST_F(CRbeamInputTest, PhotodisintegrationModeIsValidated) {
    Parse({"--particle", "nucleus", "--A", "56", "--photodisintegration", "PSB"});
    EXPECT_EQ(program->PhotodisintegrationMode(), "PSB");
    ExpectError([&] { Parse({"--photodisintegration", "unknown"}); }, "must be TALYS, PSB or Off");
    Parse({"--particle", "nucleus", "--A", "56", "--photodisintegration", "TALYS"});
}
TEST_F(CRbeamInputTest, NuclearOffAndPsbRunsProduceIsotopeFiles) {
    omp_set_num_threads(1);
    for (const std::string mode : {"Off", "PSB"}) {
        const auto dir = std::filesystem::temp_directory_path() / ("mcray-crbeam-" + mode + "-" + std::to_string(getpid()));
        std::filesystem::remove_all(dir);
        Parse({"--particle", "nucleus", "--A", "56", "--photodisintegration", mode,
               "--z", "0.0001", "--emax", "1e21", "--nparticles", "8",
               "--noEMcascade", "--background", "0", "--output", dir.string()});
        ASSERT_EQ(program->run(), 0);
        bool found = false;
        for (const auto& entry : std::filesystem::directory_iterator(dir / "z0")) {
            if (entry.path().filename().string().front() == 'A') found = true;
            EXPECT_NE(entry.path().filename().string(), "nucleus");
        }
        EXPECT_TRUE(found);
        if (mode == "Off") {
            std::ifstream file(dir / "z0" / "A56Z26"); ASSERT_TRUE(file.good());
            unsigned rows=0;
            for (std::string line; std::getline(file,line);) if (!line.empty() && line[0]!='#') ++rows;
            EXPECT_EQ(rows,8u);
        }
        std::filesystem::remove_all(dir);
    }
}
TEST_F(CRbeamInputTest, OffPropagatesExplicitLightAndHeavyIsotopesWithMassOnlyData) {
    omp_set_num_threads(1);
    const auto root=std::filesystem::temp_directory_path() / ("mcray-crbeam-masses-"+std::to_string(getpid()));
    std::filesystem::remove_all(root);std::filesystem::create_directory(root);
    std::ofstream(root/"provenance.txt")<<"mcray-talys-v1\n34b10e27375c0e9894d7c984ac15fb96e12975ef\ngrid_policy=crpropa-padding\n";
    std::ofstream(root/"manifest.json")<<"{}\n";
    std::ofstream(root/"masses.tsv")<<"26 13 24200 NIST\n238 92 221700 NIST\n";
    for (const auto& isotope : std::vector<std::pair<unsigned,unsigned>>{{26,13},{238,92}}) {
        const std::string name="A"+std::to_string(isotope.first)+"Z"+std::to_string(isotope.second);
        const auto output=root/name;
        Parse({"--particle","nucleus","--A",std::to_string(isotope.first),"--Z",std::to_string(isotope.second),
               "--talys-data",root.string(),"--photodisintegration","Off","--z","0.0001",
               "--emax","1e21","--nparticles","8","--noEMcascade","--background","0","--output",output.string()});
        ASSERT_EQ(program->run(),0);
        std::ifstream file(output/"z0"/name);ASSERT_TRUE(file.good());
        unsigned rows=0;
        for (std::string line;std::getline(file,line);) if (!line.empty() && line[0]!='#') ++rows;
        EXPECT_EQ(rows,8u);
    }
    Parse({"--particle","nucleus","--A","239","--Z","92","--talys-data",root.string()});
    ExpectError([&] { program->run(); },"Missing TALYS isotope");
    Parse({"--particle","nucleus","--A","26","--Z","13","--talys-data",root.string(),"--photodisintegration","TALYS"});
    ExpectError([&] { program->run(); },"Cannot read TALYS data");
    std::filesystem::remove_all(root);
}
TEST_F(CRbeamInputTest, NuclearDecayRequiresExplicitCompatibleInput) {
    ExpectError([&] { Parse({"--nuclear-decay","unknown"}); }, "must be On or Off");
    ExpectError([&] { Parse({"--nuclear-decay","On"}); }, "requires a nuclear primary");
    Parse({"--particle","nucleus","--A","56","--nuclear-decay","On"});
    ExpectError([&] { Parse({"--particle","nucleus","--A","56","--photodisintegration","PSB","--nuclear-decay","On","--talys-data","unused"}); }, "photodisintegration TALYS or Off");
    Parse({"--particle","nucleus","--A","56","--nuclear-decay","On","--talys-data","unused"});
}
} // namespace

namespace {
TEST_F(CRbeamInputTest, NuclearPppControlsValidateMassesCascadeAndDefaults) {
    Parse({}); EXPECT_EQ(program->NuclearPPPMode(),"Off");
    ExpectError([&]{Parse({"--nuclear-ppp","On"});},"must be Off, Losses or Full");
    ExpectError([&]{Parse({"--nuclear-ppp","Losses"});},"requires a nuclear primary");
    Parse({"--particle","nucleus","--A","56","--nuclear-ppp","Losses"});
    ExpectError([&]{Parse({"--particle","nucleus","--A","56","--talys-data","unused","--photodisintegration","PSB","--nuclear-ppp","Losses"});},"photodisintegration TALYS or Off");
    ExpectError([&]{Parse({"--particle","nucleus","--A","56","--talys-data","unused","--nuclear-ppp","Full","--noEMcascade"});},"requires EM cascade");
    ExpectError([&]{Parse({"--particle","nucleus","--A","238","--Z","92","--talys-data","unused","--nuclear-ppp","Losses"});},"A<=56");
    for (const auto mode : {"Losses","Full"}) {
        Parse({"--particle","nucleus","--A","56","--talys-data","unused","--nuclear-ppp",mode});
        EXPECT_EQ(program->NuclearPPPMode(),mode);
    }
    Parse({"--particle","nucleus","--A","56","--talys-data","unused","--nuclear-ppp","Losses","--noEMcascade"});
    EXPECT_EQ(program->NuclearPPPMode(),"Losses");
}
}

namespace {
TEST_F(CRbeamInputTest, NuclearGzkPolicyDefaultsAndValidation) {
    Parse({}); EXPECT_EQ(program->NuclearGZKMode(),"Off");EXPECT_EQ(program->GZKResidualPolicy(),"Strict");
    ExpectError([&]{Parse({"--nuclear-gzk","Full"});},"must be On or Off");
    ExpectError([&]{Parse({"--nuclear-gzk-residuals","1"});},"must be Strict or CRPropa");
    ExpectError([&]{Parse({"--nuclear-gzk","On"});},"requires a nuclear primary");
    Parse({"--particle","nucleus","--A","56","--nuclear-gzk","On"});
    ExpectError([&]{Parse({"--particle","nucleus","--A","56","--talys-data","unused","--nuclear-gzk","On","--nuclear-gzk-residuals","CRPropa"});},"--nuclear-decay On");
    ExpectError([&]{Parse({"--nuclear-gzk-residuals","CRPropa"});},"--nuclear-gzk On");
    for (const auto policy:{"Strict","CRPropa"}) {
        Parse({"--particle","nucleus","--A","56","--talys-data","unused","--nuclear-gzk","On","--nuclear-decay","On","--nuclear-gzk-residuals",policy});
        EXPECT_EQ(program->NuclearGZKMode(),"On");EXPECT_EQ(program->GZKResidualPolicy(),policy);
    }
}
}

namespace {
TEST_F(CRbeamInputTest, NuclearPhotonControlsAndYieldPolicy) {
    Parse({});EXPECT_EQ(program->PDPhotonMode(),"Off");EXPECT_EQ(program->DecayPhotonMode(),"Off");EXPECT_EQ(program->PDPhotonYieldPolicy(),"CRPropa");
    ExpectError([&]{Parse({"--photodisintegration-photons","invalid"});},"must be On or Off");
    ExpectError([&]{Parse({"--decay-photons","On"});},"--nuclear-decay On");
    ExpectError([&]{Parse({"--photodisintegration-photons","On"});},"nuclear TALYS primary");
    ExpectError([&]{Parse({"--pd-photon-yields","invalid"});},"must be Strict or CRPropa");
    for(const auto policy:{"Strict","CRPropa"}) {
        Parse({"--particle","nucleus","--A","56","--photodisintegration","TALYS","--talys-data","unused",
               "--nuclear-decay","On","--photodisintegration-photons","On","--decay-photons","On","--pd-photon-yields",policy,"--noEMcascade"});
        EXPECT_EQ(program->PDPhotonYieldPolicy(),policy);EXPECT_EQ(program->PDPhotonMode(),"On");
    }
}
}

namespace {
TEST_F(CRbeamInputTest, ExplicitSeedIsUnsignedAndOptional) {
    Parse({});EXPECT_FALSE(program->HasExplicitSeed());
    Parse({"--seed","0"});EXPECT_TRUE(program->HasExplicitSeed());EXPECT_EQ(program->ConfiguredSeed(),0u);
    Parse({"-seed","20261011"});EXPECT_TRUE(program->HasExplicitSeed());EXPECT_EQ(program->ConfiguredSeed(),20261011u);
    ExpectError([&]{Parse({"--seed","-1"});},"seed");
    ExpectError([&]{Parse({"--seed","1x"});},"seed");
}
}
