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
TEST_F(CRbeamInputTest, UnsupportedPrimariesAreNotSilentlyProjected) {
    ExpectError([&] { Parse({"--particle", "nucleus", "--A", "55"}); }, "only primary Fe56");
    ExpectError([&] { Parse({"--particle", "nucleus", "--A", "56", "--Z", "25"}); }, "only primary Fe56");
    ExpectError([&] { Parse({"--particle", "nucleus", "--A", "57", "--Z", "26"}); }, "only primary Fe56");
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
    ExpectError([&] { Parse({"--particle", "nucleus", "--A", "56", "--photodisintegration", "TALYS"}); }, "requires --talys-data");
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
TEST_F(CRbeamInputTest, NuclearDecayRequiresExplicitCompatibleInput) {
    ExpectError([&] { Parse({"--nuclear-decay","unknown"}); }, "must be On or Off");
    ExpectError([&] { Parse({"--nuclear-decay","On"}); }, "requires a nuclear primary");
    ExpectError([&] { Parse({"--particle","nucleus","--A","56","--nuclear-decay","On"}); }, "--talys-data");
    ExpectError([&] { Parse({"--particle","nucleus","--A","56","--photodisintegration","PSB","--nuclear-decay","On","--talys-data","unused"}); }, "photodisintegration TALYS or Off");
    Parse({"--particle","nucleus","--A","56","--nuclear-decay","On","--talys-data","unused"});
}
} // namespace
