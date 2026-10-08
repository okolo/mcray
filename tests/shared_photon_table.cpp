#include "shared_photon_table.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <unistd.h>
TEST(SharedPhotonTable, PhysicalUnitsLogInterpolationAndStaticMedium) {
    auto file=std::filesystem::temp_directory_path()/("mcray-shared-photons-"+std::to_string(getpid()));
    std::ofstream(file)<<"# E/eV dn/dE/(cm^-3 eV^-1)\n1 10\n4 40\n16 0\n";
    SharedPhotonTable field(file.string());std::filesystem::remove(file);
    EXPECT_NEAR(field.n(2*mcray::units.eV,0.)*mcray::units.cm3,40.,1e-12);
    EXPECT_DOUBLE_EQ(field.n(2*mcray::units.eV,0.),field.n(2*mcray::units.eV,0.01));
    EXPECT_NEAR(field.n(8*mcray::units.eV,0.)*mcray::units.cm3,160.,1e-10);
    EXPECT_DOUBLE_EQ(field.n(0.5*mcray::units.eV,0.),0.);
    EXPECT_DOUBLE_EQ(field.n(17*mcray::units.eV,0.),0.);
}
TEST(SharedPhotonTable, NonmonotonicInputIsRejected) {
    auto file=std::filesystem::temp_directory_path()/("mcray-shared-photons-bad-"+std::to_string(getpid()));
    std::ofstream(file)<<"1 10\n1 20\n";
    try {SharedPhotonTable field(file.string());FAIL()<<"Expected invalid grid";}
    catch(Utils::Exception* error){EXPECT_NE(error->Message().find("Invalid shared"),std::string::npos);delete error;}
    std::filesystem::remove(file);
}
