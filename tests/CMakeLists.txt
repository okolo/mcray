include(FetchContent)
# Pin the test-only dependency by release and archive checksum.
set(BUILD_GMOCK OFF CACHE BOOL "Build GoogleMock" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "Install GoogleTest" FORCE)
FetchContent_Declare(googletest
    URL https://codeload.github.com/google/googletest/tar.gz/refs/tags/v1.18.0
    URL_HASH SHA256=6e3191c1455468b3fc35a417fb565c1c5071aee1b7e7f85e30cf48a98d37d8b5)
FetchContent_MakeAvailable(googletest)

add_executable(mcray_unit_tests table_function_log_zeros.cpp psb_branching_bounds.cpp psb_physics.cpp nucleon_physics.cpp nuclear_gzk.cpp nuclear_photons.cpp gzk_spectral_statistics.cpp nuclear_identity.cpp talys_physics.cpp nuclear_decay.cpp shared_photon_table.cpp)
target_include_directories(mcray_unit_tests PRIVATE ../src/lib)
# Resolve main from GoogleTest before searching mcray's legacy main object.
target_link_libraries(mcray_unit_tests PRIVATE GTest::gtest_main mcray)
include(GoogleTest)

# Fixed-S secondary spectra: preserve the fast invariant rejection sampler.
add_executable(em_secondary_accuracy_tests em_secondary_accuracy.cpp)
target_include_directories(em_secondary_accuracy_tests PRIVATE ../src/lib)
target_link_libraries(em_secondary_accuracy_tests PRIVATE GTest::gtest_main mcray)
gtest_discover_tests(em_secondary_accuracy_tests TEST_PREFIX "unit.em."
    PROPERTIES TIMEOUT 30 LABELS "unit;sampling")
# The independent PPP quadrature takes about 30 s per energy on the native
# validation server. Give only these six cases additional wall-clock headroom.
set(ppp_sampling_filter "RestFrameEnergies/ProtonPppSampling.MatchesIndependentRestFrameQuadrature/*")
gtest_discover_tests(mcray_unit_tests
    TEST_PREFIX "unit."
    TEST_FILTER "-${ppp_sampling_filter}"
    DISCOVERY_TIMEOUT 30
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin"
    PROPERTIES TIMEOUT 30 LABELS "unit")
gtest_discover_tests(mcray_unit_tests
    TEST_PREFIX "unit."
    TEST_FILTER "${ppp_sampling_filter}"
    DISCOVERY_TIMEOUT 30
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin"
    PROPERTIES TIMEOUT 180 LABELS "unit")

add_test(NAME crbeam_help COMMAND CRbeam --help)
add_test(NAME z2t_smoke COMMAND z2t 0.1)
set_tests_properties(crbeam_help z2t_smoke PROPERTIES
    TIMEOUT 30 LABELS "e2e;smoke"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin")

add_executable(nuclear_scenario nuclear_scenario.cpp)
target_include_directories(nuclear_scenario PRIVATE ../src/lib)
target_link_libraries(nuclear_scenario PRIVATE mcray)

add_executable(weighted_nuclear_scenario weighted_nuclear_scenario.cpp crpropa/weighted_common.cpp)
target_include_directories(weighted_nuclear_scenario PRIVATE ../src/lib)
target_link_libraries(weighted_nuclear_scenario PRIVATE mcray)
add_executable(weighted_common_tests weighted_common_test.cpp crpropa/weighted_common.cpp)
target_link_libraries(weighted_common_tests PRIVATE GTest::gtest_main)
gtest_discover_tests(weighted_common_tests TEST_PREFIX "unit." PROPERTIES TIMEOUT 30 LABELS "unit;validation")

add_executable(crbeam_input_tests crbeam_input.cpp configuration_input.cpp
    ../src/app/crbeam/CRbeam.cpp ../src/app/crbeam/CmdLine.cpp)
target_include_directories(crbeam_input_tests PRIVATE ../src/lib ../src/app/crbeam)
target_link_libraries(crbeam_input_tests PRIVATE GTest::gtest_main mcray crbeam_configuration
    ${BOOST_SYS_LIBRARY} ${BOOST_FS_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY} ${CMAKE_DL_LIBS})
target_compile_definitions(crbeam_input_tests PRIVATE
    CRBEAM_FULL_EXAMPLE="${CMAKE_CURRENT_SOURCE_DIR}/../examples/config.full.yaml")
gtest_discover_tests(crbeam_input_tests TEST_PREFIX "integration."
    DISCOVERY_TIMEOUT 30 WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin"
    PROPERTIES TIMEOUT 30 LABELS "integration")

add_executable(nuclear_comparison nuclear_comparison.cpp)
target_link_libraries(nuclear_comparison PRIVATE GTest::gtest ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
add_test(NAME nuclear_none_run COMMAND nuclear_scenario
    --process none --batches 2 --primaries 20
    --output "${CMAKE_CURRENT_BINARY_DIR}/none.tsv")
set_tests_properties(nuclear_none_run PROPERTIES
    FIXTURES_SETUP nuclear_none TIMEOUT 30 LABELS "e2e"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin")
add_test(NAME nuclear_none_verify COMMAND nuclear_comparison
    --mc "${CMAKE_CURRENT_BINARY_DIR}/none.tsv" --batches 2 --primaries 20)
set_tests_properties(nuclear_none_verify PROPERTIES
    FIXTURES_REQUIRED nuclear_none TIMEOUT 30 LABELS "e2e")

add_executable(mcray_integration_tests psb_convolution.cpp)
target_include_directories(mcray_integration_tests PRIVATE ../src/lib)
target_link_libraries(mcray_integration_tests PRIVATE GTest::gtest_main mcray)
gtest_discover_tests(mcray_integration_tests
    TEST_PREFIX "integration." DISCOVERY_TIMEOUT 30
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin"
    PROPERTIES TIMEOUT 30 LABELS "integration")

# External real-data decay controls use GoogleTest output assertions.
add_executable(nuclear_decay_comparison nuclear_decay_comparison.cpp)
target_link_libraries(nuclear_decay_comparison PRIVATE GTest::gtest)

add_executable(nuclear_decay_probe nuclear_decay_probe.cpp)
target_include_directories(nuclear_decay_probe PRIVATE ../src/lib)
target_link_libraries(nuclear_decay_probe PRIVATE mcray)

add_executable(decay_reference_comparison decay_reference_comparison.cpp)
target_link_libraries(decay_reference_comparison PRIVATE GTest::gtest)

add_executable(photon_background_export photon_background_export.cpp)
target_include_directories(photon_background_export PRIVATE ../src/lib)
target_link_libraries(photon_background_export PRIVATE mcray)

add_executable(beta_spectrum_comparison beta_spectrum_comparison.cpp)
target_link_libraries(beta_spectrum_comparison PRIVATE GTest::gtest)

add_executable(nuclear_photon_probe nuclear_photon_probe.cpp)
target_include_directories(nuclear_photon_probe PRIVATE ../src/lib)
target_link_libraries(nuclear_photon_probe PRIVATE mcray)

# Validation-only deterministic transport coefficients; production physics unchanged.
add_executable(deterministic_pd_export deterministic_pd_export.cpp)
target_include_directories(deterministic_pd_export PRIVATE ../src/lib)
target_link_libraries(deterministic_pd_export PRIVATE mcray)
