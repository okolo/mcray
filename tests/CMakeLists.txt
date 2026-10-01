include(FetchContent)
# Pin the test-only dependency by release and archive checksum.
set(BUILD_GMOCK OFF CACHE BOOL "Build GoogleMock" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "Install GoogleTest" FORCE)
FetchContent_Declare(googletest
    URL https://codeload.github.com/google/googletest/tar.gz/refs/tags/v1.18.0
    URL_HASH SHA256=6e3191c1455468b3fc35a417fb565c1c5071aee1b7e7f85e30cf48a98d37d8b5)
FetchContent_MakeAvailable(googletest)

add_executable(mcray_unit_tests table_function_log_zeros.cpp psb_branching_bounds.cpp)
target_include_directories(mcray_unit_tests PRIVATE ../src/lib)
# Resolve main from GoogleTest before searching mcray's legacy main object.
target_link_libraries(mcray_unit_tests PRIVATE GTest::gtest_main mcray)
include(GoogleTest)
gtest_discover_tests(mcray_unit_tests
    TEST_PREFIX "unit."
    DISCOVERY_TIMEOUT 30
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin"
    PROPERTIES TIMEOUT 30 LABELS "unit")

add_test(NAME crbeam_help COMMAND CRbeam --help)
add_test(NAME z2t_smoke COMMAND z2t 0.1)
set_tests_properties(crbeam_help z2t_smoke PROPERTIES
    TIMEOUT 30 LABELS "e2e;smoke"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../bin")
