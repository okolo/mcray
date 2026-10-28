cmake_minimum_required(VERSION 3.16)
project(mcray_crpropa_reference LANGUAGES C CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CRPROPA_SOURCE_DIR)
    message(FATAL_ERROR "Set CRPROPA_SOURCE_DIR to the frozen pinned sources")
endif()
set(DOWNLOAD_DATA OFF CACHE BOOL "Use only the pinned external data" FORCE)
set(ENABLE_PYTHON OFF CACHE BOOL "C++ reference adapter only" FORCE)
set(ENABLE_TESTING OFF CACHE BOOL "Assertions live in mcray GoogleTest" FORCE)
set(ENABLE_HDF5 OFF CACHE BOOL "TSV output" FORCE)
set(ENABLE_GIT OFF CACHE BOOL "Source revision is recorded in the manifest" FORCE)
add_subdirectory("${CRPROPA_SOURCE_DIR}" reference)
add_executable(crpropa_decay_probe decay_probe.cpp)
target_include_directories(crpropa_decay_probe PRIVATE "${CRPROPA_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/reference/include" ${CRPROPA_EXTRA_INCLUDES})
target_link_libraries(crpropa_decay_probe PRIVATE crpropa)

add_executable(crpropa_weighted_probe weighted_probe.cpp weighted_common.cpp)
target_include_directories(crpropa_weighted_probe PRIVATE "${CRPROPA_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/reference/include" ${CRPROPA_EXTRA_INCLUDES})
target_link_libraries(crpropa_weighted_probe PRIVATE crpropa)

# Link the new probe to the already-built production mcray library in the base.
option(BUILD_DOCKER_MCRAY_PROBE "Link the additional mcray probe in the Docker base" ON)
if(BUILD_DOCKER_MCRAY_PROBE)
add_executable(mcray_decay_probe mcray_decay_probe.cpp)
target_include_directories(mcray_decay_probe PRIVATE /src/src/lib)
target_link_libraries(mcray_decay_probe PRIVATE /src/build/Release/mcray/libmcray.a
    gfortran gsl gslcblas boost_system boost_filesystem gomp dl)
endif()
