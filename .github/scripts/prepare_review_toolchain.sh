#!/usr/bin/env bash
# Verify installed dependencies without configuring or executing candidate code.
set -euo pipefail
review_root="${RUNNER_TEMP:?RUNNER_TEMP must identify the scratch directory}/opencode-toolchain"
mkdir -p "$review_root/source" "$review_root/tmp" "$review_root/cache"
export TMPDIR="$review_root/tmp" XDG_CACHE_HOME="$review_root/cache"
cat > "$review_root/source/CMakeLists.txt" <<'CMAKE'
cmake_minimum_required(VERSION 3.16)
project(ReviewToolchain LANGUAGES C CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(yaml-cpp 0.6 CONFIG REQUIRED)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem system)
find_package(OpenSSL 1.1.1 REQUIRED COMPONENTS Crypto)
add_executable(toolchain-check main.cpp)
if(TARGET yaml-cpp::yaml-cpp)
  target_link_libraries(toolchain-check PRIVATE yaml-cpp::yaml-cpp)
else()
  target_link_libraries(toolchain-check PRIVATE yaml-cpp)
endif()
target_link_libraries(toolchain-check PRIVATE GSL::gsl GSL::gslcblas
  Boost::filesystem Boost::system OpenSSL::Crypto)
CMAKE
cat > "$review_root/source/main.cpp" <<'CPP'
#include <yaml-cpp/yaml.h>
#include <gsl/gsl_sf_bessel.h>
#include <boost/filesystem.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <openssl/crypto.h>
int main() {
  using boost::multiprecision::cpp_dec_float_100;
  const cpp_dec_float_100 root = sqrt(cpp_dec_float_100(2));
  return YAML::Load("value: 1")["value"].as<int>() != 1 ||
         gsl_sf_bessel_J0(0) != 1 || !boost::filesystem::exists(".") ||
         OpenSSL_version_num() == 0 || abs(root*root-2) > cpp_dec_float_100("1e-90");
}
CPP
/usr/bin/cmake --version
/usr/bin/g++ --version
/usr/bin/gfortran --version
/usr/bin/cmake -S "$review_root/source" -B "$review_root/build" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/usr \
  -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran
/usr/bin/cmake --build "$review_root/build" --parallel 2
"$review_root/build/toolchain-check"
