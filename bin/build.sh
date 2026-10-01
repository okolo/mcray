#!/bin/bash
# Custom build script for Mac OS

set -euo pipefail

if [[ $# -ne 1 || ( "$1" != "Debug" && "$1" != "Release" ) ]]
then
    echo "usage: build.sh Debug | Release"
    exit 1
fi

find_tool()
{
    local tool_name
    for tool_name in "$@"
    do
        if command -v "$tool_name" >/dev/null 2>&1
        then
            command -v "$tool_name"
            return 0
        fi
    done
    return 1
}

build_config="$1"
script_dir="$(cd "$(dirname "$0")" && pwd)"
source_dir="$script_dir/../src/app/crbeam"
build_dir="$script_dir/../build/$build_config"

cmake_options=("-DCMAKE_CXX_STANDARD=14" "-DCMAKE_CXX_STANDARD_REQUIRED=ON")
if [[ "$(uname -s)" == "Darwin" ]]
then
    brew_prefix="$(brew --prefix)"
    llvm_prefix="$(brew --prefix llvm)"
    # Homebrew Boost uses libc++; select a matching C++ compiler.
    c_compiler="$llvm_prefix/bin/clang"
    cxx_compiler="$llvm_prefix/bin/clang++"
    if [[ ! -x "$c_compiler" || ! -x "$cxx_compiler" ]]
    then
        echo "error: Homebrew LLVM was not found" >&2
        exit 1
    fi
    # Keep Conda's assembler from shadowing the system tools used by gfortran.
    export PATH="$brew_prefix/bin:/usr/bin:/bin:$PATH"
    gcc_prefix="$(brew --prefix gcc)"
    cmake_options+=(
        "-DCMAKE_C_FLAGS=-I$brew_prefix/include"
        "-DCMAKE_CXX_FLAGS=-I$brew_prefix/include"
        "-DCMAKE_EXE_LINKER_FLAGS=-L$brew_prefix/lib -L$llvm_prefix/lib -L$gcc_prefix/lib/gcc/current"
        "-DCMAKE_LIBRARY_PATH=$brew_prefix/lib"
    )
else
    c_compiler="$(find_tool gcc-15 gcc-14 gcc)" || {
        echo "error: GCC was not found" >&2
        exit 1
    }
    cxx_compiler="$(find_tool g++-15 g++-14 g++)" || {
        echo "error: G++ was not found" >&2
        exit 1
    }
fi
fortran_compiler="$(find_tool gfortran-15 gfortran-14 gfortran)" || {
    echo "error: GFortran was not found" >&2
    exit 1
}

echo "Building $build_config with:"
echo "  C:       $c_compiler"
echo "  C++:     $cxx_compiler"
echo "  Fortran: $fortran_compiler"

cmake --fresh \
    -S "$source_dir" \
    -B "$build_dir" \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    -DCMAKE_C_COMPILER="$c_compiler" \
    -DCMAKE_CXX_COMPILER="$cxx_compiler" \
    -DCMAKE_Fortran_COMPILER="$fortran_compiler" \
    -DCMAKE_BUILD_TYPE="$build_config" \
    "${cmake_options[@]}"

cmake --build "$build_dir" --parallel "${MCRAY_BUILD_JOBS:-2}"

echo "Build output: $build_dir"
