cmake_minimum_required (VERSION 2.6)
project (CRbeam)

add_definitions(-DUSE_GSL)

include_directories ("../../lib" "../../OS/include")

#include_directories(../../external/include ../../external/OS/include)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG -O0 -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG -O0 -fopenmp -pg")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O2 -fopenmp")
set(CMAKE_PREFIX_PATH "external")

# -rdynamic flux is only supported on systems with ELF executable format
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -rdynamic")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -rdynamic")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -rdynamic")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -rdynamic")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

FIND_LIBRARY(C_LIBRARY c)
FIND_LIBRARY(GSL_LIBRARY gsl)
FIND_LIBRARY(GSLCBLAS_LIBRARY gslcblas)
add_subdirectory (../../lib mcray)

set(SOURCE_FILES
CmdLine.cpp
CmdLine.h
CRbeam.cpp
CRbeam.h
)

add_executable(CRbeam ${SOURCE_FILES})
target_link_libraries(CRbeam mcray ${C_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})


