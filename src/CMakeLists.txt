cmake_minimum_required (VERSION 2.6)
project (CRbeam)

add_definitions(-DUSE_GSL)

include_directories ("${PROJECT_SOURCE_DIR}/MCPropag")

include_directories(MCPropag/external/include MCPropag/OS/include)
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
find_library(FFTW_LIBRARIES NAMES fftw3 fftw3-3)
find_path(FFTW_INCLUDE_DIRS NAMES fftw3.h)
find_library(FFTW_OMP_LIBRARY NAMES fftw3_omp)
if(NOT FFTW_OMP_LIBRARY)
    message(FATAL_ERROR "FFTW OpenMP library (fftw3_omp) not found. Please install/compile FFTW with --enable-openmp.")
endif()

add_subdirectory (MCPropag)

set(SOURCE_FILES
CmdLine.cpp
CmdLine.h
CRbeam.cpp
CRbeam.h
)

add_executable(CRbeam ${SOURCE_FILES})
#target_include_directories(CRbeam PRIVATE ${FFTW_INCLUDE_DIRS})
target_link_libraries(CRbeam MCPropag ${CMAKE_CXX_FLAGS_RELEASE} ${FFTW_OMP_LIBRARY} ${FFTW_LIBRARIES} ${C_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
