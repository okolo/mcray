cmake_minimum_required (VERSION 3.1.0)
project (utils)

add_definitions(-DUSE_GSL)
include_directories(external/include)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG -O0 -pg")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O2")
set(CMAKE_PREFIX_PATH "external")

#FIND_LIBRARY(C_LIBRARY c)
#FIND_LIBRARY(GSL_LIBRARY gsl)
#FIND_LIBRARY(GSLCBLAS_LIBRARY gslcblas)
#FIND_PACKAGE(HDF5 REQUIRED)


set(SOURCE_FILES
CmdLine.cpp
CmdLine.h
Debug.cpp
Debug.h
Deflection3D.cpp
Deflection3D.h
MathUtils.cpp
MathUtils.h
Randomizer.cpp
Randomizer.h
StringUtils.cpp
StringUtils.h
TableFunction.cpp
TableFunction.h
TableReader.cpp
TableReader.h
Utils.cpp
Utils.h
Units.cpp
Units.h
Jet.cpp
Jet.h
)

add_library(utils ${SOURCE_FILES})

#target_link_libraries(utils ${C_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY} ${HDF5_CXX_LIBRARIES})
