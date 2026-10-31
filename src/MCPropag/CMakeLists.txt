cmake_minimum_required (VERSION 2.6)
project (MCPropag)
enable_language (Fortran)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS "-ffixed-line-length-200 -fno-f2c")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-fno-f2c -O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


add_definitions(-DUSE_GSL)
include_directories(external/include OS/include)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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

#FIND_LIBRARY(FORTRAN_LIBRARY gfortran PATHS OS/lib)
# FIND_LIBRARY(OMP_LIBRARY gomp PATHS OS/lib)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
FIND_LIBRARY(C_LIBRARY c)
FIND_LIBRARY(XERCES_LIBRARY xerces-c)
FIND_LIBRARY(GSL_LIBRARY gsl)
FIND_LIBRARY(GSLCBLAS_LIBRARY gslcblas)
find_library(FFTW_LIBRARIES NAMES fftw3 fftw3-3)
find_library(FFTW_OMP_LIBRARY NAMES fftw3_omp)
if(NOT FFTW_OMP_LIBRARY)
    message(FATAL_ERROR "FFTW OpenMP library (fftw3_omp) not found. Please install/compile FFTW with --enable-openmp.")
endif()
find_path(FFTW_INCLUDE_DIRS NAMES fftw3.h)

set(SOURCE_FILES
SOPHIA/eventgen.f
SOPHIA/inpoutput.f
SOPHIA/jetset74dp.f
SOPHIA/rndm.f
SOPHIA/sampling.f
sophia2cpp.f
Sophia.cpp
Sophia.h
Background.cpp
Background.h
Cosmology.cpp
Cosmology.h
Debug.cpp
Debug.h
Deflection1D.cpp
Deflection1D.h
Deflection3D.cpp
Deflection3D.h
Filter.cpp
Filter.h
GammaPP.cpp
GammaPP.h
GZK.cpp
GZK.h
ICS.cpp
ICS.h
Inoue12IROSpectrum.cpp
Inoue12IROSpectrum.h
Interaction.cpp
Interaction.h
MathUtils.cpp
MathUtils.h
NeutronDecay.cpp
NeutronDecay.h
Nucleus.cpp
Nucleus.h
Output.cpp
Output.h
PPP.cpp
PPP.h
Particle.cpp
Particle.h
ParticleStack.cpp
ParticleStack.h
PropagationEngine.cpp
PropagationEngine.h
ProtonPP.cpp
ProtonPP.h
Randomizer.cpp
Randomizer.h
TableBackgrounds.cpp
TableBackgrounds.h
TableFunction.cpp
TableFunction.h
TableReader.cpp
TableReader.h
Test.cpp
Test.h
Thinning.cpp
Thinning.h
Units.cpp
Units.h
Utils.cpp
Utils.h
main.cpp
PhotoDisintegration.cpp
PhotoDisintegration.h
SteckerEBL.cpp
SteckerEBL.h
Stecker16Background.cpp
Stecker16Background.h
Jet.cpp
Jet.h
Observer.cpp
Observer.h
        Logger.cpp Logger.h)

add_library(MCPropag ${SOURCE_FILES})
add_executable(MCPropagTest ExampleUserMain.cpp)
add_executable(z2t z2t.cpp)
target_link_libraries(MCPropag ${CMAKE_CXX_FLAGS_RELEASE} ${FFTW_OMP_LIBRARY} ${FFTW_LIBRARIES} gfortran ${C_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
target_link_libraries(MCPropagTest MCPropag)
target_link_libraries(z2t MCPropag)
#target_include_directories(MCPropag PRIVATE ${FFTW_INCLUDE_DIRS})
