cmake_minimum_required (VERSION 2.6)
project (mcray)
# uncomment the line below for verbose output
#set(CMAKE_VERBOSE_MAKEFILE ON)
enable_language (Fortran)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS "-ffixed-line-length-200 -fno-f2c -fopenmp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS "-f77rtl -fopenmp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS "-fopenmp -fno-f2c")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2  -fopenmp...")
  set (CMAKE_Fortran_FLAGS "-fopenmp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


add_definitions(-DUSE_GSL)
include_directories(../OS/include ../external)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG -O0 -pg")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O2")
set(CMAKE_PREFIX_PATH "../external")

# -rdynamic flux is only supported on systems with ELF executable format
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -rdynamic")

endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

FIND_LIBRARY(C_LIBRARY c)
FIND_LIBRARY(GSL_LIBRARY gsl)
FIND_LIBRARY(GSLCBLAS_LIBRARY gslcblas)

set(SOURCE_FILES
../external/SOPHIA/eventgen.f
        ../external/SOPHIA/inpoutput.f
        ../external/SOPHIA/jetset74dp.f
        ../external/SOPHIA/rndm.f
        ../external/SOPHIA/sampling.f
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
Logger.cpp
Logger.h)

add_link_options(-fopenmp)
add_library(mcray ${SOURCE_FILES})
add_executable(mcray_test ExampleUserMain.cpp)
add_executable(z2t z2t.cpp)

target_link_libraries(mcray gfortran ${C_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
target_link_libraries(mcray_test mcray)
target_link_libraries(z2t mcray)
