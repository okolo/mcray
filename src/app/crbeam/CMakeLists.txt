cmake_minimum_required (VERSION 3.13)
project (CRbeam)
# uncomment the line below for verbose output
#set(CMAKE_VERBOSE_MAKEFILE ON)
add_definitions(-DUSE_GSL)

include_directories ("../../lib" "../../OS/include" "../../../OS/include")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG -O0 -pg")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O2")

# For executable targets
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
# For shared library targets
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
# For module targets
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")

#set(CMAKE_PREFIX_PATH "external")

# -rdynamic flux is only supported on systems with ELF executable format
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -rdynamic")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# FIND_LIBRARY(C_LIBRARY c)
FIND_LIBRARY(GSL_LIBRARY gsl)
FIND_LIBRARY(GSLCBLAS_LIBRARY gslcblas)
FIND_LIBRARY(BOOST_FS_LIBRARY boost_filesystem)
FIND_LIBRARY(BOOST_SYS_LIBRARY boost_system)
add_subdirectory (../../lib mcray)

set(SOURCE_FILES
CmdLine.cpp
CmdLine.h
CRbeam.cpp
CRbeam.h
)
add_link_options(-fopenmp)
add_executable(CRbeam ${SOURCE_FILES})
target_link_libraries(CRbeam mcray ${BOOST_SYS_LIBRARY} ${BOOST_FS_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY} ${CMAKE_DL_LIBS})
