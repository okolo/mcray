cmake_minimum_required (VERSION 2.6)
project (savecoef)
# uncomment the line below for verbose output
#set(CMAKE_VERBOSE_MAKEFILE ON)
add_definitions(-DUSE_GSL)

include_directories ("../../lib" "../../OS/include")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_DEBUG -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O2")
set(CMAKE_PREFIX_PATH "external")

# -rdynamic flux is only supported on systems with ELF executable format
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -rdynamic")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -rdynamic")
endif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")

FIND_LIBRARY(C_LIBRARY c)
FIND_LIBRARY(GSL_LIBRARY gsl)
FIND_LIBRARY(GSLCBLAS_LIBRARY gslcblas)
add_subdirectory (../../lib mcray)

set(SOURCE_FILES
savecoef.cpp savecoef.h
)
add_link_options(-fopenmp)
add_executable(savecoef ${SOURCE_FILES})
target_link_libraries(savecoef mcray ${C_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
