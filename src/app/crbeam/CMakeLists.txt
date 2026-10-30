cmake_minimum_required (VERSION 3.16)
project (CRbeam)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(CTest)
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
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
    # For shared library targets
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
    # For module targets
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
endif()

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

# YAML is an application input dependency, not a physics-library dependency.
find_package(yaml-cpp 0.6 CONFIG REQUIRED)
find_package(OpenSSL 1.1.1 REQUIRED COMPONENTS Crypto)
if(TARGET yaml-cpp::yaml-cpp)
    set(CRBEAM_YAML_TARGET yaml-cpp::yaml-cpp)
else()
    set(CRBEAM_YAML_TARGET yaml-cpp)
endif()
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE CRBEAM_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT CRBEAM_REVISION)
    set(CRBEAM_REVISION "unknown")
endif()
execute_process(COMMAND git diff HEAD WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE CRBEAM_PATCH ERROR_QUIET)
string(SHA256 CRBEAM_PATCH_HASH "${CRBEAM_PATCH}")
set(CRBEAM_APPLICATION_SOURCES Configuration.cpp Configuration.h MagneticConfiguration.cpp MagneticConfiguration.h SpatialConfiguration.cpp SpatialConfiguration.h Sha256.cpp Sha256.h CRbeam.cpp CRbeam.h CmdLine.cpp CmdLine.h CMakeLists.txt)
set(CRBEAM_APPLICATION_HASHES "")
foreach(source IN LISTS CRBEAM_APPLICATION_SOURCES)
    file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/${source}" source_hash)
    string(APPEND CRBEAM_APPLICATION_HASHES "${source}:${source_hash}\n")
endforeach()
string(SHA256 CRBEAM_CONFIG_HASH "${CRBEAM_APPLICATION_HASHES}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CRBEAM_APPLICATION_SOURCES})
add_library(crbeam_configuration STATIC Configuration.cpp MagneticConfiguration.cpp SpatialConfiguration.cpp Sha256.cpp)
target_include_directories(crbeam_configuration PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "../../lib")
target_link_libraries(crbeam_configuration PUBLIC ${CRBEAM_YAML_TARGET} OpenSSL::Crypto mcray)
target_compile_definitions(crbeam_configuration PRIVATE
    MCRAY_REVISION="${CRBEAM_REVISION}" MCRAY_PATCH_SHA256="${CRBEAM_PATCH_HASH}"
    MCRAY_CONFIG_SHA256="${CRBEAM_CONFIG_HASH}" MCRAY_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
    MCRAY_COMPILER="${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

set(SOURCE_FILES
CmdLine.cpp
CmdLine.h
CRbeam.cpp
CRbeam.h
)
add_link_options(-fopenmp)
add_executable(CRbeam ${SOURCE_FILES})
target_link_libraries(CRbeam mcray crbeam_configuration ${BOOST_SYS_LIBRARY} ${BOOST_FS_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY} ${CMAKE_DL_LIBS})

if(BUILD_TESTING)
    add_subdirectory(../../../tests tests)
endif()
