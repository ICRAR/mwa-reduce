cmake_minimum_required(VERSION 2.6)

project(wsclean)

set(CMAKE_CXX_FLAGS "-O3 -Wall -DNDEBUG -march=native -msse4.1")

find_library(CASA_MS_LIB casa_ms REQUIRED)
find_library(CASA_TABLES_LIB casa_tables REQUIRED)
find_library(CASA_CASA_LIB casa_casa REQUIRED)
find_library(CASA_MEASURES_LIB casa_measures REQUIRED)
find_library(CASA_FITS_LIB casa_fits REQUIRED)
find_library(FFTW3_LIB fftw3 REQUIRED)
find_package(Boost COMPONENTS thread system REQUIRED)
find_library(FITSIO_LIB cfitsio REQUIRED)
find_library(DL_LIB dl REQUIRED)
find_path(FITSIO_INCLUDE_DIR NAMES fitsio.h)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)

set(CASA_LIBS ${CASA_MS_LIB} ${CASA_TABLES_LIB} ${CASA_MEASURES_LIB} ${CASA_CASA_LIB} ${CASA_FITS_LIB} ${DL_LIB})

find_path(CASA_INCLUDE_DIR NAMES ms/MeasurementSets/MeasurementSet.h PATHS /usr/local/include/casacore PATH_SUFFIXES casacore)

include_directories(${CASA_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${FITSIO_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(CMAKE_REQUIRED_FLAGS "-std=c++0x")
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
"#include \"${CMAKE_CURRENT_SOURCE_DIR}/aocommon/uvector.h\"
int main(int argc, char* argv[])
{
return 0;
}" UVECTOR_COMPILES)
check_cxx_source_compiles(
"#include \"${CMAKE_CURRENT_SOURCE_DIR}/aocommon/lane_11.h\"
int main(int argc, char* argv[])
{
return 0;
}" LANE11_COMPILES)

if(UVECTOR_COMPILES)
  add_definitions(-DHAVE_UVECTOR)
endif(UVECTOR_COMPILES)
if(LANE11_COMPILES)
  add_definitions(-DHAVE_LANE11)
endif(LANE11_COMPILES)

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

add_executable(wsclean wscleanmain.cpp wsclean.cpp beamevaluator.cpp fftresampler.cpp fitsreader.cpp fitswriter.cpp imageweights.cpp layeredimager.cpp model.cpp modelrenderer.cpp progressbar.cpp stopwatch.cpp tilebeam.cpp wsinversion.cpp cleanalgorithms/cleanalgorithm.cpp cleanalgorithms/joinedpolclean.cpp cleanalgorithms/simpleclean.cpp msprovider/contiguousms.cpp msprovider/msprovider.cpp msprovider/partitionedms.cpp)

set_target_properties(wsclean PROPERTIES COMPILE_FLAGS "-std=c++0x")

target_link_libraries(wsclean ${CASA_LIBS} ${FFTW3_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${FITSIO_LIB})

install(TARGETS wsclean DESTINATION bin)
