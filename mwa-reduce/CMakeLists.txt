cmake_minimum_required(VERSION 2.8.8)

project(mwa-reduce)

# Add CMake subdirectory to module search path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

set(CMAKE_CXX_FLAGS "-O3 -ggdb -Wvla -Wall -DNDEBUG -march=native -std=c++0x")

SET(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits scimath scimath_f)

find_library(FFTW3_LIB fftw3 REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)
find_package(Boost COMPONENTS date_time filesystem system thread REQUIRED)
find_library(FITSIO_LIB cfitsio REQUIRED)
find_library(DL_LIB dl REQUIRED)
find_path(FITSIO_INCLUDE_DIR NAMES fitsio.h)
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
find_library(GSL_CBLAS_LIB NAMES gslcblas)
find_path(GTS_INCLUDE_DIR NAMES gts.h)
find_path(GLIB_INCLUDE_DIR NAMES glib.h PATH_SUFFIXES glib-2.0)
find_package(LibXml2)

find_package(HDF5 COMPONENTS CXX REQUIRED)

set(CASA_LIBS ${CASACORE_LIBRARIES})

find_path(CASA_INCLUDE_DIR NAMES ms/MeasurementSets/MeasurementSet.h PATHS /usr/local/include/casacore PATH_SUFFIXES casacore)

include_directories(${CASACORE_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${FITSIO_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})

set(BEAM_LIBS mwareducebeam ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${HDF5_LIBRARIES})

if(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  include_directories(${GSL_INCLUDE_DIR})
endif(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)

if(GTS_INCLUDE_DIR AND GLIB_INCLUDE_DIR)
  include_directories(${GTS_INCLUDE_DIR} ${GLIB_INCLUDE_DIR})
endif(GTS_INCLUDE_DIR AND GLIB_INCLUDE_DIR)

set(CMAKE_REQUIRED_FLAGS "-std=c++0x")
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
"#include \"${CMAKE_CURRENT_SOURCE_DIR}/aocommon/uvector_11.h\"
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

add_library(mwareducebeam fitsiochecker.cpp beam/joneslookupdipole.cpp beam/system.cpp beam/tilebeambase.cpp beam/tilebeam2013.cpp beam/tilebeam2014.cpp beam/tilebeam2016.cpp beam/beam2016implementation.cpp beam/tileimpedance.cpp beam/lnaimpedance.cpp alglib/ap.cpp alglib/alglibinternal.cpp alglib/alglibmisc.cpp alglib/integration.cpp alglib/interpolation.cpp alglib/linalg.cpp alglib/optimization.cpp alglib/solvers.cpp alglib/specialfunctions.cpp)
set(BEAMFILES "")

add_executable(addimg addimg.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp)
add_executable(aegean2model aegean2model.cpp fitsiochecker.cpp fitsreader.cpp model/model.cpp nlplfitter.cpp)
add_executable(apparently apparently.cpp model/model.cpp nlplfitter.cpp fitsiochecker.cpp fitsreader.cpp)
add_executable(applybeam applybeam.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp)
add_executable(applyion applyion.cpp fftconvolver.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp model/model.cpp modelrenderer.cpp nlplfitter.cpp)
add_executable(applysolutions applysolutions.cpp)
add_executable(areas areas.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp)
add_executable(autoprocess autoprocess.cpp calibrator.cpp model/model.cpp beamevaluator.cpp peeler.cpp calibrationmethod.cpp mspredicter.cpp nlplfitter.cpp predicter.cpp progressbar.cpp ${BEAMFILES})
add_executable(beam beam.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp mwa/metafitsfile.cpp mwa/mwaconfig.cpp progressbar.cpp ${BEAMFILES})
add_executable(calibrate calibrate.cpp calibrator.cpp beamevaluator.cpp calibrationmethod.cpp model/model.cpp predicter.cpp mspredicter.cpp nlplfitter.cpp ${BEAMFILES})
add_executable(collectionsolutions collectionsolutions.cpp beamevaluator.cpp model/model.cpp nlplfitter.cpp)
add_executable(combinespectra combinespectra.cpp beamevaluator.cpp imageweights.cpp model/model.cpp mspredicter.cpp multibanddata.cpp nlplfitter.cpp predicter.cpp progressbar.cpp spectrummaker.cpp ${BEAMFILES})
add_executable(cluster cluster.cpp model/model.cpp nlplfitter.cpp)
add_executable(fitsmodel fitsmodel.cpp fitsiochecker.cpp fitsreader.cpp image.cpp nlplfitter.cpp polynomialfitter.cpp deconvolution/deconvolutionalgorithm.cpp deconvolution/spectralfitter.cpp model/model.cpp wsclean/logger.cpp)
add_executable(flagantennae flagantennae.cpp)
add_executable(flagbaselines flagbaselines.cpp)
add_executable(flagionsolutions flagionsolutions.cpp)
add_executable(flagmwa flagmwa.cpp)
add_executable(flagsolutions flagsolutions.cpp)
add_executable(flagsubbands flagsubbands.cpp)
add_executable(ionpeel ionpeel.cpp imageweights.cpp progressbar.cpp beamevaluator.cpp ionpeeler.cpp model/model.cpp multibanddata.cpp nlplfitter.cpp predicter.cpp ${BEAMFILES})
add_executable(ionsol ionsol.cpp model/model.cpp nlplfitter.cpp ${BEAMFILES})
add_executable(ionspectrum ionspectrum.cpp beamevaluator.cpp fluxaccumulator.cpp imageweights.cpp model/model.cpp multibanddata.cpp nlplfitter.cpp ${BEAMFILES})
add_executable(matchsources matchsources.cpp model/model.cpp nlplfitter.cpp progressbar.cpp)
add_executable(model2kvis model2kvis.cpp model/model.cpp nlplfitter.cpp)
add_executable(mrc2model model/model.cpp mrc2model.cpp mrccatalogue.cpp nlplfitter.cpp)
add_executable(mwafinfo mwafinfo.cpp)
add_executable(pbaddimg pbaddimg.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp)
add_executable(pbcorrect pbcorrect.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp)
add_executable(peel peel.cpp peeler.cpp beamevaluator.cpp calibrationmethod.cpp model/model.cpp predicter.cpp mspredicter.cpp nlplfitter.cpp ${BEAMFILES})
add_executable(phasecal beamevaluator.cpp phasecal.cpp model/model.cpp nlplfitter.cpp predicter.cpp ${BEAMFILES})
add_executable(regridimg regridimg.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp)
add_executable(render render.cpp modelrenderer.cpp fftconvolver.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp model/model.cpp nlplfitter.cpp)
add_executable(editmodel editmodel.cpp beamevaluator.cpp fitsreader.cpp fitsiochecker.cpp model/model.cpp nlplfitter.cpp polynomialfitter.cpp rmsynthesis.cpp spectrumft.cpp deconvolution/spectralfitter.cpp ${BEAMFILES})
add_executable(simcompress simcompress.cpp gausencoder.cpp)
add_executable(storetime storetime.cpp stopwatch.cpp)
add_executable(subtrmodel beamevaluator.cpp subtrmodel.cpp model/model.cpp nlplfitter.cpp predicter.cpp ${BEAMFILES} mspredicter.cpp progressbar.cpp)
#add_executable(subtrspectrum subtrspectrum.cpp spectrumsubtractor.cpp beamevaluator.cpp model/model.cpp predicter.cpp ${BEAMFILES})
#add_executable(spectrum beamevaluator.cpp spectrum.cpp spectrummaker.cpp imageweights.cpp model/model.cpp predicter.cpp mspredicter.cpp progressbar.cpp ${BEAMFILES})
add_executable(tests tests.cpp nlplfitter.cpp ${BEAMFILES})

if(LIBXML2_FOUND)
	include_directories(${LIBXML2_INCLUDE_DIR})
	
	add_executable(vo2model vo2model.cpp model/model.cpp nlplfitter.cpp parser/votableparser.cpp)
	target_link_libraries(vo2model ${LIBXML2_LIBRARIES})

	add_executable(scaleimage scaleimage.cpp fitsiochecker.cpp fitsreader.cpp fitswriter.cpp parser/votableparser.cpp)
	target_link_libraries(scaleimage ${LIBXML2_LIBRARIES} ${FITSIO_LIB} ${CASA_LIBS})
	message("Libxml2 found: will build vo2model and scaleimage")
else(LIBXML2_FOUND)
  message("Libxml2 NOT found: will NOT build vo2model and scaleimage")
endif(LIBXML2_FOUND)

target_link_libraries(addimg ${FITSIO_LIB} ${CASA_LIBS} ${Boost_DATE_TIME_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(aegean2model ${FITSIO_LIB} ${CASA_LIBS} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(applybeam ${FITSIO_LIB} ${CASA_LIBS})
target_link_libraries(apparently ${FITSIO_LIB} ${CASA_LIBS} ${GSL_LIB} ${GSL_CBLAS_LIB} )
target_link_libraries(areas ${FITSIO_LIB} ${CASA_LIBS})

target_link_libraries(applyion ${FITSIO_LIB} ${CASA_LIBS} ${FFTW3_LIB} ${GSL_LIB} ${GSL_CBLAS_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(applysolutions ${CASA_LIBS})
target_link_libraries(autoprocess ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${CASA_MEASURES_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(beam ${BEAM_LIBS} ${CASA_LIBS} ${CASA_MEASURES_LIB} ${FITSIO_LIB} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(calibrate ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(collectionsolutions ${BEAM_LIBS} ${FITSIO_LIB} ${CASA_LIBS} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(combinespectra ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(editmodel ${BEAM_LIBS} ${FITSIO_LIB} ${CASA_LIBS} ${FFTW3_LIB})
target_link_libraries(fitsmodel ${FITSIO_LIB} ${CASA_LIBS} ${GSL_LIB} ${GSL_CBLAS_LIB} ${Boost_DATE_TIME_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(flagantennae ${CASA_LIBS})
target_link_libraries(flagbaselines ${CASA_LIBS})
target_link_libraries(flagmwa ${CASA_LIBS})
target_link_libraries(flagsubbands ${CASA_LIBS})
target_link_libraries(flagsolutions ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(ionpeel ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(ionspectrum ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(mrc2model ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(pbaddimg ${FITSIO_LIB} ${CASA_LIBS} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(pbcorrect ${FITSIO_LIB} ${CASA_LIBS} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(peel ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(phasecal ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(regridimg ${CASA_LIBS} ${FITSIO_LIB})
target_link_libraries(render ${FITSIO_LIB} ${CASA_LIBS} ${FFTW3_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(simcompress ${CASA_LIBS})
#target_link_libraries(spectrum ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(subtrmodel ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
#target_link_libraries(subtrspectrum ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(tests ${BEAM_LIBS} ${CASA_LIBS} ${FITSIO_LIB})

add_executable(beam2016test EXCLUDE_FROM_ALL beam/beam2016implementation.cpp beam/beam2016test.cpp beam/system.cpp )
target_link_libraries(beam2016test ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${HDF5_LIBRARIES})

if(NOT GSL_CBLAS_LIB)
  message("GSL CBLAS lib was not found. GSL needs CBLAS: disabling GSL.")
endif(NOT GSL_CBLAS_LIB)

if(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  add_definitions(-DHAVE_GSL)
  target_link_libraries(cluster ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(ionpeel ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(ionsol ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(tests ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(editmodel ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(simcompress ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(model2kvis ${GSL_LIB} ${GSL_CBLAS_LIB} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
  target_link_libraries(matchsources ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(render ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(vo2model ${GSL_LIB} ${GSL_CBLAS_LIB})
  
	add_executable(sedcombine sedcombine.cpp sedanalyser.cpp model/model.cpp nlplfitter.cpp rmsynthesis.cpp progressbar.cpp spectrumft.cpp subbandpassband.cpp)
  target_link_libraries(sedcombine ${FFTW3_LIB} ${GSL_LIB} ${GSL_CBLAS_LIB})
  
else(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)
  message("GNU scientific library (gsl) not found: some features will be disabled (this is not an error).")
endif(GSL_LIB AND GSL_CBLAS_LIB AND GSL_INCLUDE_DIR)

if(UVECTOR_COMPILES)
  add_definitions(-DHAVE_UVECTOR)
else(UVECTOR_COMPILES)
  message("ao::uvector does not compile: will use slower std::vector (this is not an error).")
endif(UVECTOR_COMPILES)
if(LANE11_COMPILES)
  add_definitions(-DHAVE_LANE11)
else(LANE11_COMPILES)
  message("ao::lane does not compile with C++11 features: will only use C++03 (this is not an error).")
endif(LANE11_COMPILES)

