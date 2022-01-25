cmake_minimum_required(VERSION 2.8.8)

project(mwa-reduce)

# enable CUDA version of the MWA-2016 beam model
if(CUDA_SUPPORT)
  enable_language(CUDA)
endif(CUDA_SUPPORT)

if(PORTABLE)
  add_compile_options(-std=c++11 -O3 -ggdb -Wvla -Wall -DNDEBUG)
else()
  add_compile_options(
    -std=c++11
    -O3
    -ggdb
    -Wvla
    -Wall
    -Wextra
    -DNDEBUG
    -march=native)
endif(PORTABLE)

# Find and include git submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
              --depth 1
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
          "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules"
      )
    endif()
  endif()
endif()

# Include aocommon
include_directories("${CMAKE_SOURCE_DIR}/external/aocommon/include/")

# Add CMake subdirectory to module search path
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

find_package(
  HDF5
  COMPONENTS C CXX
  REQUIRED)

set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(
  Casacore REQUIRED
  COMPONENTS casa
             ms
             tables
             measures
             fits)
set(CASA_LIBS ${CASACORE_LIBRARIES})

find_package(CFITSIO REQUIRED)

find_package(Threads REQUIRED)
find_library(FFTW3_LIB fftw3 REQUIRED)
find_library(FFTW3_THREADS_LIB fftw3_threads REQUIRED)
find_library(FFTW3F_LIB fftw3f REQUIRED)
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h)
find_package(
  Boost
  COMPONENTS date_time filesystem system thread
  REQUIRED)
find_library(DL_LIB dl REQUIRED)
find_library(GSL_LIB NAMES gsl)
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_version.h)
find_library(GSL_CBLAS_LIB NAMES gslcblas)
find_path(GTS_INCLUDE_DIR NAMES gts.h)
find_path(
  GLIB_INCLUDE_DIR
  NAMES glib.h
  PATH_SUFFIXES glib-2.0)
find_package(LibXml2)

find_path(
  CASA_INCLUDE_DIR
  NAMES ms/MeasurementSets/MeasurementSet.h
  PATHS /usr/local/include/casacore
  PATH_SUFFIXES casacore)

include_directories(${CASACORE_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${CFITSIO_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIRS})

set(BEAM_LIBS mwareducebeam ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY}
              ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})

if(GSL_LIB
   AND GSL_CBLAS_LIB
   AND GSL_INCLUDE_DIR)
  include_directories(${GSL_INCLUDE_DIR})
endif(
  GSL_LIB
  AND GSL_CBLAS_LIB
  AND GSL_INCLUDE_DIR)

if(GTS_INCLUDE_DIR AND GLIB_INCLUDE_DIR)
  include_directories(${GTS_INCLUDE_DIR} ${GLIB_INCLUDE_DIR})
endif(GTS_INCLUDE_DIR AND GLIB_INCLUDE_DIR)

# The following stuff will set the "rpath" correctly, so that LD_LIBRARY_PATH
# doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# when building, don't use the install RPATH already (but later on when
# installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH which point to directories
# outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

link_libraries(${CMAKE_THREAD_LIBS_INIT})

add_library(
  mwareducebeam
  beam/joneslookupdipole.cpp
  beam/system.cpp
  beam/tilebeambase.cpp
  beam/tilebeam2013.cpp
  beam/tilebeam2014.cpp
  beam/tilebeam2016.cpp
  beam/beam2016implementation.cpp
  beam/tileimpedance.cpp
  beam/lnaimpedance.cpp
  alglib/ap.cpp
  alglib/alglibinternal.cpp
  alglib/alglibmisc.cpp
  alglib/integration.cpp
  alglib/interpolation.cpp
  alglib/linalg.cpp
  alglib/optimization.cpp
  alglib/solvers.cpp
  alglib/specialfunctions.cpp)
set(BEAMFILES "")

add_executable(addimg addimg.cpp)
add_executable(aegean2model aegean2model.cpp model/model.cpp
                            math/nlplfitter.cpp)
add_executable(apparently apparently.cpp model/model.cpp math/nlplfitter.cpp)
add_executable(applybeam applybeam.cpp)
add_executable(
  applyion
  applyion.cpp
  fftconvolver.cpp
  fftkernels.cpp
  fftwmanager.cpp
  model/model.cpp
  math/modelrenderer.cpp
  math/nlplfitter.cpp)
add_executable(applysolutions applysolutions.cpp)
add_executable(
  autoprocess
  autoprocess.cpp
  calibrator.cpp
  model/model.cpp
  beamevaluator.cpp
  peeler.cpp
  calibrationmethod.cpp
  mspredicter.cpp
  math/nlplfitter.cpp
  predicter.cpp
  progressbar.cpp
  ${BEAMFILES})
add_executable(bbs2model bbs2model.cpp model/model.cpp math/nlplfitter.cpp)
add_executable(beam beam.cpp mwa/metafitsfile.cpp mwa/mwaconfig.cpp
                    progressbar.cpp ${BEAMFILES})
add_executable(
  calibrate
  calibrate.cpp
  calibrator.cpp
  beamevaluator.cpp
  calibrationmethod.cpp
  fftwmanager.cpp
  model/model.cpp
  predicter.cpp
  mspredicter.cpp
  math/nlplfitter.cpp
  progressbar.cpp
  ${BEAMFILES})
add_executable(cluster cluster.cpp model/model.cpp math/nlplfitter.cpp)
add_executable(
  fitsmodel
  fitsmodel.cpp
  image.cpp
  math/nlplfitter.cpp
  math/polynomialfitter.cpp
  deconvolution/spectralfitter.cpp
  model/model.cpp
  wsclean/logger.cpp)
add_executable(flagantennae flagantennae.cpp)
add_executable(flagbaselines flagbaselines.cpp)
add_executable(flagionsolutions flagionsolutions.cpp)
add_executable(flagmwa flagmwa.cpp)
add_executable(flagsolutions flagsolutions.cpp)
add_executable(flagsubbands flagsubbands.cpp)
add_executable(
  ionpeel
  ionpeel.cpp
  imageweights.cpp
  progressbar.cpp
  beamevaluator.cpp
  ionpeeler.cpp
  model/model.cpp
  math/nlplfitter.cpp
  predicter.cpp
  wsclean/logger.cpp
  ${BEAMFILES})
add_executable(matchsources matchsources.cpp model/model.cpp
                            math/nlplfitter.cpp progressbar.cpp)
add_executable(mrc2model model/model.cpp mrc2model.cpp mrccatalogue.cpp
                         math/nlplfitter.cpp)
add_executable(mwafinfo mwafinfo.cpp)
add_executable(pbaddimg pbaddimg.cpp)
add_executable(pbcorrect pbcorrect.cpp)
add_executable(
  peel
  peel.cpp
  peeler.cpp
  beamevaluator.cpp
  calibrationmethod.cpp
  model/model.cpp
  predicter.cpp
  mspredicter.cpp
  math/nlplfitter.cpp
  ${BEAMFILES})
add_executable(phasecal beamevaluator.cpp phasecal.cpp model/model.cpp
                        math/nlplfitter.cpp predicter.cpp ${BEAMFILES})
add_executable(regridimg regridimg.cpp)
add_executable(
  render
  fftwmanager.cpp
  render.cpp
  math/modelrenderer.cpp
  fftconvolver.cpp
  fftkernels.cpp
  image.cpp
  wsclean/logger.cpp
  model/model.cpp
  math/nlplfitter.cpp
  render/interpolatingrenderer.cpp)
add_executable(
  editmodel
  editmodel.cpp
  beamevaluator.cpp
  fftwmanager.cpp
  image.cpp
  model/model.cpp
  math/nlplfitter.cpp
  math/polynomialfitter.cpp
  rmsynthesis.cpp
  spectrumft.cpp
  deconvolution/spectralfitter.cpp
  ${BEAMFILES})
add_executable(storetime storetime.cpp stopwatch.cpp)
add_executable(
  subtrmodel
  beamevaluator.cpp
  fftwmanager.cpp
  subtrmodel.cpp
  model/model.cpp
  math/nlplfitter.cpp
  predicter.cpp
  ${BEAMFILES}
  mspredicter.cpp
  progressbar.cpp)
add_executable(solutiontool solutiontool.cpp)

if(CUDA_SUPPORT)
  target_compile_definitions(calibrate PUBLIC -DCUDA_SUPPORT)
  target_compile_definitions(mwareducebeam PUBLIC -DCUDA_SUPPORT)
  set_property(TARGET calibrate PROPERTY CUDA_SEPARABLE_COMPILATION ON)
  add_library(jones_cuda beam/gpu/factorial.c beam/gpu/jones_gpu.cu
                         beam/gpu/legendre.c)
  target_compile_definitions(jones_cuda PUBLIC -DCUDA_SUPPORT)
  set_target_properties(jones_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_link_libraries(mwareducebeam jones_cuda)
  set_target_properties(mwareducebeam PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_target_properties(beam PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif(CUDA_SUPPORT)

# Experiments
add_executable(
  ionspectrum EXCLUDE_FROM_ALL
  experiments/ionspectrum/ionspectrum.cpp
  beamevaluator.cpp
  fluxaccumulator.cpp
  imageweights.cpp
  model/model.cpp
  math/nlplfitter.cpp
  ${BEAMFILES})

add_executable(
  combinespectra EXCLUDE_FROM_ALL
  experiments/spectrummaker/combinespectra.cpp
  beamevaluator.cpp
  imageweights.cpp
  model/model.cpp
  mspredicter.cpp
  math/nlplfitter.cpp
  predicter.cpp
  progressbar.cpp
  experiments/spectrummaker/spectrummaker.cpp
  ${BEAMFILES})

add_executable(
  ionsol EXCLUDE_FROM_ALL experiments/ionspectrum/ionsol.cpp model/model.cpp
                          math/nlplfitter.cpp ${BEAMFILES})

add_executable(
  collectionsolutions EXCLUDE_FROM_ALL
  experiments/ionspectrum/collectionsolutions.cpp beamevaluator.cpp
  model/model.cpp math/nlplfitter.cpp)

if(LIBXML2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})

  add_executable(vo2model vo2model.cpp model/model.cpp math/nlplfitter.cpp
                          parser/votableparser.cpp)
  target_link_libraries(vo2model ${LIBXML2_LIBRARIES})

  add_executable(scaleimage scaleimage.cpp parser/votableparser.cpp)
  target_link_libraries(scaleimage ${LIBXML2_LIBRARIES} ${CFITSIO_LIBRARY}
                        ${CASA_LIBS})
  message("Libxml2 found: will build vo2model and scaleimage")
else(LIBXML2_FOUND)
  message("Libxml2 NOT found: will NOT build vo2model and scaleimage")
endif(LIBXML2_FOUND)

target_link_libraries(addimg ${CFITSIO_LIBRARY} ${CASA_LIBS}
                      ${Boost_DATE_TIME_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(aegean2model ${CFITSIO_LIBRARY} ${CASA_LIBS} ${GSL_LIB}
                      ${GSL_CBLAS_LIB})
target_link_libraries(applybeam ${CFITSIO_LIBRARY} ${CASA_LIBS})
target_link_libraries(apparently ${CFITSIO_LIBRARY} ${CASA_LIBS} ${GSL_LIB}
                      ${GSL_CBLAS_LIB})

target_link_libraries(
  applyion
  ${CFITSIO_LIBRARY}
  ${CASA_LIBS}
  ${FFTW3_LIB}
  ${FFTW3F_LIB}
  ${FFTW3_THREADS_LIB}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB}
  ${Boost_SYSTEM_LIBRARY})
target_link_libraries(applysolutions ${CASA_LIBS})
target_link_libraries(
  autoprocess
  ${BEAM_LIBS}
  ${CASA_LIBS}
  ${CFITSIO_LIBRARY}
  ${CASA_MEASURES_LIB}
  ${Boost_SYSTEM_LIBRARY}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB})
target_link_libraries(bbs2model ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(
  beam
  ${BEAM_LIBS}
  ${CASA_LIBS}
  ${CASA_MEASURES_LIB}
  ${CFITSIO_LIBRARY}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB})
target_link_libraries(
  calibrate
  ${BEAM_LIBS}
  ${CASA_LIBS}
  ${FFTW3_LIB}
  ${FFTW3_THREADS_LIB}
  ${CFITSIO_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB})
target_link_libraries(collectionsolutions ${BEAM_LIBS} ${CFITSIO_LIBRARY}
                      ${CASA_LIBS} ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(
  combinespectra
  ${BEAM_LIBS}
  ${CASA_LIBS}
  ${CFITSIO_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB})
target_link_libraries(
  editmodel
  ${BEAM_LIBS}
  ${FFTW3_LIB}
  ${FFTW3_THREADS_LIB}
  ${CFITSIO_LIBRARY}
  ${CASA_LIBS}
  ${FFTW3_LIB})
target_link_libraries(
  fitsmodel
  ${CFITSIO_LIBRARY}
  ${CASA_LIBS}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB}
  ${Boost_DATE_TIME_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY})
target_link_libraries(flagantennae ${CASA_LIBS})
target_link_libraries(flagbaselines ${CASA_LIBS})
target_link_libraries(flagmwa ${CASA_LIBS})
target_link_libraries(flagsubbands ${CASA_LIBS})
target_link_libraries(flagsolutions ${Boost_SYSTEM_LIBRARY}
                      ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(
  ionpeel
  ${BEAM_LIBS}
  ${CASA_LIBS}
  ${CFITSIO_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB}
  ${Boost_DATE_TIME_LIBRARY})
target_link_libraries(
  ionspectrum
  ${BEAM_LIBS}
  ${CASA_LIBS}
  ${CFITSIO_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB})
target_link_libraries(mrc2model ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(pbaddimg ${CFITSIO_LIBRARY} ${CASA_LIBS}
                      ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(pbcorrect ${CFITSIO_LIBRARY} ${CASA_LIBS}
                      ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(
  peel
  ${BEAM_LIBS}
  ${CASA_LIBS}
  ${CFITSIO_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB})
target_link_libraries(phasecal ${BEAM_LIBS} ${CASA_LIBS} ${CFITSIO_LIBRARY}
                      ${GSL_LIB} ${GSL_CBLAS_LIB})
target_link_libraries(regridimg ${CASA_LIBS} ${CFITSIO_LIBRARY})
target_link_libraries(
  render
  ${CFITSIO_LIBRARY}
  ${CASA_LIBS}
  ${FFTW3_LIB}
  ${FFTW3F_LIB}
  ${FFTW3_THREADS_LIB}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY})
target_link_libraries(
  subtrmodel
  ${BEAM_LIBS}
  ${CASA_LIBS}
  ${FFTW3_LIB}
  ${FFTW3_THREADS_LIB}
  ${CFITSIO_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${GSL_LIB}
  ${GSL_CBLAS_LIB})

add_executable(
  beam2016test EXCLUDE_FROM_ALL beam/beam2016implementation.cpp
                                beam/beam2016test.cpp beam/system.cpp)
target_link_libraries(beam2016test ${Boost_SYSTEM_LIBRARY}
                      ${Boost_FILESYSTEM_LIBRARY} ${HDF5_LIBRARIES})

if(NOT GSL_CBLAS_LIB)
  message("GSL CBLAS lib was not found. GSL needs CBLAS: disabling GSL.")
endif(NOT GSL_CBLAS_LIB)

if(GSL_LIB
   AND GSL_CBLAS_LIB
   AND GSL_INCLUDE_DIR)
  add_definitions(-DHAVE_GSL)
  target_link_libraries(cluster ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(ionpeel ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(ionsol ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(editmodel ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(matchsources ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(render ${GSL_LIB} ${GSL_CBLAS_LIB})
  target_link_libraries(vo2model ${GSL_LIB} ${GSL_CBLAS_LIB})

  add_executable(
    sedcombine
    experiments/sed/sedcombine.cpp
    experiments/sed/sedanalyser.cpp
    model/model.cpp
    math/nlplfitter.cpp
    image.cpp
    rmsynthesis.cpp
    progressbar.cpp
    spectrumft.cpp
    subbandpassband.cpp)
  target_link_libraries(
    sedcombine
    ${CASA_LIBS}
    ${FFTW3_LIB}
    ${GSL_LIB}
    ${GSL_CBLAS_LIB}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY})

else(
  GSL_LIB
  AND GSL_CBLAS_LIB
  AND GSL_INCLUDE_DIR)
  message(
    "GNU scientific library (gsl) not found: some features will be disabled (this is not an error)."
  )
endif(
  GSL_LIB
  AND GSL_CBLAS_LIB
  AND GSL_INCLUDE_DIR)

message(STATUS "Flags passed to C++ compiler: " ${CMAKE_CXX_FLAGS})
