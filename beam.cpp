#include "beam/tilebeam.h"
#ifndef CUDA_SUPPORT
#include "beam/tilebeam2013.h"
#include "beam/tilebeam2014.h"
#endif
#include "beam/tilebeam2016.h"

// 202008 : GPU version of beam2016 model added (developed by Maciej Cytowski
// and Di Pietrantonio from PAWSEY)

#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/banddata.h>
#include <aocommon/matrix2x2.h>
#include "numberlist.h"
#include "progressbar.h"

#include <aocommon/imagecoordinates.h>

#include "mwa/metafitsfile.h"
#include "mwa/mwaconfig.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ScalarColumn.h>

#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>

#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>

#include <stdexcept>

template <typename TileImplementation, bool DoSquare>
void MakeBeam(double** imgPtr, size_t width, size_t height, double pixelSizeX,
              double pixelSizeY, double refRA, double refDec,
              double arrLatitude, double zenithHa, double zenithDec,
              double centralFrequency, const double* delays,
              const casacore::MeasFrame& frame, double dl, double dm,
              bool freqInterpolation) {
  const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC, frame);
  const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO,
                                             frame);
  const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
  casacore::MDirection::Convert j2000ToHaDecRef(j2000Ref, hadecRef),
      j2000ToAzelGeoRef(j2000Ref, azelgeoRef);

  TileBeamBase<TileImplementation> tilebeam(delays, freqInterpolation, "");
  ProgressBar progressBar("Constructing beam");
  for (size_t y = 0; y != height; ++y) {
    for (size_t x = 0; x != width; ++x) {
      double l, m, ra, dec;
      aocommon::ImageCoordinates::XYToLM(x, y, pixelSizeX, pixelSizeY, width,
                                         height, l, m);
      l += dl;
      m += dm;
      aocommon::ImageCoordinates::LMToRaDec(l, m, refRA, refDec, ra, dec);

      std::complex<double> gain[4];
      tilebeam.ArrayResponse(ra, dec, j2000Ref, j2000ToHaDecRef,
                             j2000ToAzelGeoRef, arrLatitude, zenithHa,
                             zenithDec, centralFrequency, gain);
      if (DoSquare) {
        std::complex<double> gainSq[4];
        aocommon::Matrix2x2::ATimesHermB(gainSq, gain, gain);
        aocommon::Matrix2x2::Assign(gain, gainSq);
      }

      for (size_t i = 0; i != 4; ++i) {
        *imgPtr[i * 2] = gain[i].real();
        *imgPtr[i * 2 + 1] = gain[i].imag();
        ++imgPtr[i * 2];
        ++imgPtr[i * 2 + 1];
      }
    }
    progressBar.SetProgress(y, height);
  }
  progressBar.SetProgress(height, height);
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Syntax: beam [-2013 | -2014 | -2014i | -2016] [-allsky] "
                 "[-square] [-proto <input fitsfile>] [-name <prefix>] [-ms "
                 "<measurementset>] -[m <metafits>] [-delays <0,0,..>]\n"
                 "  -2013 selects the 'analytical' beam\n"
                 "  -2014 selected the beam that takes tile impedance etc. "
                 "into account\n"
                 "  -2014i selects the same beam but also enables frequency "
                 "interpolation\n"
                 "  -2016 selects the improved beam described by M. Sokolowski "
                 "(2017)\n";
    return 0;
  }

  int argi = 1;
  const char* inpFitsname = 0;
  const char* msName = 0;
  const char* metafitsName = 0;
  string prefixName = "beam";
  double delays[16];
  bool doSquare = false, hasDelays = false;
  bool use2013 = false;
  bool use2016 = false;
  bool doInterpolate = false;
  for (size_t i = 0; i != 16; ++i) delays[i] = 0.0;
  while (argi < argc && argv[argi][0] == '-') {
    std::string param(&argv[argi][1]);
    if (param == "proto") {
      ++argi;
      inpFitsname = argv[argi];
    } else if (param == "ms") {
      ++argi;
      msName = argv[argi];
    } else if (param == "m") {
      ++argi;
      metafitsName = argv[argi];
    } else if (param == "name") {
      ++argi;
      prefixName = argv[argi];
    } else if (param == "2013")
      use2013 = true;
    else if (param == "2014") {
      use2013 = false;
      doInterpolate = false;
    } else if (param == "2014i") {
      use2013 = false;
      doInterpolate = true;
    } else if (param == "2016") {
      use2016 = true;
      doInterpolate = false;
    } else if (param == "allsky") {
    } else if (param == "delays") {
      ++argi;
      aocommon::UVector<int> list = NumberList::ParseIntList(argv[argi]);
      if (list.size() != 16) {
        std::cerr << "Need 16 delays\n";
        exit(1);
      }
      for (size_t i = 0; i != 16; ++i) delays[i] = list[i];
      hasDelays = true;
    } else if (param == "square") {
      doSquare = true;
    } else
      throw std::runtime_error(std::string("Invalid param: ") + param);
    ++argi;
  }

  casacore::MPosition arrayPos;
  casacore::MEpoch time;
  double centralFrequency;
  bool haveTime = false;

  /** If nothing else is read in, use some sensible values:*/
  arrayPos = casacore::MPosition(casacore::MVPosition(
      -2.55952e+06, 5.09585e+06, -2.84899e+06));  // pos of tile 011
  centralFrequency = 150000000.0;

  if (metafitsName != 0) {
    MetaFitsFile mfFile(metafitsName);
    MWAHeader header;
    MWAHeaderExt headerExt;
    mfFile.ReadHeader(header, headerExt);
    header.Validate(false);
    std::cout << "Start=" << header.GetDateFirstScanFromFields()
              << ", end=" << header.GetDateLastScanMJD() << '\n';
    time = casacore::MEpoch(casacore::MVEpoch(
        (header.GetDateFirstScanFromFields() + header.GetDateLastScanMJD()) *
        0.5));
    centralFrequency = header.centralFrequencyMHz * 1e6;
    haveTime = true;
    if (!hasDelays) {
      for (size_t i = 0; i != 16; ++i) delays[i] = headerExt.delays[i];
      hasDelays = true;
    }
    std::cout << "Metafits specifies mid time = " << time << '\n';
    std::cout << "Frequency = " << centralFrequency << '\n';
    std::cout << " ***\n"
                 " *** WARNING: using time from metafits file, but these times "
                 "can be off by a few seconds!\n"
                 " ***\n";
  }

  if (msName != 0) {
    /**
     * Read some meta data from the measurement set
     */
    casacore::MeasurementSet ms(msName);
    if (ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");

    casacore::MSAntenna aTable = ms.antenna();
    if (aTable.nrow() == 0) throw std::runtime_error("No antennae in set");
    casacore::MPosition::ROScalarColumn antPosColumn(
        aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
    arrayPos = antPosColumn(0);
    casacore::MEpoch::ROScalarColumn timeColumn(
        ms, ms.columnName(casacore::MSMainEnums::TIME));
    // Find the mid time step
    casacore::MEpoch firstTime = timeColumn(0);
    casacore::MEpoch lastTime = timeColumn(ms.nrow() - 1);
    std::cout << "Start time = " << firstTime << ", end time = " << lastTime
              << '\n';
    time =
        casacore::MEpoch(casacore::MVEpoch(0.5 * (firstTime.getValue().get() +
                                                  lastTime.getValue().get())),
                         firstTime.getRef());
    std::cout << "Using mid time = " << time << '\n';
    haveTime = true;

    casacore::Table mwaTilePointing =
        ms.keywordSet().asTable("MWA_TILE_POINTING");
    casacore::ROArrayColumn<int> delaysCol(mwaTilePointing, "DELAYS");
    casacore::Array<int> delaysArr = delaysCol(0);
    casacore::Array<int>::contiter delaysArrPtr = delaysArr.cbegin();
    if (!hasDelays) {
      for (int i = 0; i != 16; ++i) delays[i] = delaysArrPtr[i];
      hasDelays = true;
    }

    aocommon::BandData bandData(ms.spectralWindow());
    centralFrequency = bandData.CentreFrequency();
    std::cout << "Central frequency of measurement set: "
              << round(centralFrequency * 1e-6 * 10.0) * 0.1 << " MHz.\n";
  }

  std::cout << "Delays: [";
  for (int i = 0; i != 16; ++i) {
    std::cout << delays[i];
    if (i != 15) std::cout << ',';
  }
  std::cout << "]\n";

  size_t width, height;
  double pixelSizeX, pixelSizeY, refRA, refDec, phaseCentreDL, phaseCentreDM;
  aocommon::FitsWriter writer;
  if (inpFitsname == 0) {
    // All sky
    width = 512;
    height = 512;
    pixelSizeX = 2.0 / (double)width;
    pixelSizeY = 2.0 / (double)height;
    refRA = 0.0;  // will be set later
    refDec = 0.0;
    phaseCentreDL = 0.0;
    phaseCentreDM = 0.0;
    double bandWidth = 1000000.0;

    if (!haveTime) {
      time = casacore::MEpoch(
          casacore::MVEpoch(casacore::Quantity(4.88193e+09, "s")));
      std::cout << "Warning: using random *wrong* start time, because no ms, "
                   "fits file or metafits file was specified.\n";
    }

    writer.SetFrequency(centralFrequency, bandWidth);
    writer.SetDate(time.getValue().get());
  } else {
    aocommon::FitsReader reader(inpFitsname);
    writer.SetMetadata(reader);

    width = reader.ImageWidth();
    height = reader.ImageHeight();
    pixelSizeX = reader.PixelSizeX();
    pixelSizeY = reader.PixelSizeY();
    refRA = reader.PhaseCentreRA();
    refDec = reader.PhaseCentreDec();
    phaseCentreDL = reader.LShift();
    phaseCentreDM = reader.MShift();
    centralFrequency = reader.Frequency();
    std::cout << "Using frequency from fits file: "
              << round(centralFrequency * 1e-6 * 10.0) * 0.1 << " MHz.\n";
    if (!haveTime) {
      time = casacore::MEpoch(casacore::MVEpoch(
          reader.DateObs()));  // casacore::Quantity(reader.DateObs(), "s")));
      std::cout << "Using time from Fits file: " << time << '\n';
      std::cout << "WARNING: Fits file records START of observation instead of "
                   "(as preferred) the MIDDLE of the observation.\n";
      std::cout << "WARNING: Suggest using -ms to get more accurate timing.\n";
    }
    if (!hasDelays) {
      int dDelays[16];
      MetaFitsFile::GetDelays(reader.FitsHandle(), dDelays);
      for (size_t i = 0; i != 16; ++i) delays[i] = dDelays[i];
      hasDelays = true;
    }
  }

  casacore::MeasFrame frame(arrayPos, time);
  const casacore::MDirection::Ref hadecRef(casacore::MDirection::HADEC, frame);
  const casacore::MDirection::Ref azelgeoRef(casacore::MDirection::AZELGEO,
                                             frame);
  const casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
  casacore::MPosition wgs =
      casacore::MPosition::Convert(arrayPos, casacore::MPosition::WGS84)();
  double arrLatitude = wgs.getValue().getLat();

  casacore::MDirection zenith(casacore::MVDirection(0.0, 0.0, 1.0), azelgeoRef);
  casacore::MDirection zenithHaDec =
      casacore::MDirection::Convert(zenith, hadecRef)();
  double zenithHa = zenithHaDec.getAngle().getValue()[0];
  double zenithDec = zenithHaDec.getAngle().getValue()[1];
  std::cout << "Zenith: "
            << (casacore::MDirection::Convert(zenith, j2000Ref)())
                       .getAngle()
                       .getValue()[0] *
                   180.0 / M_PI
            << " RA, " << zenithDec * 180.0 / M_PI << " dec, "
            << zenithHa * 180.0 / M_PI << " HA.\n";

  if (inpFitsname == 0) {
    refRA = (casacore::MDirection::Convert(zenith, j2000Ref)())
                .getAngle()
                .getValue()[0];
    refDec = zenithDec;
    writer.SetImageDimensions(width, height, refRA, refDec, pixelSizeX,
                              pixelSizeY);
  }

  std::cout << "Reference dir: " << refRA * 180.0 / M_PI << " RA, "
            << refDec * 180.0 / M_PI << " dec.\n";

  std::vector<double> outImage[8];
  double* imgPtr[8];
  for (size_t i = 0; i != 8; ++i) {
    outImage[i].resize(width * height);
    imgPtr[i] = &outImage[i][0];
  }

  if (doSquare) {
#ifndef CUDA_SUPPORT
    if (use2016) {
      MakeBeam<TileBeam2016, true>(
          imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec,
          arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame,
          phaseCentreDL, phaseCentreDM, doInterpolate);
    } else if (use2013) {
      MakeBeam<TileBeam2013, true>(
          imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec,
          arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame,
          phaseCentreDL, phaseCentreDM, doInterpolate);
    } else {
      MakeBeam<TileBeam2014, true>(
          imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec,
          arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame,
          phaseCentreDL, phaseCentreDM, doInterpolate);
    }
#else
    MakeBeam<TileBeam2016, true>(imgPtr, width, height, pixelSizeX, pixelSizeY,
                                 refRA, refDec, arrLatitude, zenithHa,
                                 zenithDec, centralFrequency, delays, frame,
                                 phaseCentreDL, phaseCentreDM, doInterpolate);
#endif
  } else {
#ifndef CUDA_SUPPORT
    if (use2016) {
      MakeBeam<TileBeam2016, false>(
          imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec,
          arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame,
          phaseCentreDL, phaseCentreDM, doInterpolate);
    } else if (use2013) {
      MakeBeam<TileBeam2013, false>(
          imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec,
          arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame,
          phaseCentreDL, phaseCentreDM, doInterpolate);
    } else {
      MakeBeam<TileBeam2014, false>(
          imgPtr, width, height, pixelSizeX, pixelSizeY, refRA, refDec,
          arrLatitude, zenithHa, zenithDec, centralFrequency, delays, frame,
          phaseCentreDL, phaseCentreDM, doInterpolate);
    }
#else
    MakeBeam<TileBeam2016, false>(imgPtr, width, height, pixelSizeX, pixelSizeY,
                                  refRA, refDec, arrLatitude, zenithHa,
                                  zenithDec, centralFrequency, delays, frame,
                                  phaseCentreDL, phaseCentreDM, doInterpolate);
#endif
  }

  std::cout << "\nWriting...\n";

  const string names[8] = {prefixName + "-xxr.fits", prefixName + "-xxi.fits",
                           prefixName + "-xyr.fits", prefixName + "-xyi.fits",
                           prefixName + "-yxr.fits", prefixName + "-yxi.fits",
                           prefixName + "-yyr.fits", prefixName + "-yyi.fits"};

  aocommon::PolarizationEnum linPols[4] = {
      aocommon::Polarization::XX, aocommon::Polarization::XY,
      aocommon::Polarization::YX, aocommon::Polarization::YY};
  for (size_t i = 0; i != 8; ++i) {
    writer.SetPolarization(linPols[i / 2]);
    writer.Write<double>(names[i], &outImage[i][0]);
  }
}
