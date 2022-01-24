#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/uvector.h>

#include <boost/filesystem/operations.hpp>

#include <stdexcept>
#include <string>
#include <iostream>
#include <complex>
#include <limits>

void read(const aocommon::FitsReader& templateReader,
          const std::string& filename, aocommon::UVector<double>& data) {
  std::cout << "- " << filename << '\n';
  if (boost::filesystem::exists(filename)) {
    aocommon::FitsReader inpReader(filename);
    size_t width = inpReader.ImageWidth(), height = inpReader.ImageHeight();

    if (width != templateReader.ImageWidth() ||
        height != templateReader.ImageHeight())
      throw std::runtime_error("Not all images had same size");
    data.resize(width * height);
    inpReader.Read<double>(&data[0]);
  } else {
    std::cout << "Warning: file '" << filename
              << "' did not exist, assuming it is zero!\n";
    data.assign(templateReader.ImageWidth() * templateReader.ImageHeight(),
                0.0);
  }
}

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::cout
        << "Syntax:\n"
           "  pbaddimg [options] <out-prefix> {<image-prefix> <image-postfix> "
           "<beam-prefix> [..]}\n\n"
           "Options:\n"
           "-withp\n"
           "  Also save total polarization image.\n"
           "-intermediate <nopb-image-prefix> <nopb-weight-prefix>"
           "  Save the intermediate images and image weights before beam "
           "correction. These can be added to other uncorrected images.\n"
           "\n"
           "Example:\n"
           "  pbaddimg integrated-stokes wsclean-A image.fits beam-A wsclean-B "
           "image.fits beam-B\n"
           "This command will look for image names wsclean-A-XX-image.fits, "
           "beam-A-xxr.fits, ... and save to integrated-stokes-i.fits, ...\n";
    return -1;
  }

  std::string intermImagePrefix, intermWeightPrefix;
  bool withStokesP = false;

  size_t argi = 1;
  while (argv[argi][0] == '-') {
    const std::string param(&argv[argi][1]);
    if (param == "withp") {
      withStokesP = true;
    } else if (param == "intermediate") {
      ++argi;
      intermImagePrefix = argv[argi];
      ++argi;
      intermWeightPrefix = argv[argi];
    } else
      throw std::runtime_error("Invalid command line parameter: " + param);
    ++argi;
  }

  const std::string outPrefix = argv[argi];
  ++argi;

  size_t inputCount = (argc - argi) / 3;

  aocommon::UVector<std::complex<double>> outputJonesData, beamJonesData;

  std::string firstFilename(std::string(argv[argi]) + "-XX-" + argv[argi + 1]);
  aocommon::FitsReader firstImage(firstFilename);
  const size_t width = firstImage.ImageWidth(),
               height = firstImage.ImageHeight(), imgSize = width * height;

  for (size_t imageSetIndex = 0; imageSetIndex != inputCount; ++imageSetIndex) {
    const std::string imagePrefix = argv[argi], imagePostfix = argv[argi + 1],
                      beamPrefix = argv[argi + 2];
    argi += 3;
    std::cout << "Processing image set with prefix '" << imagePrefix << "' ("
              << (imageSetIndex + 1) << " / " << inputCount << ")...\n";

    std::string inpFilenames[4] = {imagePrefix + "-XX-" + imagePostfix,
                                   imagePrefix + "-XY-" + imagePostfix,
                                   imagePrefix + "-XYi-" + imagePostfix,
                                   imagePrefix + "-YY-" + imagePostfix},
                beamFilename[8] = {
                    beamPrefix + "-xxr.fits", beamPrefix + "-xxi.fits",
                    beamPrefix + "-xyr.fits", beamPrefix + "-xyi.fits",
                    beamPrefix + "-yxr.fits", beamPrefix + "-yxi.fits",
                    beamPrefix + "-yyr.fits", beamPrefix + "-yyi.fits"};
    aocommon::UVector<double> inputData[4], beamData[8];
    for (size_t i = 0; i != 4; ++i)
      read(firstImage, inpFilenames[i], inputData[i]);
    for (size_t i = 0; i != 8; ++i)
      read(firstImage, beamFilename[i], beamData[i]);

    if (outputJonesData.empty()) {
      outputJonesData.assign(imgSize * 4, std::complex<double>(0.0));
      beamJonesData.assign(imgSize * 4, std::complex<double>(0.0));
    } else if (imgSize * 4 != outputJonesData.size())
      throw std::runtime_error("Image sizes did not match!");

    aocommon::UVector<std::complex<double>>::iterator
        jonesIter = outputJonesData.begin(),
        beamIter = beamJonesData.begin();
    for (size_t i = 0; i != imgSize; ++i) {
      std::complex<double> beamValues[4];
      beamValues[0] = std::complex<double>(beamData[0][i], beamData[1][i]);
      beamValues[1] = std::complex<double>(beamData[2][i], beamData[3][i]);
      beamValues[2] = std::complex<double>(beamData[4][i], beamData[5][i]);
      beamValues[3] = std::complex<double>(beamData[6][i], beamData[7][i]);

      std::complex<double> imgValues[4];
      imgValues[0] = inputData[0][i];
      imgValues[1] = std::complex<double>(inputData[1][i], inputData[2][i]);
      imgValues[2] = std::conj(imgValues[1]);
      imgValues[3] = inputData[3][i];

      // Calculate Flux += w B* V B  (from: w (B* B) B^-1 V B*^-1 (B* B))
      // 'w' is 1 for now.
      std::complex<double> tempValues[4];
      aocommon::Matrix2x2::HermATimesB(tempValues, beamValues, imgValues);
      aocommon::Matrix2x2::PlusATimesB(&*jonesIter, tempValues, beamValues);

      // Calculate Weight += w (B* B) (B* B)
      aocommon::Matrix2x2::HermATimesB(tempValues, beamValues, beamValues);
      aocommon::Matrix2x2::HermATimesB(beamValues, tempValues, tempValues);
      aocommon::Matrix2x2::Add(&*beamIter, beamValues);

      jonesIter += 4;
      beamIter += 4;
    }
  }

  aocommon::FitsWriter writer(firstImage);

  if (!intermImagePrefix.empty()) {
    aocommon::UVector<double> intermediateImages[2];
    intermediateImages[0].resize(imgSize);
    intermediateImages[1].resize(imgSize);

    std::cout << "Saving intermediate files...\n";
    std::string outImageFilenames[4] = {intermImagePrefix + "-XX",
                                        intermImagePrefix + "-XY",
                                        intermImagePrefix + "-YX",
                                        intermImagePrefix + "-YY"},
                outWeightFilenames[4] = {
                    intermWeightPrefix + "-XX", intermWeightPrefix + "-XY",
                    intermWeightPrefix + "-YX", intermWeightPrefix + "-YY"};
    aocommon::PolarizationEnum intermPols[4] = {
        aocommon::Polarization::XX, aocommon::Polarization::XY,
        aocommon::Polarization::YX, aocommon::Polarization::YY};

    for (size_t p = 0; p != 4; ++p) {
      writer.SetPolarization(intermPols[p]);

      aocommon::UVector<std::complex<double>>::iterator jonesIter =
          outputJonesData.begin() + p;
      for (size_t i = 0; i != imgSize; ++i) {
        intermediateImages[0][i] = jonesIter->real();
        intermediateImages[1][i] = jonesIter->imag();
        jonesIter += 4;
      }

      writer.Write<double>(outImageFilenames[p] + ".fits",
                           &intermediateImages[0][0]);
      writer.Write<double>(outImageFilenames[p] + "i.fits",
                           &intermediateImages[1][0]);

      aocommon::UVector<std::complex<double>>::iterator beamIter =
          beamJonesData.begin() + p;
      for (size_t i = 0; i != imgSize; ++i) {
        intermediateImages[0][i] = beamIter->real();
        intermediateImages[1][i] = beamIter->imag();
        beamIter += 4;
      }
      writer.Write<double>(outWeightFilenames[p] + ".fits",
                           &intermediateImages[0][0]);
      writer.Write<double>(outWeightFilenames[p] + "i.fits",
                           &intermediateImages[1][0]);
    }
  }

  std::cout << "Dividing by weights and primary beams...\n";

  aocommon::UVector<double> stokesOutputImages[4];
  for (size_t p = 0; p != 4; ++p) stokesOutputImages[p].resize(imgSize);

  aocommon::UVector<std::complex<double>>::iterator jonesIter =
                                                        outputJonesData.begin(),
                                                    beamIter =
                                                        beamJonesData.begin();
  for (size_t i = 0; i != imgSize; ++i) {
    double stokesMatrix[4];

    // Calculate: D sum(W*W)^-1
    if (aocommon::Matrix2x2::Invert(&*beamIter)) {
      std::complex<double> correctedLinear[4];
      aocommon::Matrix2x2::ATimesB(correctedLinear, &*beamIter, &*jonesIter);
      aocommon::Polarization::LinearToStokes(correctedLinear, stokesMatrix);
    } else {
      for (size_t p = 0; p != 4; ++p)
        stokesMatrix[p] = std::numeric_limits<double>::quiet_NaN();
    }

    for (size_t p = 0; p != 4; ++p) stokesOutputImages[p][i] = stokesMatrix[p];

    jonesIter += 4;
    beamIter += 4;
  }

  std::string outFilenames[4] = {outPrefix + "-I.fits", outPrefix + "-Q.fits",
                                 outPrefix + "-U.fits", outPrefix + "-V.fits"};
  aocommon::PolarizationEnum stokesPols[4] = {
      aocommon::Polarization::StokesI, aocommon::Polarization::StokesQ,
      aocommon::Polarization::StokesU, aocommon::Polarization::StokesV};

  std::cout << "Saving files...\n";
  for (size_t p = 0; p != 4; ++p) {
    writer.SetPolarization(stokesPols[p]);
    writer.Write<double>(outFilenames[p], &stokesOutputImages[p][0]);
  }

  if (withStokesP) {
    for (size_t i = 0; i != imgSize; ++i) {
      double q = stokesOutputImages[1][i], u = stokesOutputImages[2][i];
      stokesOutputImages[1][i] = sqrt(q * q + u * u);
    }
    writer.SetPolarization(aocommon::Polarization::StokesI);
    writer.Write<double>(outPrefix + "-P.fits", &stokesOutputImages[1][0]);
  }
}
