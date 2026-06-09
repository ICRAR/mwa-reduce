#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/matrix2x2.h>

#include <boost/filesystem/operations.hpp>

#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <complex>
#include <limits>

void read(const aocommon::FitsReader& templateReader,
          const std::string& filename, std::vector<double>& data) {
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
    data.resize(templateReader.ImageWidth() * templateReader.ImageHeight(),
                0.0);
  }
}

void correctStokesI(const std::string& imagePrefix,
                    const std::string& imagePostfix,
                    const std::string& beamPrefix,
                    const std::string& outPrefix) {
  const std::string beamFilename[8] = {
      beamPrefix + "-xxr.fits", beamPrefix + "-xxi.fits",
      beamPrefix + "-xyr.fits", beamPrefix + "-xyi.fits",
      beamPrefix + "-yxr.fits", beamPrefix + "-yxi.fits",
      beamPrefix + "-yyr.fits", beamPrefix + "-yyi.fits"};

  std::string filename(imagePrefix + "-" + imagePostfix);
  aocommon::FitsReader inputReader(filename);

  const size_t width = inputReader.ImageWidth(),
               height = inputReader.ImageHeight(), imgSize = width * height;

  std::vector<double> inputData, weightData(imgSize), beamData[8];
  read(inputReader, filename, inputData);
  for (size_t i = 0; i != 8; ++i)
    read(inputReader, beamFilename[i], beamData[i]);

  for (size_t i = 0; i != imgSize; ++i) {
    aocommon::MC2x2 beamValues{
    std::complex<double>(beamData[0][i], beamData[1][i]),
    std::complex<double>(beamData[2][i], beamData[3][i]),
    std::complex<double>(beamData[4][i], beamData[5][i]),
    std::complex<double>(beamData[6][i], beamData[7][i])};

    aocommon::MC2x2 squared = beamValues.MultiplyHerm(beamValues);
    if (squared.Invert()) {
      double factor = 0.5 * (squared.Get(0).real() + squared.Get(3).real());
      inputData[i] = inputData[i] * factor;
      weightData[i] = 1.0 / (factor * factor);
    } else {
      inputData[i] = std::numeric_limits<double>::quiet_NaN();
      weightData[i] = 0.0;
    }
  }

  aocommon::FitsWriter writer(inputReader);
  writer.Write<double>(outPrefix + "-I.fits", inputData.data());
  writer.Write<double>(outPrefix + "-weight.fits", weightData.data());
}

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::cout
        << "Syntax:\n"
           "pbcorrect [-from-stokes/-stokesi] [-uncorrect] <image-prefix> "
           "<image-postfix> <beam-prefix> <out-prefix>\nFor "
           "example:\n\tpbcorrect wsclean image.fits beam stokes\n"
           "will use images like wsclean-XX-image.fits, beam-xxr.fits and save "
           "to stokes-i.fits, ...\n"
           "Uncorrect will do the inverse operation, so it will read the "
           "images indicated by out-prefix,\n"
           "and write the linear images indicated by the image pre- & "
           "postfixes.\n";
    return -1;
  }

  bool doUncorrect = false, fromStokes = false, stokesIOnly = false;
  size_t argi = 1;
  while (argv[argi][0] == '-') {
    const std::string param(&argv[argi][1]);
    if (param == "uncorrect")
      doUncorrect = true;
    else if (param == "from-stokes")
      fromStokes = true;
    else if (param == "stokesi")
      stokesIOnly = true;
    else
      throw std::runtime_error("Invalid command line parameter: " + param);
    ++argi;
  }
  const std::string imagePrefix = argv[argi], imagePostfix = argv[argi + 1],
                    beamPrefix = argv[argi + 2], outPrefix = argv[argi + 3];

  if (stokesIOnly)
    correctStokesI(imagePrefix, imagePostfix, beamPrefix, outPrefix);
  else {
    std::string linFilenames[4];
    if (fromStokes) {
      linFilenames[0] = imagePrefix + "-I-" + imagePostfix;
      linFilenames[1] = imagePrefix + "-Q-" + imagePostfix;
      linFilenames[2] = imagePrefix + "-U-" + imagePostfix;
      linFilenames[3] = imagePrefix + "-V-" + imagePostfix;
    } else {
      linFilenames[0] = imagePrefix + "-XX-" + imagePostfix;
      linFilenames[1] = imagePrefix + "-XY-" + imagePostfix;
      linFilenames[2] = imagePrefix + "-XYi-" + imagePostfix;
      linFilenames[3] = imagePrefix + "-YY-" + imagePostfix;
    }

    std::string beamFilename[8] = {beamPrefix + "-xxr.fits",
                                   beamPrefix + "-xxi.fits",
                                   beamPrefix + "-xyr.fits",
                                   beamPrefix + "-xyi.fits",
                                   beamPrefix + "-yxr.fits",
                                   beamPrefix + "-yxi.fits",
                                   beamPrefix + "-yyr.fits",
                                   beamPrefix + "-yyi.fits"},
                stokesFilenames[4] = {
                    outPrefix + "-I.fits", outPrefix + "-Q.fits",
                    outPrefix + "-U.fits", outPrefix + "-V.fits"};
    std::string *inpFilenames = doUncorrect ? stokesFilenames : linFilenames,
                *outFilenames = doUncorrect ? linFilenames : stokesFilenames;
    aocommon::FitsReader firstImage(inpFilenames[0]);
    std::vector<double> inputData[4], beamData[8];
    for (size_t i = 0; i != 4; ++i)
      read(firstImage, inpFilenames[i], inputData[i]);
    for (size_t i = 0; i != 8; ++i)
      read(firstImage, beamFilename[i], beamData[i]);

    const size_t width = firstImage.ImageWidth(),
                 height = firstImage.ImageHeight(), imgSize = width * height;

    for (size_t i = 0; i != imgSize; ++i) {
      std::complex<double> beamValues[4];
      beamValues[0] = std::complex<double>(beamData[0][i], beamData[1][i]);
      beamValues[1] = std::complex<double>(beamData[2][i], beamData[3][i]);
      beamValues[2] = std::complex<double>(beamData[4][i], beamData[5][i]);
      beamValues[3] = std::complex<double>(beamData[6][i], beamData[7][i]);

      if (doUncorrect) {
        double stokes[4];
        for (size_t p = 0; p != 4; ++p) stokes[p] = inputData[p][i];

        std::complex<double> tempValues[4], out[4];
        aocommon::Polarization::StokesToLinear(stokes, out);
        aocommon::Matrix2x2::ATimesB(tempValues, beamValues, out);
        aocommon::Matrix2x2::ATimesHermB(out, tempValues, beamValues);

        inputData[0][i] = out[0].real();
        inputData[1][i] = 0.5 * (out[1].real() + out[2].real());
        inputData[2][i] = 0.5 * (out[1].imag() - out[2].imag());
        inputData[3][i] = out[3].real();

        for (size_t p = 0; p != 4; ++p)
          if (!std::isfinite(inputData[p][i])) inputData[p][i] = 0.0;
      } else {
        std::complex<double> imgValues[4];

        if (fromStokes) {
          double temp[4];
          temp[0] = inputData[0][i];
          temp[1] = inputData[1][i];
          temp[2] = inputData[2][i];
          temp[3] = inputData[3][i];
          aocommon::Polarization::StokesToLinear(temp, imgValues);
        } else {
          imgValues[0] = inputData[0][i];
          imgValues[1] = std::complex<double>(inputData[1][i], inputData[2][i]);
          imgValues[2] = std::conj(imgValues[1]);
          imgValues[3] = inputData[3][i];
        }

        if (aocommon::Matrix2x2::Invert(beamValues)) {
          std::complex<double> tempValues[4];
          aocommon::Matrix2x2::ATimesB(tempValues, beamValues, imgValues);
          aocommon::Matrix2x2::ATimesHermB(imgValues, tempValues, beamValues);

          double outputValues[4];
          aocommon::Polarization::LinearToStokes(imgValues, outputValues);
          for (size_t p = 0; p != 4; ++p) inputData[p][i] = outputValues[p];
        } else {
          for (size_t p = 0; p != 4; ++p)
            inputData[p][i] = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }

    aocommon::PolarizationEnum linPols[4] = {aocommon::Polarization::XX,
                                             aocommon::Polarization::XY,
                                             aocommon::Polarization::XY,
                                             aocommon::Polarization::YY},
                               stokesPols[4] = {
                                   aocommon::Polarization::StokesI,
                                   aocommon::Polarization::StokesQ,
                                   aocommon::Polarization::StokesU,
                                   aocommon::Polarization::StokesV};

    aocommon::FitsWriter writer(firstImage);

    for (size_t p = 0; p != 4; ++p) {
      writer.SetPolarization(doUncorrect ? linPols[p] : stokesPols[p]);
      writer.Write<double>(outFilenames[p], &inputData[p][0]);
    }
  }
}
