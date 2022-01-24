#include <iostream>
#include <stdexcept>
#include <vector>
#include <random>

#include <boost/optional/optional.hpp>

#include "image.h"
#include "progressbar.h"

#include "math/modelrenderer.h"

#include "render/interpolatingrenderer.h"

#include "model/bbsmodel.h"
#include "model/model.h"
#include "model/modelparser.h"

#include <aocommon/banddata.h>
#include <aocommon/parallelfor.h>
#include <aocommon/threadpool.h>
#include <aocommon/lane.h>
#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>
#include <aocommon/units/angle.h>

using aocommon::units::Angle;

struct ImageSettings {
  double startFrequency;
  double endFrequency;
  size_t width;
  size_t height;
  double ra;
  double dec;
  double pixelSizeX;
  double pixelSizeY;
  double pdl;
  double pdm;
};

void renderSources(aocommon::Lane<ModelSource>* sourceLane, Image* image,
                   size_t windowSize, const ImageSettings& settings) {
  InterpolatingRenderer renderer(windowSize);
  // ModelRenderer renderer(ra, dec, *pixelSizeX, *pixelSizeY, pdl, pdm );
  ModelSource source;
  while (sourceLane->read(source)) {
    for (const ModelComponent& comp : source) {
      const float flux = comp.SED().IntegratedFlux(
          settings.startFrequency, settings.endFrequency,
          aocommon::Polarization::StokesI);
      if (!std::isfinite(flux)) {
        std::cout << "Evaluating the spectral for source " + source.Name()
                  << " resulted in a non-finite value.\n";
        throw std::runtime_error("Evaluating the spectral for source " +
                                 source.Name() +
                                 " resulted in a non-finite value");
      }
      if (comp.Type() != ModelComponent::PointSource) {
        ModelRenderer::RenderGaussianComponent(
            image->data(), settings.width, settings.height, settings.ra,
            settings.dec, settings.pixelSizeX, settings.pixelSizeY,
            settings.pdl, settings.pdm, comp.PosRA(), comp.PosDec(),
            comp.MajorAxis(), comp.MinorAxis(), comp.PositionAngle(), flux);
      } else {
        double l, m;
        float x, y;
        aocommon::ImageCoordinates::RaDecToLM<double>(
            comp.PosRA(), comp.PosDec(), settings.ra, settings.dec, l, m);
        l += settings.pdl;
        m += settings.pdm;
        aocommon::ImageCoordinates::LMToXYfloat<float>(
            l, m, settings.pixelSizeX, settings.pixelSizeY, settings.width,
            settings.height, x, y);
        if (windowSize)
          renderer.RenderWindowedSource(image->data(), settings.width,
                                        settings.height, flux, x, y);
        else
          renderer.RenderSource(image->data(), settings.width, settings.height,
                                flux, x, y);
      }
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc == 1)
    std::cout
        << "syntax: render [options] <model>\n"
           "Options:\n"
           "  [-n <noiselevel>]"
           "    Add noise to the image.\n"
           "  [-t templatefits]"
           "    Use the image size and coordinate system from this file.\n"
           "  [-o <outputfits>]"
           "    Save the output to the given file.\n"
           "  [-r [-beam <maj> <min> <pa>]]"
           "    Use a restoring beam with optionally a manually given size.\n"
           "  [-a]"
           "    Add the sources to the given template file.\n"
           "  [-centre <ra> <dec>]"
           "    Change the centre of the fits image.\n"
           "  [-size <width> <height>]"
           "    Set the size of the output image.\n"
           "  [-scale <scale>]"
           "    Set the pixel scale of the output image.\n"
           "  [-frequency <valueHz>]\n"
           "    Set the frequency of the output image.\n"
           "  [-i [-window <size>]]\n"
           "    Interpolate sources that fall between pixels (ignored when -r "
           "is also given)\n";
  else {
    std::string templateFits;
    std::string outputFitsName;
    std::string ionOutPrefix;
    bool restore = false, addToTemplate = false, hasManualBeam = false;
    int argi = 1;
    double ra = 0.0, dec = 0.0, pdl = 0.0, pdm = 0.0;
    boost::optional<double> pixelSizeX, pixelSizeY;
    double noise = 0.0;
    double beamMaj = 2.0 * (M_PI / 180.0 / 60.0),
           beamMin = 2.0 * (M_PI / 180.0 / 60.0), beamPA = 0.0;
    size_t sizeWidth = 0, sizeHeight = 0;
    double setFrequency = 0.0;
    bool interpolate = false;
    size_t windowSize = 0;
    while (argi < argc && argv[argi][0] == '-') {
      std::string param(&argv[argi][1]);
      if (param == "t") {
        ++argi;
        templateFits = argv[argi];
      } else if (param == "r") {
        restore = true;
      } else if (param == "n") {
        ++argi;
        noise = atof(argv[argi]);
      } else if (param == "beam") {
        hasManualBeam = true;
        ++argi;
        beamMaj =
            Angle::Parse(argv[argi], "beam major axis", Angle::kArcseconds);
        ++argi;
        beamMin =
            Angle::Parse(argv[argi], "beam minor axis", Angle::kArcseconds);
        ++argi;
        beamPA =
            Angle::Parse(argv[argi], "beam position angle", Angle::kDegrees);
      } else if (param == "a") {
        addToTemplate = true;
      } else if (param == "o") {
        ++argi;
        outputFitsName = argv[argi];
      } else if (param == "centre") {
        ++argi;
        ra = RaDecCoord::ParseRA(argv[argi]);
        ++argi;
        dec = RaDecCoord::ParseDec(argv[argi]);
      } else if (param == "size") {
        sizeWidth = atoi(argv[argi + 1]);
        sizeHeight = atoi(argv[argi + 2]);
        argi += 2;
      } else if (param == "scale") {
        ++argi;
        pixelSizeX = Angle::Parse(argv[argi], "scale", Angle::kDegrees);
        pixelSizeY = *pixelSizeX;
      } else if (param == "frequency") {
        ++argi;
        setFrequency = atof(argv[argi]);
      } else if (param == "i") {
        interpolate = true;
      } else if (param == "window") {
        ++argi;
        windowSize = atoi(argv[argi]);
      } else
        throw std::runtime_error("Invalid param");
      ++argi;
    }

    size_t inputWidth = 4096;
    size_t inputHeight = 4096;
    double bandwidth = 1000000.0, dateObs = 0.0, frequency = 150000000.0;

    std::unique_ptr<aocommon::FitsWriter> writer;
    std::unique_ptr<aocommon::FitsReader> reader;
    Image image;
    if (!templateFits.empty()) {
      double wscImgWeight = 0.0;
      reader.reset(new aocommon::FitsReader(templateFits));
      inputWidth = reader->ImageWidth();
      inputHeight = reader->ImageHeight();
      image = Image(inputWidth, inputHeight, 0.0);
      ra = reader->PhaseCentreRA();
      dec = reader->PhaseCentreDec();
      pdl = reader->PhaseCentreDL();
      pdm = reader->PhaseCentreDM();
      if (!pixelSizeX) pixelSizeX = reader->PixelSizeX();
      if (!pixelSizeY) pixelSizeY = reader->PixelSizeY();
      bandwidth = reader->Bandwidth();
      dateObs = reader->DateObs();
      frequency = reader->Frequency();
      if (reader->HasBeam() && !hasManualBeam) {
        beamMaj = reader->BeamMajorAxisRad();
        beamMin = reader->BeamMinorAxisRad();
        beamPA = reader->BeamPositionAngle();
      }
      reader->ReadDoubleKeyIfExists("WSCIMGWG", wscImgWeight);
      if (addToTemplate) reader->Read(&image[0]);

      writer.reset(new aocommon::FitsWriter(*reader));
      if (wscImgWeight != 0.0)
        writer->SetExtraKeyword("WSCIMGWG", wscImgWeight);
    } else {
      inputWidth = sizeWidth;
      inputHeight = sizeHeight;
      image = Image(inputWidth, inputHeight, 0.0);
      writer.reset(new aocommon::FitsWriter());
    }

    if (sizeWidth != 0 && sizeHeight != 0) {
      if (sizeWidth > inputWidth && sizeHeight > inputHeight) {
        image = image.Untrim(sizeWidth, sizeHeight);
        inputWidth = sizeWidth;
        inputHeight = sizeHeight;
      }
    }

    if (setFrequency != 0.0) {
      frequency = setFrequency;
      writer->SetFrequency(setFrequency, writer->Bandwidth());
    }

    ModelRenderer renderer(ra, dec, *pixelSizeX, *pixelSizeY, pdl, pdm);
    if (noise != 0.0) {
      std::random_device rd;
      std::mt19937 rnd(rd());
      std::normal_distribution<double> dist(0.0, noise);
      for (size_t i = 0; i != inputWidth * inputHeight; ++i)
        image[i] += dist(rnd);
    }
    if (restore) {
      const Model model(argv[argi]);
      renderer.Restore(&image[0], inputWidth, inputHeight, model, beamMaj,
                       beamMin, beamPA, frequency - bandwidth * 0.5,
                       frequency + bandwidth * 0.5,
                       aocommon::Polarization::StokesI);
    } else if (interpolate) {
      std::cout << "Rendering sources...\n";
      const size_t ncpus = aocommon::system::ProcessorCount();
      std::vector<Image> images(ncpus);
      std::vector<std::thread> threads;
      ImageSettings settings;
      settings.startFrequency = frequency - bandwidth * 0.5;
      settings.endFrequency = frequency + bandwidth * 0.5;
      settings.width = inputWidth;
      settings.height = inputHeight;
      settings.ra = ra;
      settings.dec = dec;
      settings.pixelSizeX = *pixelSizeX;
      settings.pixelSizeY = *pixelSizeY;
      settings.pdl = pdl;
      settings.pdm = pdm;
      aocommon::Lane<ModelSource> sourceQueue(ncpus);
      for (Image& threadImage : images) {
        threadImage = Image(inputWidth, inputHeight, 0.0);
        threads.emplace_back(renderSources, &sourceQueue, &threadImage,
                             windowSize, settings);
      }
      std::ifstream file(argv[argi]);
      const bool isModelFormat = ModelParser::IsInModelFormat(file);
      file.close();
      auto processFunction = [&sourceQueue](const ModelSource& s) {
        sourceQueue.write(s);
      };
      if (isModelFormat) {
        std::ifstream modelFile(argv[argi]);
        ModelParser parser;
        parser.Stream(modelFile, processFunction);
      } else {
        BBSModel::Read(argv[argi], processFunction);
      }
      std::cout << "Finishing...\n";
      sourceQueue.write_end();
      for (std::thread& t : threads) t.join();
      for (const Image& threadImage : images) {
        for (size_t i = 0; i != threadImage.size(); ++i)
          image[i] += threadImage[i];
      }
    } else {
      const Model model(argv[argi]);
      renderer.RenderModel(image.data(), inputWidth, inputHeight, model,
                           frequency - bandwidth * 0.5,
                           frequency + bandwidth * 0.5,
                           aocommon::Polarization::StokesI);
    }

    if (!outputFitsName.empty()) {
      writer->SetImageDimensions(inputWidth, inputHeight, ra, dec, *pixelSizeX,
                                 *pixelSizeY);
      writer->SetFrequency(frequency, bandwidth);
      writer->SetDate(dateObs);
      writer->Write(outputFitsName, &image[0]);
    }
  }
}
