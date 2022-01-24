#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <set>

#include "model/model.h"

#include <aocommon/fits/fitsreader.h>
#include <aocommon/imagecoordinates.h>
#include <aocommon/uvector.h>

using aocommon::ImageCoordinates;

void GetModelFromImage(Model &model, const double *image, size_t width,
                       size_t height, double phaseCentreRA,
                       double phaseCentreDec, double pixelSizeX,
                       double pixelSizeY, double phaseCentreDL,
                       double phaseCentreDM, double spectralIndex,
                       double refFreq) {
  for (size_t y = 0; y != height; ++y) {
    for (size_t x = 0; x != width; ++x) {
      double value = image[y * width + x];
      if (value != 0.0 && std::isfinite(value)) {
        long double l, m;
        ImageCoordinates::XYToLM<long double>(x, y, pixelSizeX, pixelSizeY,
                                              width, height, l, m);
        l += phaseCentreDL;
        m += phaseCentreDM;
        ModelComponent component;
        long double ra, dec;
        ImageCoordinates::LMToRaDec<long double>(l, m, phaseCentreRA,
                                                 phaseCentreDec, ra, dec);
        std::stringstream nameStr;
        nameStr << "component" << model.SourceCount();
        component.SetSED(MeasuredSED(value, refFreq, spectralIndex,
                                     aocommon::Polarization::StokesI));
        component.SetPosRA(ra);
        component.SetPosDec(dec);

        ModelSource source;
        source.SetName(nameStr.str());
        source.AddComponent(component);
        model.AddSource(source);
      }
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage: fitsmodel [options] <output model> <fitsfile> "
                 "[spectral index] [ref freq MHz]\n"
                 "Turns components in fitsfile into a model.\nOptions:\n"
                 "\t-a <output areafile>\n"
                 "\t-d <merge distance>\n"
                 "\t-l <lower limit>\n"
                 "\t-in <areafile>\n"
                 "\t-out <areafile>\n";
  } else {
    size_t argi = 1;
    bool merge = false;
    double mergeDistance = 0.0;
    double limit = 0.0;
    std::string outputAreaFilename, insideFilename, outsideFilename;
    while (argv[argi][0] == '-') {
      std::string option(&argv[argi][1]);
      if (option == "a") {
        ++argi;
        outputAreaFilename = argv[argi];
      } else if (option == "d") {
        ++argi;
        merge = true;
        mergeDistance = atof(argv[argi]) * (M_PI / 180.0);
      } else if (option == "l") {
        ++argi;
        limit = atof(argv[argi]);
      } else if (option == "in") {
        ++argi;
        insideFilename = argv[argi];
      } else if (option == "out") {
        ++argi;
        outsideFilename = argv[argi];
      } else
        throw std::runtime_error("Invalid param");
      ++argi;
    }
    const char *modelFilename = argv[argi];
    ++argi;
    const char *fitsFilename = argv[argi];
    long double spectralIndex, refFreq;
    if (argc - argi > 1)
      spectralIndex = atof(argv[argi + 1]);
    else
      spectralIndex = 0.0;
    aocommon::FitsReader fitsReader(fitsFilename);
    if (argc - argi > 2)
      refFreq = atof(argv[argi + 2]) * 1000000;
    else
      refFreq = fitsReader.Frequency();
    const size_t width = fitsReader.ImageWidth(),
                 height = fitsReader.ImageHeight();
    aocommon::UVector<double> image(width * height);
    fitsReader.Read(&image[0]);

    Model model;
    GetModelFromImage(model, &image[0], width, height,
                      fitsReader.PhaseCentreRA(), fitsReader.PhaseCentreDec(),
                      fitsReader.PixelSizeX(), fitsReader.PixelSizeY(),
                      fitsReader.PhaseCentreDL(), fitsReader.PhaseCentreDM(),
                      spectralIndex, refFreq);

    model.Sort();

    if (merge) {
      std::set<size_t> sourcesToRemove;
      for (size_t i = 0; i != model.SourceCount(); ++i) {
        ModelSource &refSource = model.Source(i);
        for (size_t j = i + 1; j != model.SourceCount(); ++j) {
          const ModelSource &source2 = model.Source(j);
          double distance = aocommon::ImageCoordinates::AngularDistance(
              refSource.Peak().PosRA(), refSource.Peak().PosDec(),
              source2.Peak().PosRA(), source2.Peak().PosDec());
          if (distance <= mergeDistance) {
            sourcesToRemove.insert(j);
          }
        }
      }
      for (std::set<size_t>::reverse_iterator i = sourcesToRemove.rbegin();
           i != sourcesToRemove.rend(); ++i) {
        model.RemoveSource(*i);
      }
    }

    if (limit > 0.0) {
      std::cout << "Thresholding sources... " << std::flush;
      double totalFlux = 0.0, removedFlux = 0.0;
      for (size_t i = model.SourceCount(); i > 0; --i) {
        size_t sIndex = i - 1;
        double flux = model.Source(sIndex).TotalFlux(
            refFreq, aocommon::Polarization::StokesI);
        if (flux < limit) {
          model.RemoveSource(sIndex);
          removedFlux += flux;
        }
        totalFlux += flux;
      }
      std::cout << "DONE (removed: " << removedFlux << " / " << totalFlux
                << " Jy)\n";
    }

    std::cout << "Writing model with " << model.SourceCount() << " sources... "
              << std::flush;
    model.Save(modelFilename);
    std::cout << "DONE\n";
  }
}
