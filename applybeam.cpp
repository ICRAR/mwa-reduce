#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>

#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>

using namespace aocommon;

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout << "Syntax: applybeam <inpfits> <beamfits> <outfits>\n";
    return 0;
  }
  const char *inpFits = argv[1];
  const char *beamFits = argv[2];
  const char *outFits = argv[3];

  FitsReader inpReader(inpFits);
  size_t width = inpReader.ImageWidth(), height = inpReader.ImageHeight();
  FitsReader beamReader(beamFits);
  if (beamReader.ImageWidth() != width || beamReader.ImageHeight() != height)
    throw std::runtime_error("Beam and image do not have same size!");

  std::vector<double> inpImage(width * height), beamImage(width * height);

  inpReader.Read<double>(&inpImage[0]);
  beamReader.Read<double>(&beamImage[0]);

  std::vector<double>::iterator beamIter = beamImage.begin();
  for (std::vector<double>::iterator i = inpImage.begin(); i != inpImage.end();
       ++i) {
    *i /= *beamIter * *beamIter;
    ++beamIter;
  }

  FitsWriter writer(inpReader);
  writer.Write<double>(outFits, &inpImage[0]);
}
