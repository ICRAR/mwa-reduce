#include "imageadder.h"

#include <string>
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc < 5) {
    std::cerr << "Syntax: addimg [options] <outimage> <outweights> <inpimage1> "
                 "<inpbeam1> [<inpimage2> <inpbeam2> ...]\n"
                 "All images should be fits files. Use -1 as inpweight for "
                 "unity weight, or \"-c inpbeam-real inpbeam-imag\" for "
                 "complex beam. Options:\n"
                 "-threshold <value>\n"
                 "  Apply threshold\n"
                 "-no-weighting\n"
                 "  Don't use WSCIMGWG values.\n"
                 "-square\n"
                 "  Square before adding. Handy e.g. for creating a Stokes P "
                 "image.\n";
  } else {
    int argi = 1;
    ImageAdder adder;
    while (argi < argc && argv[argi][0] == '-') {
      std::string p(&argv[argi][1]);
      if (p == "threshold") {
        ++argi;
        adder.SetThreshold(atof(argv[argi]));
      } else if (p == "square") {
        adder.SetSquare(true);
      } else if (p == "no-weighting") {
        adder.SetWeighting(false);
      } else
        throw std::runtime_error("Invalid parameter");
      ++argi;
    }
    const char *outImageName = argv[argi];
    const char *outWeightName = argv[argi + 1];
    argi += 2;
    for (; argi + 1 < argc; argi += 2) {
      const char *inpImageName = argv[argi];
      const char *inpWeightName = argv[argi + 1];
      if (std::string(inpWeightName) == "-c") {
        adder.Add(inpImageName, argv[argi + 1], argv[argi + 2]);
        argi += 2;
      } else {
        adder.Add(inpImageName, inpWeightName);
      }
    }
    std::cout << '\n';
    adder.Finish(outImageName, outWeightName);
  }
}
