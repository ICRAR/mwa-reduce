#include <iostream>

#include "model/bbsmodel.h"

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "Syntax: bbs2model [options] <input.txt> <output.txt>\n"
                 "Options:\n"
                 " -source <name>\n"
                 "    Let all components belong to a single source with the "
                 "given name.\n";
  } else {
    int argi = 1;
    std::string sourcename;
    while (argv[argi][0] == '-') {
      std::string p(&argv[argi][1]);
      if (p == "source") {
        ++argi;
        sourcename = argv[argi];
      } else
        throw std::runtime_error("Unknown parameter");
      ++argi;
    }
    BBSModel::Read(argv[argi], sourcename).Save(argv[argi + 1]);
  }
  return 0;
}
