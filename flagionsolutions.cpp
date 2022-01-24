#include <iostream>

#include "ionsolutionfile.h"

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage: flagionsolutions <ionsolutions.bin>\n"
                 "Gives flag statistics.\n";
  } else {
    size_t argi = 1;
    while (argv[argi][0] == '-') {
      // std::string p(&argv[argi][1]);
      throw std::runtime_error("Bad option");
      ++argi;
    }

    IonSolutionFile file;
    file.OpenForReading(argv[1]);

    for (size_t i = 0; i != file.DirectionCount(); ++i) {
      std::string clusterName;
      std::vector<std::string> sourceNames;
      file.ReadClusterMetaInfo(clusterName, sourceNames);
    }

    size_t count = 0, finiteCount = 0;
    for (size_t direction = 0; direction != file.DirectionCount();
         ++direction) {
      for (size_t interval = 0; interval != file.IntervalCount(); ++interval) {
        for (size_t cb = 0; cb != file.ChannelBlockCount(); ++cb) {
          for (size_t p = 0; p != file.PolarizationCount(); ++p) {
            double gain = file.ReadSolution(IonSolutionFile::GainSolution,
                                            interval, cb, p, direction);
            ++count;
            if (std::isfinite(gain)) ++finiteCount;
          }
        }
      }
    }
    double perc = 100.0 * double(count - finiteCount) / double(count);
    std::cout << "Flagged count: " << (count - finiteCount) << " / " << count
              << " (" << perc << "%)\n";
  }
}
