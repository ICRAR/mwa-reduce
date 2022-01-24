#include "solutionflagfile.h"
#include "solutionfile.h"

#include <iostream>
#include <limits>
#include <boost/filesystem/operations.hpp>

double phaseDist(double phase1, double phase2) {
  double dist = fmod(fabs(phase2 - phase1), 2.0 * M_PI);
  if (dist > M_PI) dist = M_PI - dist;
  return dist;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: flagsolutions [-c <second-sol.bin>] [-a ant] [-sum] "
                 "<solutions.bin> [<solutions-flag-file.txt>]\n"
                 "Will flag solutions in the solution file as specified by the "
                 "flag file.\n";
  } else {
    size_t argi = 1;
    const char* comparisonFilename = 0;
    int displayAntenna = -1;
    bool sum = false;
    while (argv[argi][0] == '-') {
      std::string p(&argv[argi][1]);
      if (p == "c") {
        ++argi;
        comparisonFilename = argv[argi];
      } else if (p == "a") {
        ++argi;
        displayAntenna = atoi(argv[argi]);
      } else if (p == "sum") {
        sum = true;
      } else
        throw std::runtime_error("Bad option");
      ++argi;
    }
    const char* solutionFilename = argv[argi];

    SolutionFile solutionFile;
    solutionFile.OpenForReading(solutionFilename);

    if (displayAntenna >= 0 || sum) {
      for (size_t interval = 0; interval != solutionFile.IntervalCount();
           ++interval) {
        std::cout << interval << ' ';
        std::complex<double> sumValue = 0.0;
        size_t sumCount = 0;
        for (size_t a = 0; a != solutionFile.AntennaCount(); ++a) {
          for (size_t ch = 0; ch != solutionFile.ChannelCount(); ++ch) {
            for (size_t p = 0; p != 4; ++p) {
              std::complex<double> val = solutionFile.ReadNextSolution();
              if (sum) {
                if (std::isfinite(val.real()) && val != 0.0) {
                  sumValue += val;
                  ++sumCount;
                }
              } else if (int(a) == displayAntenna) {
                std::cout << val.real() << ' ' << val.imag() << ' '
                          << std::abs(val) << ' ' << std::arg(val) << ' ';
              }
            }
            if (!sum && int(a) == displayAntenna) std::cout << '\n';
          }
        }
        if (sum) {
          std::cout << sumValue.real() / sumCount << ' '
                    << sumValue.imag() / sumCount << ' '
                    << std::abs(sumValue) / sumCount << ' '
                    << std::arg(sumValue) << '\n';
        }
      }
    } else {
      std::cout << solutionFile.IntervalCount() << " intervals, "
                << solutionFile.AntennaCount() << " antennas, "
                << solutionFile.ChannelCount() << " channels, "
                << solutionFile.PolarizationCount()
                << " polarizations in solution file.\n";

      std::unique_ptr<SolutionFlagFile> flagFile;

      std::unique_ptr<SolutionFile> newFile;
      std::string newFilename(std::string(argv[argi]) + "-tmp");
      bool hasFlagfile = argi + 1 < size_t(argc);
      if (hasFlagfile) {
        flagFile.reset(new SolutionFlagFile(argv[argi + 1]));
        if (solutionFile.IntervalCount() != flagFile->IntervalCount())
          throw std::runtime_error(
              "solutionFile.IntervalCount() != flagFile.IntervalCount()");
        if (solutionFile.AntennaCount() != flagFile->AntennaCount())
          throw std::runtime_error(
              "solutionFile.AntennaCount() != flagFile.AntennaCount()");
        if (solutionFile.ChannelCount() != flagFile->ChannelCount())
          throw std::runtime_error(
              "solutionFile.ChannelCount() != flagFile.ChannelCount()");

        newFile->SetIntervalCount(solutionFile.IntervalCount());
        newFile->SetAntennaCount(solutionFile.AntennaCount());
        newFile->SetChannelCount(solutionFile.ChannelCount());
        newFile->SetPolarizationCount(solutionFile.PolarizationCount());
        newFile->OpenForWriting(newFilename.c_str());
      }

      size_t totalSolutions = solutionFile.IntervalCount() *
                              solutionFile.AntennaCount() *
                              solutionFile.ChannelCount() * 4;
      size_t alreadyFlaggedCount = 0, flaggedInFlagFile = 0, flagsChanged = 0,
             totalFlags = 0;
      for (size_t interval = 0; interval != solutionFile.IntervalCount();
           ++interval) {
        for (size_t a = 0; a != solutionFile.AntennaCount(); ++a) {
          for (size_t ch = 0; ch != solutionFile.ChannelCount(); ++ch) {
            for (size_t p = 0; p != 4; ++p) {
              std::complex<double> val = solutionFile.ReadNextSolution();
              bool alreadyFlagged =
                  !std::isfinite(val.real()) || !std::isfinite(val.imag());
              bool isFlaggedInFile = (flagFile == 0)
                                         ? false
                                         : flagFile->IsFlagged(interval, a, ch);
              if (alreadyFlagged) ++alreadyFlaggedCount;
              if (isFlaggedInFile) {
                val = std::complex<double>(
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN());
                ++flaggedInFlagFile;
                if (!alreadyFlagged) ++flagsChanged;
              }
              if (alreadyFlagged || isFlaggedInFile) ++totalFlags;
              if (hasFlagfile) newFile->WriteSolution(val, interval, a, ch, p);
            }
          }
        }
      }
      std::cout << "Already flagged:     " << alreadyFlaggedCount << " ("
                << round(1000.0 * double(alreadyFlaggedCount) /
                         totalSolutions) /
                       10.0
                << "%)\n"
                << "Flagged by flagfile: " << flaggedInFlagFile << '\n'
                << "Flags changed:       " << flagsChanged << '\n'
                << "Total flags now:     " << totalFlags << '\n';

      if (hasFlagfile)
        boost::filesystem::rename(newFilename, std::string(argv[argi]));

      if (comparisonFilename != 0) {
        SolutionFile compSol1, compSol2;
        compSol1.OpenForReading(comparisonFilename);
        compSol2.OpenForReading(solutionFilename);

        aocommon::UVector<double> sumXX(compSol1.AntennaCount(), 0.0),
            sumYY(compSol1.AntennaCount(), 0.0), referencePhase[4];
        aocommon::UVector<size_t> countXX(compSol1.AntennaCount(), 0),
            countYY(compSol1.AntennaCount(), 0);
        for (size_t i = 0; i != 4; ++i)
          referencePhase[i].resize(compSol1.ChannelCount());

        for (size_t a = 0; a != compSol1.AntennaCount(); ++a) {
          for (size_t ch = 0; ch != compSol1.ChannelCount(); ++ch) {
            double vals[4];
            vals[0] = std::arg(compSol1.ReadNextSolution());
            compSol1.ReadNextSolution();
            compSol1.ReadNextSolution();
            vals[1] = std::arg(compSol1.ReadNextSolution());

            vals[2] = std::arg(compSol2.ReadNextSolution());
            compSol2.ReadNextSolution();
            compSol2.ReadNextSolution();
            vals[3] = std::arg(compSol2.ReadNextSolution());

            if (a == 0) {
              for (size_t i = 0; i != 4; ++i) referencePhase[i][ch] = vals[i];
            }

            double xxDist = phaseDist(vals[0] - referencePhase[0][ch],
                                      vals[2] - referencePhase[2][ch]),
                   yyDist = phaseDist(vals[1] - referencePhase[1][ch],
                                      vals[3] - referencePhase[3][ch]);
            if (std::isfinite(xxDist)) {
              sumXX[a] += xxDist * xxDist;
              countXX[a]++;
            }
            if (std::isfinite(yyDist)) {
              sumYY[a] += yyDist * yyDist;
              countYY[a]++;
            }
          }
        }

        aocommon::UVector<std::pair<double, size_t>> xxToAntenna(
            compSol1.AntennaCount()),
            yyToAntenna(compSol1.AntennaCount());
        for (size_t a = 0; a != compSol1.AntennaCount(); ++a) {
          xxToAntenna[a] =
              std::make_pair(sqrt(sumXX[a] / countXX[a]) * 180.0 / M_PI, a);
          yyToAntenna[a] =
              std::make_pair(sqrt(sumYY[a] / countYY[a]) * 180.0 / M_PI, a);
        }
        std::sort(xxToAntenna.rbegin(), xxToAntenna.rend());
        std::sort(yyToAntenna.rbegin(), yyToAntenna.rend());

        std::cout << "XX\tRMS\tYY\tRM\n";
        for (size_t a = 0; a != compSol1.AntennaCount(); ++a) {
          std::cout << xxToAntenna[a].second << '\t' << xxToAntenna[a].first
                    << '\t' << yyToAntenna[a].second << '\t'
                    << yyToAntenna[a].first << '\n';
        }
      }
    }
  }
}
