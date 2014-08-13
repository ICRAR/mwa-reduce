#include "solutionflagfile.h"
#include "solutionfile.h"

#include <iostream>
#include <limits>
#include <boost/filesystem/operations.hpp>

int main(int argc, char **argv)
{
  if(argc < 3)
    {
      std::cout << "Usage: flagsolutions <solutions.bin> <solutions-flag-file.txt>\n"
				"Will flag solutions in the solution file as specified by the flag file.\n";
    } else {
		size_t argi = 1;
		while(argv[argi][0] == '-')
		{
			//std::string p(&argv[argi][1]);
			throw std::runtime_error("Bad option");
			++argi;
		}
		SolutionFile solutionFile;
		solutionFile.OpenForReading(argv[argi]);
		
		std::cout << solutionFile.IntervalCount() << " intervals, " << solutionFile.AntennaCount() << " antennas, " << solutionFile.ChannelCount() << " channels, " << solutionFile.PolarizationCount() << " polarizations in solution file.\n";
		
		SolutionFlagFile flagFile(argv[argi+1]);
		
		if(solutionFile.IntervalCount() != flagFile.IntervalCount())
			throw std::runtime_error("solutionFile.IntervalCount() != flagFile.IntervalCount()");
		if(solutionFile.AntennaCount() != flagFile.AntennaCount())
			throw std::runtime_error("solutionFile.AntennaCount() != flagFile.AntennaCount()");
		if(solutionFile.ChannelCount() != flagFile.ChannelCount())
			throw std::runtime_error("solutionFile.ChannelCount() != flagFile.ChannelCount()");
		
		SolutionFile newFile;
		std::string newFilename(std::string(argv[argi]) + "-tmp");
		newFile.SetIntervalCount(solutionFile.IntervalCount());
		newFile.SetAntennaCount(solutionFile.AntennaCount());
		newFile.SetChannelCount(solutionFile.ChannelCount());
		newFile.SetPolarizationCount(solutionFile.PolarizationCount());
		newFile.OpenForWriting(newFilename.c_str());
		
		size_t alreadyFlaggedCount = 0, flaggedInFlagFile = 0, flagsChanged = 0, totalFlags = 0;
		for(size_t interval=0; interval!=solutionFile.IntervalCount(); ++interval) {
			for(size_t a = 0; a!=solutionFile.AntennaCount(); ++a) {
				for(size_t ch = 0; ch!=solutionFile.ChannelCount(); ++ch) {
					for(size_t p = 0; p!=4; ++p) {
						std::complex<double> val = solutionFile.ReadNextSolution();
						bool alreadyFlagged = !std::isfinite(val.real()) || !std::isfinite(val.imag());
						if(alreadyFlagged)
							++alreadyFlaggedCount;
						if(flagFile.IsFlagged(interval, a, ch))
						{
							val = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
							++flaggedInFlagFile;
							if(!alreadyFlagged)
								++flagsChanged;
						}
						if(alreadyFlagged || (flagFile.IsFlagged(interval, a, ch)))
							++totalFlags;
						newFile.WriteSolution(val, interval, a, ch, p);
					}
				}
			}
		}
		std::cout <<
			"Already flagged:     " << alreadyFlaggedCount << '\n' <<
			"Flagged by flagfile: " << flaggedInFlagFile << '\n' <<
			"Flags changed:       " << flagsChanged << '\n' <<
			"Total flags now:     " << totalFlags << '\n';
		boost::filesystem::rename(newFilename, std::string(argv[argi]));
	}
}
