#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include <stdint.h>

struct MWAFHeader {
  char fileIdentifier[4];
  uint16_t versionMinor, versionMajor;
  uint32_t timestepCount, antennaCount, channelCount, polarizationCount,
      gpuBoxIndex;
  uint64_t gpsTime;
  char baselineSelection;
};

/**
 * Convert an array of bits to an array of bytes containing '1' or '0'.
 * Each bit thus represents a boolean.
 * @param output Destination array. Should have room to store "count" bytes.
 * @param input Input array, should contain ceil(count/8) bytes.
 * @param count Number of bits in the input array.
 */
void unpackBitsToBytes(unsigned char* output, const unsigned char* input,
                       size_t count) {
  size_t countWhole = count / 8;
  for (size_t i = 0; i < countWhole; ++i) {
    char val = *input;
    for (size_t b = 0; b != 8; ++b) {
      *output = ((val & 1) == 0) ? 0 : 1;
      val >>= 1;
      ++output;
    }
    ++input;
  }
  size_t remainder = count % 8;
  char val = *input;
  for (size_t b = 0; b != remainder; ++b) {
    *output = ((val & 1) == 0) ? 0 : 1;
    val >>= 1;
    ++output;
  }
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Expecting filename: mwafinfo <filename>\n";
    exit(1);
  }

  std::ifstream file(argv[1]);
  MWAFHeader header;
  file.read(reinterpret_cast<char*>(&header), sizeof(header));
  const char* ident = header.fileIdentifier;
  if (ident[0] != 'M' || ident[1] != 'W' || ident[2] != 'A' ||
      ident[3] != 'F') {
    std::cerr << "Error: file does not appear to be a MWA flag file "
                 "(identifier missing at start of file)\n";
    exit(1);
  }

  if (header.baselineSelection != 'B') {
    std::cerr << "Error in file: baseline selection not set to 'both'\n";
    exit(1);
  }

  std::cout << "Version:        " << header.versionMajor << "."
            << header.versionMinor
            << "\n"
               "Timestep count: "
            << header.timestepCount
            << "\n"
               "Antenna count:  "
            << header.antennaCount
            << "\n"
               "Channel count:  "
            << header.channelCount
            << "\n"
               "Polarizations:  "
            << header.polarizationCount
            << "\n"
               "GPU box index:  "
            << header.gpuBoxIndex
            << "\n"
               "GPS time:       "
            << header.gpsTime
            << "\n"
               "Baselines       ALL\n";

  if (header.polarizationCount != 1) {
    std::cerr << "Error: I only know how to read single-polarization files.\n";
    exit(1);
  }

  // Calculate nr bytes per "row"
  // (A row is a scan in one baseline, i.e. nr channels consecutive samples)
  size_t stride = (header.channelCount + 7) / 8;
  std::vector<unsigned char> readBuffer(stride),
      unpackBuffer(header.channelCount);

  std::vector<size_t> antennaFlagCount(header.antennaCount, 0),
      channelFlagCount(header.channelCount, 0),
      timestepFlagCount(header.timestepCount, 0);
  size_t totalCount = 0, totalFlaggedCount = 0;

  for (size_t timestep = 0; timestep != header.timestepCount; ++timestep) {
    for (size_t ant1 = 0; ant1 != header.antennaCount; ++ant1) {
      for (size_t ant2 = ant1; ant2 != header.antennaCount; ++ant2) {
        file.read(reinterpret_cast<char*>(&readBuffer[0]), stride);
        if (!file.good()) {
          std::cerr << "Error reading file.\n";
          exit(1);
        }

        unpackBitsToBytes(&unpackBuffer[0], &readBuffer[0],
                          header.channelCount);

        if (ant1 != ant2) {
          for (size_t ch = 0; ch != header.channelCount; ++ch) {
            ++totalCount;
            if (unpackBuffer[ch] != 0) {
              ++totalFlaggedCount;
              ++antennaFlagCount[ant1];
              ++antennaFlagCount[ant2];
              ++channelFlagCount[ch];
              ++timestepFlagCount[timestep];
            }
          }
        }
      }
    }
  }

  std::vector<size_t>::const_iterator
      bestAntenna =
          std::min_element(antennaFlagCount.begin(), antennaFlagCount.end()),
      bestChannel =
          std::min_element(channelFlagCount.begin(), channelFlagCount.end()),
      bestTimestep =
          std::min_element(timestepFlagCount.begin(), timestepFlagCount.end()),
      worstAntenna =
          std::max_element(antennaFlagCount.begin(), antennaFlagCount.end()),
      worstChannel =
          std::max_element(channelFlagCount.begin(), channelFlagCount.end()),
      worstTimestep =
          std::max_element(timestepFlagCount.begin(), timestepFlagCount.end());

  // These statistics are not extremely useful, because they do not distinguish
  // between actual RFI or subband/dc channels etc., but useful for validating
  // file format.
  std::cout << "Total flags:    " << totalFlaggedCount << " ("
            << round(1000.0 * totalFlaggedCount / totalCount) / 10.0
            << "%)\n"
               "Best antenna:   "
            << (bestAntenna - antennaFlagCount.begin()) << " (" << *bestAntenna
            << " flags)\n"
               "Best channel:   "
            << (bestChannel - channelFlagCount.begin()) << " (" << *bestChannel
            << " flags)\n"
               "Best timestep:  "
            << (bestTimestep - timestepFlagCount.begin()) << " ("
            << *bestTimestep
            << " flags)\n"
               "Worst antenna:  "
            << (worstAntenna - antennaFlagCount.begin()) << " ("
            << *worstAntenna
            << " flags)\n"
               "Worst channel:  "
            << (worstChannel - channelFlagCount.begin()) << " ("
            << *worstChannel
            << " flags)\n"
               "Worst timestep: "
            << (worstTimestep - timestepFlagCount.begin()) << " ("
            << *worstTimestep << " flags)\n";

  std::streampos filePos = file.tellg();
  file.seekg(0, std::ios::end);
  if (filePos != file.tellg()) {
    std::cerr << "ERROR: after reading all data, file position was not at the "
                 "end of the file!\n";
  }
}
