#include <iostream>
#include <fstream>
#include <string>

#include "model/model.h"

#include <aocommon/fits/fitsreader.h>
#include <aocommon/units/angle.h>

#include <boost/tokenizer.hpp>

using namespace std;
using namespace aocommon;

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "Syntax: aegean2model [-fitsfreq <fitsfile>] [-use-peak] [-correct "
            "beam <maj> <min>] [-ncp] <aegean-file> <output-model> "
            "<src-prefix>\n";
  } else {
    std::string freqfitsfile;
    size_t argi = 1;
    bool usePeak = false, correctBeam = false, ncp = false;
    double beamMaj = 0.0, beamMin = 0.0, beamIntegral = 0.0;
    while (argv[argi][0] == '-') {
      std::string p(&argv[argi][1]);
      if (p == "fitsfreq") {
        ++argi;
        freqfitsfile = argv[argi];
      } else if (p == "use-peak") {
        usePeak = true;
      } else if (p == "correct-beam") {
        correctBeam = true;
        beamMaj = atof(argv[argi + 1]);
        beamMin = atof(argv[argi + 2]);
        // Using the FWHM formula for a Gaussian:
        long double sigmaMaj = beamMaj / (2.0L * sqrtl(2.0L * logl(2.0L)));
        long double sigmaMin = beamMin / (2.0L * sqrtl(2.0L * logl(2.0L)));
        beamIntegral = 2.0L * M_PI * sigmaMaj * sigmaMin;
        argi += 2;
      } else if (p == "ncp") {
        ncp = true;
      } else {
        throw std::runtime_error("Unknown parameter");
      }
      ++argi;
    }
    ifstream aegeanFile(argv[argi]);
    std::string outputFilename(argv[argi + 1]);
    string srcPrefix(argv[argi + 2]);
    double frequency = 150.0e6;
    if (!freqfitsfile.empty()) {
      FitsReader reader(freqfitsfile);
      frequency = reader.Frequency();
    }
    boost::char_delimiters_separator<char> sep(false, "", " \t");
    Model model;
    size_t srcIndex = 0;
    while (aegeanFile.good()) {
      string line;
      getline(aegeanFile, line);
      if (!line.empty() && line[0] != '#') {
        boost::tokenizer<boost::char_delimiters_separator<char> > tok(line,
                                                                      sep);
        boost::tokenizer<boost::char_delimiters_separator<char> >::iterator
            beg = tok.begin();
        for (size_t i = 0; i != 5; ++i) ++beg;
        double ra = atof(beg->c_str());
        if (ncp) ra += 180.0;
        ++beg;
        ++beg;
        double dec = atof(beg->c_str());
        ++beg;
        ++beg;
        double flux;
        if (usePeak) {
          flux = atof(beg->c_str());
          ++beg;
          ++beg;
        } else {
          ++beg;
          ++beg;
          flux = atof(beg->c_str());
        }
        ++beg;
        ++beg;
        double major = atof(beg->c_str());
        ++beg;
        ++beg;
        double minor = atof(beg->c_str());
        ++beg;
        ++beg;
        double pa = atof(beg->c_str());
        ModelComponent component;
        component.SetPosRA(ra * M_PI / 180.0);
        component.SetPosDec(dec * M_PI / 180.0);
        component.SetSED(MeasuredSED(flux, frequency));
        if (correctBeam) {
          long double sigmaMaj = major / (2.0L * sqrtl(2.0L * logl(2.0L)));
          long double sigmaMin = minor / (2.0L * sqrtl(2.0L * logl(2.0L)));
          long double sourceIntegral = 2.0L * M_PI * sigmaMaj * sigmaMin;
          long double correction = sqrt(beamIntegral / sourceIntegral);
          long double asec = M_PI / 180.0 / 60.0 / 60.0;
          std::cout << "Correcting " << flux << " Jy source from "
                    << units::Angle::ToNiceString(major * asec) << " x "
                    << units::Angle::ToNiceString(minor * asec) << " to "
                    << units::Angle::ToNiceString(major * correction * asec)
                    << " x "
                    << units::Angle::ToNiceString(minor * correction * asec)
                    << "\n";
          major *= correction;
          minor *= correction;
        }
        if (usePeak) {
          component.SetType(ModelComponent::GaussianSource);
          component.SetMajorAxis(major * M_PI / (180.0 * 60.0 * 60.0));
          component.SetMinorAxis(minor * M_PI / (180.0 * 60.0 * 60.0));
          component.SetPositionAngle(pa * M_PI / 180.0);
        }
        ModelSource source;
        ostringstream srcName;
        ++srcIndex;
        srcName << srcPrefix << srcIndex;
        source.SetName(srcName.str());
        source.AddComponent(component);
        model.AddSource(source);
      }
    }
    model.Save(outputFilename);
  }
}
