#include <iostream>
#include <fstream>
#include <string>

#include "model.h"
#include "fitsreader.h"

#include <boost/tokenizer.hpp>

using namespace std;

int main(int argc, char* argv[])
{
	if(argc < 4) {
		cout << "Syntax: aegean2model [-fitsfreq <fitsfile>] <aegean-file> <output-model> <src-prefix>\n";
	}
	else {
		std::string freqfitsfile;
		size_t argi = 1;
		while(argv[argi][0] == '-')
		{
			std::string p(&argv[argi][1]);
			if(p == "fitsfreq")
			{
				++argi;
				freqfitsfile = argv[argi];
			}
			else {
				throw std::runtime_error("Unknown parameter");
			}
			++argi;
		}
		ifstream aegeanFile(argv[argi]);
		std::string outputFilename(argv[argi+1]);
		string srcPrefix(argv[argi+2]);
		double frequency = 150.0e6;
		if(!freqfitsfile.empty())
		{
			FitsReader reader(freqfitsfile);
			frequency = reader.Frequency();
		}
		boost::char_delimiters_separator<char> sep(false, "", " \t");
		Model model;
		size_t srcIndex = 0;
		while(aegeanFile.good())
		{
			string line;
			getline(aegeanFile, line);
			if(!line.empty() && line[0] != '#')
			{
				boost::tokenizer<boost::char_delimiters_separator<char> > tok(line, sep);
				boost::tokenizer<boost::char_delimiters_separator<char> >::iterator beg=tok.begin();
				for(size_t i=0; i!=5; ++i) ++beg;
				double ra = atof(beg->c_str());
				++beg; ++beg;
				double dec = atof(beg->c_str());
				++beg; ++beg; ++beg; ++beg;
				double flux = atof(beg->c_str());
				ModelComponent component;
				component.SetPosRA(ra * M_PI / 180.0);
				component.SetPosDec(dec * M_PI / 180.0);
				component.SetSED(SpectralEnergyDistribution(flux, frequency));
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
