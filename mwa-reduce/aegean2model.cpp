#include <iostream>
#include <fstream>
#include <string>

#include "model.h"

#include <boost/tokenizer.hpp>

using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 4) {
		cout << "Syntax: aegean2model <aegean-file> <output-model> <src-prefix>\n";
	}
	else {
		Model model;
		ifstream aegeanFile(argv[1]);
		string srcPrefix(argv[3]);
		boost::char_delimiters_separator<char> sep(false, "", " \t");
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
				component.SetSED(SpectralEnergyDistribution(flux, 150000000.0));
				ModelSource source;
				ostringstream srcName;
				++srcIndex;
				srcName << srcPrefix << srcIndex;
				source.SetName(srcName.str());
				source.AddComponent(component);
				model.AddSource(source);
			}
		}
		model.Save(argv[2]);
	}
}
