#ifndef MRC_CATALOGUE_H
#define MRC_CATALOGUE_H

#include <string>
#include <fstream>

#include "modelsource.h"

class MRCCatalogue
{
	public:
		MRCCatalogue(const std::string &filename);
		
		bool ReadNext(ModelSource &source);
	private:
		int digToVal(char dig)
		{
			if(dig == 32) return 0;
			else return dig-'0';
		}
		
		std::string _filename;
		std::ifstream _file;
};

#endif
