#include "model.h"
#include "radeccoord.h"

#include "parser/modelparser.h"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <fstream>
#include <stdexcept>

Model::Model(const char *filename)
{
	ModelParser parser;
	std::ifstream stream(filename);
	if(!stream.good())
		throw std::runtime_error("Could not open model");
	parser.Parse(*this, stream);
}

void Model::Optimize()
{
	Model copy(*this);
	_sources.clear();
	for(const_iterator i = copy.begin(); i!=copy.end(); ++i)
		addOptimized(*i);
}

void Model::addOptimized(const ModelSource& source)
{
	for(iterator i = begin(); i!=end(); ++i)
	{
		if(source.PosDec() == i->PosDec() && source.PosRA() == i->PosRA())
		{
			/* merge */
			return;
		}
	}
	_sources.push_back(source);
}

void Model::Save(const char* filename)
{
	std::ofstream file(filename);
	file << "skymodel fileformat 1.0\n";
	for(const_iterator i=begin(); i!=end(); ++i)
	{
		file << i->ToString();
	}
}

void Model::SortOnBrightness()
{
	std::sort(_sources.rbegin(), _sources.rend());
}
