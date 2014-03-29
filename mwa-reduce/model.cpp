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

void Model::operator+=(const Model& rhs)
{
	if(Empty())
		(*this) = rhs;
	else {
		for(const_iterator i = rhs.begin(); i!=rhs.end(); ++i)
			add(*i);
	}
}

void Model::Optimize()
{
	Model copy(*this);
	_sources.clear();
	for(const_iterator i = copy.begin(); i!=copy.end(); ++i)
		addOptimized(*i);
}

void Model::add(const ModelSource& source)
{
	for(iterator i = begin(); i!=end(); ++i)
	{
		if(source.Peak().PosDec() == i->Peak().PosDec() && source.Peak().PosRA() == i->Peak().PosRA())
		{
			(*i) += source;
			return;
		}
	}
	_sources.push_back(source);
}

void Model::addOptimized(const ModelSource& source)
{
	for(iterator i = begin(); i!=end(); ++i)
	{
		if(source.Peak().PosDec() == i->Peak().PosDec() && source.Peak().PosRA() == i->Peak().PosRA())
		{
			/* merge */
			return;
		}
	}
	_sources.push_back(source);
}

void Model::combineMeasurements(const ModelSource& source)
{
	for(iterator i = begin(); i!=end(); ++i)
	{
		if(source.Peak().PosDec() == i->Peak().PosDec() && source.Peak().PosRA() == i->Peak().PosRA())
		{
			i->CombineMeasurements(source);
			return;
		}
	}
	throw std::runtime_error("Combining measurements while not same sources were measured!");
}

void Model::Save(const char* filename) const
{
	std::ofstream stream(filename);
	Save(stream);
}

void Model::Save(std::ostream& stream) const
{
	stream << "skymodel fileformat 1.1\n";
	for(const_iterator i=begin(); i!=end(); ++i)
	{
		stream << i->ToString();
	}
}

void Model::SortOnBrightness()
{
	std::sort(_sources.rbegin(), _sources.rend());
}
