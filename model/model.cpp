#include "model.h"

#include "modelparser.h"

#include "../units/radeccoord.h"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <fstream>
#include <stdexcept>

size_t Model::npos = std::numeric_limits<size_t>::max();

void Model::read(const char* filename)
{
	ModelParser parser;
	std::ifstream stream(filename);
	if(!stream.good())
		throw std::runtime_error(std::string("Could not open model ") + filename);
	parser.Parse(*this, stream);
}

void Model::operator+=(const Model& rhs)
{
	if(Empty())
		(*this) = rhs;
	else {
		for(const_iterator i = rhs.begin(); i!=rhs.end(); ++i)
			add(*i);
		for(std::map<std::string,ModelCluster>::const_iterator c = rhs._clusters.begin();
				c!=rhs._clusters.end(); ++c)
			_clusters.insert(*c);
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
	if(source.ComponentCount()!=0)
	{
		for(iterator i = begin(); i!=end(); ++i)
		{
			if(i->ComponentCount()!=0 && source.Peak().PosDec() == i->Peak().PosDec() && source.Peak().PosRA() == i->Peak().PosRA())
			{
				(*i) += source;
				return;
			}
		}
	}
	_sources.push_back(source);
}

void Model::addOptimized(const ModelSource& source)
{
	if(source.ComponentCount()!=0)
	{
		for(iterator i = begin(); i!=end(); ++i)
		{
			if(i->ComponentCount()!=0 && source.Peak().PosDec() == i->Peak().PosDec() && source.Peak().PosRA() == i->Peak().PosRA())
			{
				/* merge */
				return;
			}
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
