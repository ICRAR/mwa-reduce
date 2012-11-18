#ifndef MODEL_H
#define MODEL_H

#include <vector>

#include "modelsource.h"

class Model
{
	public:
		typedef std::vector<ModelSource>::iterator iterator;
		typedef std::vector<ModelSource>::const_iterator const_iterator;
		
		Model()
		{
		}
		
		Model(const char *filename);
		
		size_t SourceCount() { return _sources.size(); }
		ModelSource &Source(size_t index) { return _sources[index]; }
		const ModelSource &Source(size_t index) const { return _sources[index]; }
		
		const_iterator begin() const { return _sources.begin(); }
		const_iterator end() const { return _sources.end(); }
		
		iterator begin() { return _sources.begin(); }
		iterator end() { return _sources.end(); }
		
	private:
		std::vector<ModelSource> _sources;
		
		static bool isCommentSymbol(char c) { return c=='#'; }
		static bool isDelimiter(char c) { return c==' ' || c=='\t' || c=='\r' || c=='\n';	}
};

#endif
