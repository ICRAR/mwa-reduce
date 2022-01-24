#include "model.h"

#include "bbsmodel.h"
#include "modelparser.h"

#include <aocommon/radeccoord.h>

#include <algorithm>
#include <cstdlib>
#include <string>
#include <fstream>
#include <stdexcept>

size_t Model::npos = std::numeric_limits<size_t>::max();

void Model::read(const char* filename) {
  std::ifstream stream(filename);
  if (!stream.good())
    throw std::runtime_error(std::string("Could not open model ") + filename);
  if (ModelParser::IsInModelFormat(stream)) {
    stream.seekg(0);
    ModelParser().Parse(*this, stream);
  } else {
    stream.close();
    *this = BBSModel::Read(filename);
  }
}

void Model::operator+=(const Model& rhs) {
  if (Empty())
    (*this) = rhs;
  else {
    for (const ModelSource& s : rhs) add(s);
    for (const std::pair<const std::string, ModelCluster>& c : rhs._clusters)
      _clusters.emplace(c);
  }
}

void Model::add(const ModelSource& source) { _sources.push_back(source); }

void Model::combineMeasurements(const ModelSource& source) {
  for (iterator i = begin(); i != end(); ++i) {
    if (source.Peak().PosDec() == i->Peak().PosDec() &&
        source.Peak().PosRA() == i->Peak().PosRA()) {
      i->CombineMeasurements(source);
      return;
    }
  }
  throw std::runtime_error(
      "Combining measurements while not same sources were measured!");
}

void Model::Save(const char* filename) const {
  std::ofstream stream(filename);
  Save(stream);
}

void Model::Save(std::ostream& stream) const {
  stream << "skymodel fileformat 1.1\n";
  for (const_iterator i = begin(); i != end(); ++i) {
    stream << i->ToString();
  }
}
