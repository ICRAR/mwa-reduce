#ifndef MODELSOURCE_H
#define MODELSOURCE_H

#include <string>
#include <sstream>

#include "sourcesdf.h"
#include "sourcesdfwithsamples.h"
#include "radeccoord.h"
#include "spectralenergydistribution.h"

class ModelSource
{
	public:
		enum Type { PointSource };
		
		ModelSource() : _name(), _type(PointSource), _posRA(0.0), _posDec(0.0), _sed(), _l(0.0), _m(0.0)
		{
		}
		
		ModelSource(const ModelSource &source) : _name(source._name), _type(source._type), _posRA(source._posRA), _posDec(source._posDec), _sed(source._sed)
		{
		}
		
		~ModelSource()
		{
		}
		
		ModelSource& operator=(const ModelSource &source)
		{
			_name = source._name;
			_type = source._type;
			_posRA = source._posRA;
			_posDec = source._posDec;
			_sed = source._sed;
			return *this;
		}
		
		const std::string &Name() const { return _name; }
		enum Type Type() const { return _type; }
		long double PosRA() const { return _posRA; }
		long double PosDec() const { return _posDec; }
		const SpectralEnergyDistribution &SED() const { return _sed; }
		
		void SetName(const std::string &name) { _name = name; }
		void SetType(enum Type type) { _type = type; }
		void SetPosRA(long double posRA) { _posRA = posRA; }
		void SetPosDec(long double posDec) { _posDec = posDec; }
		SpectralEnergyDistribution &SED() { return _sed; }
		void SetSED(const SpectralEnergyDistribution &sed) {
			_sed = sed;
		}
		
		long double L() const { return _l; }
		long double M() const { return _m; }
		void *UserData() const { return _userdata; }
		
		void SetL(long double l) { _l = l; }
		void SetM(long double m) { _m = m; }
		void SetUserData(void *userData) { _userdata = userData; }
		
		std::string ToString() const
		{
			std::stringstream s;
			s << "source {\n  name \"" << _name << "\"\n"
				"  component {\n"
				"    type point\n"
				"    position " << RaDecCoord::RAToString(_posRA) << ' ' << RaDecCoord::DecToString(_posDec) << '\n' <<
				_sed.ToString() << "  }\n}\n";
			return s.str();
		}
		
		bool operator<(const ModelSource &rhs) const
		{
			return _sed < rhs._sed;
		}
	private:
		std::string _name;
		enum Type _type;
		long double _posRA, _posDec;
		SpectralEnergyDistribution _sed;
		long double _l, _m;
		void *_userdata;
};

#endif
