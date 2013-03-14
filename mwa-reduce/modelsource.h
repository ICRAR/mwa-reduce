#ifndef MODELSOURCE_H
#define MODELSOURCE_H

#include <string>
#include <sstream>

#include "sourcesdf.h"
#include "sourcesdfwithsamples.h"
#include "radeccoord.h"

class ModelSource
{
	public:
		enum Type { PointSource };
		
		ModelSource() : _name(), _type(PointSource), _posRA(0.0), _posDec(0.0), _brightness(new SourceSDFWithSI<long double>()), _l(0.0), _m(0.0)
		{
		}
		
		ModelSource(const ModelSource &source) : _name(source._name), _type(source._type), _posRA(source._posRA), _posDec(source._posDec), _brightness(source._brightness->Copy())
		{
		}
		
		~ModelSource()
		{
			delete _brightness;
		}
		
		ModelSource& operator=(const ModelSource &source)
		{
			_name = source._name;
			_type = source._type;
			_posRA = source._posRA;
			_posDec = source._posDec;
			_brightness = source._brightness->Copy();
			return *this;
		}
		
		const std::string &Name() const { return _name; }
		enum Type Type() const { return _type; }
		long double PosRA() const { return _posRA; }
		long double PosDec() const { return _posDec; }
		const SourceSDF<long double> &Brightness() const { return *_brightness; }
		
		void SetName(const std::string &name) { _name = name; }
		void SetType(enum Type type) { _type = type; }
		void SetPosRA(long double posRA) { _posRA = posRA; }
		void SetPosDec(long double posDec) { _posDec = posDec; }
		SourceSDF<long double> &Brightness() { return *_brightness; }
		void SetBrightness(const SourceSDF<long double> &sdf) {
			delete _brightness; 
			_brightness = sdf.Copy();
		}
		
		long double L() const { return _l; }
		long double M() const { return _m; }
		void *UserData() const { return _userdata; }
		
		void SetL(long double l) { _l = l; }
		void SetM(long double m) { _m = m; }
		void SetUserData(void *userData) { _userdata = userData; }
		
		std::string ToStringLine() const {
			std::stringstream s;
			s << _name << " point " << RaDecCoord::RAToString(_posRA) << ' ' << RaDecCoord::DecToString(_posDec) << ' ' << _brightness->ToString();
			return s.str();
		}
	private:
		std::string _name;
		enum Type _type;
		long double _posRA, _posDec;
		SourceSDF<long double> *_brightness;
		long double _l, _m;
		void *_userdata;
};

#endif
