#ifndef MODEL_COMPONENT_H
#define MODEL_COMPONENT_H

#include "spectralenergydistribution.h"
#include "measuredsed.h"

#include "../units/imagecoordinates.h"
#include "../units/radeccoord.h"

class ModelComponent
{
	public:
		enum Type { PointSource, GaussianSource };
		
		ModelComponent() : _type(PointSource), _posRA(0.0), _posDec(0.0), _sed(), _l(0.0), _m(0.0), _positionAngle(0.0), _majorAxis(0.0), _minorAxis(0.0), _userdata(0)
		{
		}
		
		ModelComponent(const ModelComponent& source) :
		 _type(source._type),
		 _posRA(source._posRA), _posDec(source._posDec),
		 _sed((source._sed == 0) ? 0 : source._sed->Clone()),
		 _l(source._l), _m(source._m),
		 _positionAngle(source._positionAngle), _majorAxis(source._majorAxis), _minorAxis(source._minorAxis),
			_userdata(source._userdata)
		{
		}
		
		ModelComponent& operator=(const ModelComponent& source)
		{
			_type=source._type;
			_posRA=source._posRA; _posDec=source._posDec;
			if(source._sed == 0)
				_sed.reset();
			else
				_sed.reset(source._sed->Clone());
			_l=source._l; _m=source._m;
			_positionAngle=source._positionAngle; _majorAxis=source._majorAxis; _minorAxis=source._minorAxis;
			_userdata=source._userdata;
			return *this;
		}
		
		enum Type Type() const { return _type; }
		long double PosRA() const { return _posRA; }
		long double PosDec() const { return _posDec; }
		bool HasSED() const { return _sed != 0; }
		SpectralEnergyDistribution &SED() { return *_sed; }
		const SpectralEnergyDistribution &SED() const { return *_sed; }
		bool HasMeasuredSED() const { return dynamic_cast<MeasuredSED*>(_sed.get())!=0; }
		MeasuredSED& MSED() { return static_cast<MeasuredSED&>(*_sed); }
		const MeasuredSED& MSED() const { return static_cast<const MeasuredSED&>(*_sed); }
		long double L() const { return _l; }
		long double M() const { return _m; }
		/** PA in radians. */
		long double PositionAngle() const { return _positionAngle; }
		/** Major Gaussian axis in radians. */
		long double MajorAxis() const { return _majorAxis; }
		/** Minor Gaussian axis in radians. */
		long double MinorAxis() const { return _minorAxis; }
		
		void* UserData() const { return _userdata; }
		
		void SetType(enum Type type) { _type = type; }
		void SetPosRA(long double posRA) { _posRA = posRA; }
		void SetPosDec(long double posDec) { _posDec = posDec; }
		void SetSED(const SpectralEnergyDistribution& sed) {
			_sed.reset(sed.Clone());
		}
		void SetL(long double l) { _l = l; }
		void SetM(long double m) { _m = m; }
		void SetPositionAngle(long double pa) { _positionAngle = pa; }
		void SetMajorAxis(long double majorAxis) { _majorAxis = majorAxis; }
		void SetMinorAxis(long double minorAxis) { _minorAxis = minorAxis; }
		void SetUserData(void *userData) { _userdata = userData; }
		
		std::string ToString() const
		{
			std::stringstream s;
			s.precision(15);
			s << "  component {\n";
			switch(_type)
			{
				case PointSource:
					s << "    type point";
					break;
				case GaussianSource:
					s << "    type gaussian\n    shape "
						<< _majorAxis*(60.0*60.0*180.0/M_PI) << ' '
						<< _minorAxis*(60.0*60.0*180.0/M_PI) << ' '
						<< _positionAngle*(180.0/M_PI);
					break;
			}
				s << "\n    position "
					<< RaDecCoord::RAToString(_posRA) << ' '
					<< RaDecCoord::DecToString(_posDec) << '\n';
			if(_sed)
				s << _sed->ToString();
			s << "  }\n";
			return s.str();
		}
		
		bool HasValidMeasurement() const { return HasMeasuredSED() && MSED().HasValidMeasurement(); }
		
		bool operator<(const ModelComponent& rhs) const
		{
			return (*_sed) < (*rhs._sed);
		}
		
		void operator*=(double factor)
		{
			(*_sed) *= factor;
		}
	private:
		enum Type _type;
		long double _posRA, _posDec;
		std::unique_ptr<SpectralEnergyDistribution> _sed;
		long double _l, _m;
		long double _positionAngle, _majorAxis, _minorAxis;
		void *_userdata;
};

#endif
