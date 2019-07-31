#ifndef MS_SELECTION
#define MS_SELECTION

#include <cstring>
#include <casacore/casa/Arrays/Vector.h>

class MSSelection
{
public:
	enum EvenOddSelection { AllTimesteps, EvenTimesteps, OddTimesteps };
	
	MSSelection() :
		_fieldId(0),
		_bandId(0),
		_startChannel(0), _endChannel(0),
		_startTimestep(0), _endTimestep(0),
		_minUVWInM(0.0), _maxUVWInM(0.0),
		_autoCorrelations(false),
		_evenOddSelection(AllTimesteps)
	{ }
	
	bool HasChannelRange() const { return _endChannel != 0; }
	bool HasInterval() const { return _endTimestep != 0; }
	bool HasMinUVWInM() const { return _minUVWInM != 0.0; }
	bool HasMaxUVWInM() const { return _maxUVWInM != 0.0; }
	
	size_t ChannelRangeStart() const { return _startChannel; }
	size_t ChannelRangeEnd() const { return _endChannel; }
	
	size_t IntervalStart() const { return _startTimestep; }
	size_t IntervalEnd() const { return _endTimestep; }
	
	size_t FieldId() const { return _fieldId; }
	
	double MinUVWInM() const { return _minUVWInM; }
	double MaxUVWInM() const { return _maxUVWInM; }
	
	bool IsSelected(size_t fieldId, size_t timestep, size_t antenna1, size_t antenna2, const casacore::Vector<double>& uvw) const
	{
		if(HasMinUVWInM() || HasMaxUVWInM())
		{
			double u = uvw(0), v = uvw(1), w = uvw(2);
			return IsSelected(fieldId, timestep, antenna1, antenna2, sqrt(u*u + v*v + w*w));
		}
		else {
			return IsSelected(fieldId, timestep, antenna1, antenna2, 0.0);
		}
	}
	
	bool IsSelected(size_t fieldId, size_t timestep, size_t antenna1, size_t antenna2, double uvwInMeters) const
	{
		if(fieldId != _fieldId)
			return false;
		else if(HasInterval() && (timestep < _startTimestep || timestep >= _endTimestep))
			return false;
		else if(!_autoCorrelations && (antenna1 == antenna2))
			return false;
		else if(HasMinUVWInM() && uvwInMeters < _minUVWInM)
			return false;
		else if(HasMaxUVWInM() && uvwInMeters > _maxUVWInM)
			return false;
		else if(_evenOddSelection!=AllTimesteps)
		{
			if(_evenOddSelection==EvenTimesteps && timestep%2!=0)
				return false;
			else if(_evenOddSelection==OddTimesteps && timestep%2!=1)
				return false;
		}
		return true;
	}
	
	bool IsFieldSelected(size_t fieldId) const
	{
		return fieldId == _fieldId;
	}
	
	bool IsTimeSelected(size_t timestep)
	{
		if(HasInterval() && (timestep < _startTimestep || timestep >= _endTimestep))
			return false;
		else if(_evenOddSelection!=AllTimesteps)
		{
			if(_evenOddSelection==EvenTimesteps && timestep%2!=0)
				return false;
			else if(_evenOddSelection==OddTimesteps && timestep%2!=1)
				return false;
		}
		return true;
	}
	
	void SetFieldId(size_t fieldId)
	{ 
		_fieldId = fieldId; 
	}
	void SetBandId(size_t bandId)
	{
		_bandId = bandId;
	}
	void SetChannelRange(size_t startChannel, size_t endChannel)
	{
		_startChannel = startChannel;
		_endChannel = endChannel;
	}
	void SetNoChannelRange()
	{
		_startChannel = 0;
		_endChannel = 0;
	}
	void SetInterval(size_t startTimestep, size_t endTimestep)
	{
		_startTimestep = startTimestep;
		_endTimestep = endTimestep;
	}
	void SetMinUVWInM(double minUVW)
	{
		_minUVWInM = minUVW;
	}
	void SetMaxUVWInM(double maxUVW)
	{
		_maxUVWInM = maxUVW;
	}
	void SetEvenOrOddTimesteps(EvenOddSelection evenOrOdd)
	{
		_evenOddSelection = evenOrOdd;
	}
	static MSSelection Everything() { return MSSelection(); }
private:
	size_t _fieldId, _bandId;
	size_t _startChannel, _endChannel;
	size_t _startTimestep, _endTimestep;
	double _minUVWInM, _maxUVWInM;
	bool _autoCorrelations;
	enum EvenOddSelection _evenOddSelection;
};

#endif
