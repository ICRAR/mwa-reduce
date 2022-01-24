#ifndef MS_SELECTION
#define MS_SELECTION

#include <cstring>
#include <limits>

#include <casacore/casa/Arrays/Vector.h>

#include <vector>

class MSSelection {
 public:
  enum EvenOddSelection { AllTimesteps, EvenTimesteps, OddTimesteps };

  const static size_t ALL_FIELDS = std::numeric_limits<size_t>::max();

  MSSelection()
      : _fieldIds{0},
        _bandId(0),
        _startChannel(0),
        _endChannel(0),
        _startTimestep(0),
        _endTimestep(0),
        _minUVWInM(0.0),
        _maxUVWInM(0.0),
        _autoCorrelations(false),
        _evenOddSelection(AllTimesteps) {}

  bool HasChannelRange() const { return _endChannel != 0; }
  bool HasInterval() const { return _endTimestep != 0; }
  bool HasMinUVWInM() const { return _minUVWInM != 0.0; }
  bool HasMaxUVWInM() const { return _maxUVWInM != 0.0; }

  size_t BandId() const { return _bandId; }

  size_t ChannelRangeStart() const { return _startChannel; }
  size_t ChannelRangeEnd() const { return _endChannel; }

  size_t IntervalStart() const { return _startTimestep; }
  size_t IntervalEnd() const { return _endTimestep; }

  const std::vector<size_t>& FieldIds() const { return _fieldIds; }

  double MinUVWInM() const { return _minUVWInM; }
  double MaxUVWInM() const { return _maxUVWInM; }

  bool IsSelected(size_t fieldId, size_t timestep, size_t antenna1,
                  size_t antenna2, const casacore::Vector<double>& uvw) const {
    if (HasMinUVWInM() || HasMaxUVWInM()) {
      double u = uvw(0), v = uvw(1), w = uvw(2);
      return IsSelected(fieldId, timestep, antenna1, antenna2,
                        std::sqrt(u * u + v * v + w * w));
    } else {
      return IsSelected(fieldId, timestep, antenna1, antenna2, 0.0);
    }
  }

  bool IsSelected(size_t fieldId, size_t timestep, size_t antenna1,
                  size_t antenna2, double uvwInMeters) const {
    if (!IsFieldSelected(fieldId))
      return false;
    else if (HasInterval() &&
             (timestep < _startTimestep || timestep >= _endTimestep))
      return false;
    else if (!_autoCorrelations && (antenna1 == antenna2))
      return false;
    else if (HasMinUVWInM() && uvwInMeters < _minUVWInM)
      return false;
    else if (HasMaxUVWInM() && uvwInMeters > _maxUVWInM)
      return false;
    else if (_evenOddSelection != AllTimesteps) {
      if (_evenOddSelection == EvenTimesteps && timestep % 2 != 0)
        return false;
      else if (_evenOddSelection == OddTimesteps && timestep % 2 != 1)
        return false;
    }
    return true;
  }

  bool IsFieldSelected(size_t fieldId) const {
    return std::find(_fieldIds.begin(), _fieldIds.end(), fieldId) !=
               _fieldIds.end() ||
           _fieldIds[0] == ALL_FIELDS;
  }

  bool IsTimeSelected(size_t timestep) {
    if (HasInterval() &&
        (timestep < _startTimestep || timestep >= _endTimestep))
      return false;
    else if (_evenOddSelection != AllTimesteps) {
      if (_evenOddSelection == EvenTimesteps && timestep % 2 != 0)
        return false;
      else if (_evenOddSelection == OddTimesteps && timestep % 2 != 1)
        return false;
    }
    return true;
  }

  void SetFieldIds(const std::vector<size_t>& fieldIds) {
    _fieldIds = fieldIds;
  }
  void SetBandId(size_t bandId) { _bandId = bandId; }
  void SetChannelRange(size_t startChannel, size_t endChannel) {
    _startChannel = startChannel;
    _endChannel = endChannel;
  }
  void SetNoChannelRange() {
    _startChannel = 0;
    _endChannel = 0;
  }
  void SetInterval(size_t startTimestep, size_t endTimestep) {
    _startTimestep = startTimestep;
    _endTimestep = endTimestep;
  }
  void SetMinUVWInM(double minUVW) { _minUVWInM = minUVW; }
  void SetMaxUVWInM(double maxUVW) { _maxUVWInM = maxUVW; }
  void SetEvenOrOddTimesteps(EvenOddSelection evenOrOdd) {
    _evenOddSelection = evenOrOdd;
  }
  EvenOddSelection EvenOrOddTimesteps() const { return _evenOddSelection; }
  static MSSelection Everything() { return MSSelection(); }

  void Serialize(class SerialOStream& stream) const;
  void Unserialize(class SerialIStream& stream);

 private:
  std::vector<size_t> _fieldIds;
  size_t _bandId;
  size_t _startChannel, _endChannel;
  size_t _startTimestep, _endTimestep;
  double _minUVWInM, _maxUVWInM;
  bool _autoCorrelations;
  enum EvenOddSelection _evenOddSelection;
};

#endif
