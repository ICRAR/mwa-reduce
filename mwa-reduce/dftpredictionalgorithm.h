#ifndef DFT_PREDICTION_ALGORITHM_H
#define DFT_PREDICTION_ALGORITHM_H

#include "imagebufferallocator.h"
#include "polarizationenum.h"
#include "uvector.h"
#include "matrix2x2.h"

#include "lofar/lbeamevaluator.h"
#include "banddata.h"

#include <vector>
#include <complex>

/**
 * Structure:
 * - PredictionImage: images[4] -- collects the model images.
 * - PredictionInput: components[nComponents] -- made from image, used as input for prediction.
 * - PredictionComponent: l, m, flux[nChannel x 4], antennaBeamValues[nAntenna] (these are updated per timestep)
 * - DFTAntennaInfo: beamValuesPerChannel[nChannel] of Matrix2x2
 */

class DFTAntennaInfo
{
public:
	const MC2x2& BeamValue(size_t channelIndex) const { return _beamValuesPerChannel[channelIndex]; }
	MC2x2& BeamValue(size_t channelIndex) { return _beamValuesPerChannel[channelIndex]; }
	
	std::vector<MC2x2>::iterator begin() { return _beamValuesPerChannel.begin(); }
	std::vector<MC2x2>::iterator end() { return _beamValuesPerChannel.end(); }
	size_t ChannelCount() const { return _beamValuesPerChannel.size(); }
	void InitializeChannelBuffers(size_t channelCount) { _beamValuesPerChannel.resize(channelCount); }
private:
	std::vector<MC2x2> _beamValuesPerChannel;
};

class DFTPredictionComponent
{
public:
	DFTPredictionComponent(double ra, double dec, double l, double m, double fluxLinear[4]) :
		_ra(ra), _dec(dec), _l(l), _m(m), _lmSqrt(sqrt(1.0 - l*l - m*m))
	{
		for(size_t p=0; p!=4; ++p) _flux[p] = fluxLinear[p];
	}
	double L() const { return _l; }
	double M() const { return _m; }
	double RA() const { return _ra; }
	double Dec() const { return _dec; }
	double LMSqrt() const { return _lmSqrt; }
	const DFTAntennaInfo& AntennaInfo(size_t antennaIndex) const { return _beamValuesPerAntenna[antennaIndex]; }
	DFTAntennaInfo& AntennaInfo(size_t antennaIndex) { return _beamValuesPerAntenna[antennaIndex]; }
	const MC2x2& LinearFlux() const { return _flux; }
	size_t AntennaCount() const { return _beamValuesPerAntenna.size(); }
	void InitializeBeamBuffers(size_t antennaCount, size_t channelCount)
	{
		_beamValuesPerAntenna.resize(antennaCount);
		for(DFTAntennaInfo& a : _beamValuesPerAntenna)
			a.InitializeChannelBuffers(channelCount);
	}
private:
	double _ra, _dec, _l, _m, _lmSqrt;
	MC2x2 _flux;
	std::vector<DFTAntennaInfo> _beamValuesPerAntenna;
};

class DFTPredictionInput
{
public:
	void AddComponent(const DFTPredictionComponent& component)
	{
		_components.push_back(component);
	}
	size_t ComponentCount() const { return _components.size(); }
	void InitializeBeamBuffers(size_t antennaCount, size_t channelCount) {
		for(DFTPredictionComponent& c : _components) c.InitializeBeamBuffers(antennaCount, channelCount);
	}
	std::vector<DFTPredictionComponent>::const_iterator begin() const { return _components.begin(); }
	std::vector<DFTPredictionComponent>::const_iterator end() const { return _components.end(); }
	std::vector<DFTPredictionComponent>::iterator begin() { return _components.begin(); }
	std::vector<DFTPredictionComponent>::iterator end() { return _components.end(); }
private:
	std::vector<DFTPredictionComponent> _components;
};

class DFTPredictionImage
{
public:
	DFTPredictionImage(size_t width, size_t height, ImageBufferAllocator<double>& allocator);
	
	void Add(PolarizationEnum polarization, const double* image);
	void Add(PolarizationEnum polarization, const double* real, const double* imaginary);
	
	void FindComponents(DFTPredictionInput& destination, double phaseCentreRA, double phaseCentreDec, double pixelSizeX, double pixelSizeY, double dl, double dm);
private:
	size_t _width, _height;
	ImageBufferAllocator<double>* _allocator;
	ImageBufferAllocator<double>::Ptr _images[4];
	std::vector<PolarizationEnum> _pols;
};

class DFTPredictionAlgorithm
{
public:
	DFTPredictionAlgorithm(DFTPredictionInput& input, const BandData& band) : _input(input), _band(band)
	{ }
	
	void Predict(MC2x2& dest, double u, double v, double w, size_t channelIndex, size_t a1, size_t a2);

	void UpdateBeam(LBeamEvaluator& beamEvaluator);
	
private:
	void predict(MC2x2& dest, double u, double v, double w, size_t channelIndex, size_t a1, size_t a2, DFTPredictionComponent& component);
	
	DFTPredictionInput& _input;
	BandData _band;
};

#endif
