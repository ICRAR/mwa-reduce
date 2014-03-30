#ifndef FLUX_SPECTRUM_ACCUMULATOR_H
#define FLUX_SPECTRUM_ACCUMULATOR_H

#include "fluxaccumulator.h"
#include "banddata.h"
#include "imagecoordinates.h"

#include <vector>

class FluxSpectrumAccumulator
{
public:
	FluxSpectrumAccumulator(const ModelComponent& component, const BandData* bandData, size_t channelBlocks, double phaseCentreRA, double phaseCentreDec) :
		_accumulators(bandData->ChannelCount()),
		_component(component),
		_bandData(bandData),
		_channelBlocks(channelBlocks)
	{
		long double l, m;
		ImageCoordinates::RaDecToLM<long double>(component.PosRA(), component.PosDec(), phaseCentreRA, phaseCentreDec, l, m);
		for(size_t ch=0; ch!=bandData->ChannelCount(); ++ch)
			_accumulators[ch] = new FluxAccumulator(l, m, bandData->ChannelWavelength(ch));
	}
	
	void UpdateBeam(BeamEvaluator& evaluator, double* ionG, double* ionDL, double* ionDM)
	{
		BeamEvaluator::PrecalcPosInfo posInfo;
		evaluator.PrecalculatePositionInfo(posInfo, _component.PosRA(), _component.PosDec());
		std::complex<double> beamGains[4];
		for(size_t ch=0; ch!=_bandData->ChannelCount(); ++ch)
		{
			size_t channelBlockIndex = ch * _channelBlocks / _bandData->ChannelCount();
			const double frequency = _bandData->ChannelFrequency(ch);
			evaluator.EvaluateApparentToAbsGain(posInfo, frequency, beamGains);
			_accumulators[ch]->UpdateBeam(beamGains, ionG[channelBlockIndex], ionDL[channelBlockIndex], ionDM[channelBlockIndex]);
		}
	}
	
	void Accumulate(std::complex<double>* data, const double weight, size_t ch, double u, double v, double w)
	{
		_accumulators[ch]->Add(data, weight, u, v, w);
	}
private:
	std::vector<FluxAccumulator*> _accumulators;
	ModelComponent _component;
	const BandData* _bandData;
	size_t _channelBlocks;
};

#endif
