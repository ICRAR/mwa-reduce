#ifndef FLUX_SPECTRUM_ACCUMULATOR_H
#define FLUX_SPECTRUM_ACCUMULATOR_H

#include "fluxaccumulator.h"
#include "banddata.h"
#include "imagecoordinates.h"

#include <vector>

class FluxSpectrumAccumulator : public Serializable
{
public:
	/**
	 * Constructor for accumulating visibilities
	 */
	FluxSpectrumAccumulator(const ModelComponent& component, const BandData* bandData, size_t channelBlocks, double phaseCentreRA, double phaseCentreDec) :
		_accumulators(bandData->ChannelCount()),
		_component(component),
		_bandData(bandData),
		_channelBlocks(channelBlocks)
	{
		long double l, m;
		ImageCoordinates::RaDecToLM<long double>(component.PosRA(), component.PosDec(), phaseCentreRA, phaseCentreDec, l, m);
		double stokesValues[4];
		std::complex<double> linearValues[4];
		for(size_t ch=0; ch!=bandData->ChannelCount(); ++ch)
		{
			for(size_t p=0; p!=4; ++p)
				stokesValues[p] = component.SED().FluxAtFrequency(bandData->ChannelFrequency(ch), Polarization::IndexToStokes(p));
			Polarization::StokesToLinear(stokesValues, linearValues);
			_accumulators[ch] = new FluxAccumulator(l, m, bandData->ChannelWavelength(ch), linearValues);
		}
	}
	
	/**
	 * Constructor for accumulating from files
	 */
	FluxSpectrumAccumulator(const ModelComponent& component, const BandData* bandData) :
		_accumulators(bandData->ChannelCount()),
		_component(component),
		_bandData(bandData),
		_channelBlocks(0)
	{
		for(size_t ch=0; ch!=bandData->ChannelCount(); ++ch)
		{
			_accumulators[ch] = new FluxAccumulator(bandData->ChannelWavelength(ch));
		}
	}
	
	~FluxSpectrumAccumulator()
	{
		for(size_t ch=0; ch!=_accumulators.size(); ++ch)
			delete _accumulators[ch];
	}
	
	void UpdateBeam(BeamEvaluator& evaluator, const double* ionG, const double* ionDL, const double* ionDM)
	{
		BeamEvaluator::PrecalcPosInfo posInfo;
		evaluator.PrecalculatePositionInfo(posInfo, _component.PosRA(), _component.PosDec());
		std::complex<double> beamGains[4];
		for(size_t ch=0; ch!=_bandData->ChannelCount(); ++ch)
		{
			size_t channelBlockIndex = ch * _channelBlocks / _bandData->ChannelCount();
			const double frequency = _bandData->ChannelFrequency(ch);
			evaluator.EvaluateAbsToApparentGain(posInfo, frequency, beamGains);
			_accumulators[ch]->UpdateBeam(beamGains, ionG[channelBlockIndex], ionDL[channelBlockIndex], ionDM[channelBlockIndex]);
		}
	}
	
	void Accumulate(std::complex<double>* data, const double weight, size_t ch, double u, double v, double w)
	{
		_accumulators[ch]->Add(data, weight, u, v, w);
	}
	
	void Finish()
	{
		for(std::vector<FluxAccumulator*>::iterator i=_accumulators.begin(); i!=_accumulators.end(); ++i)
			(*i)->Finish();
	}
	
	/**
	 * Only the accumulators are serialized, it is assumed the other values are
	 * set by initialization.
	 */
	virtual void Serialize(std::ostream& stream) const
	{
		SerializeToUInt64(stream, _accumulators.size());
		for(size_t i=0; i!=_accumulators.size(); ++i)
			_accumulators[i]->Serialize(stream);
	}
	
	virtual void Unserialize(std::istream& stream)
	{
		size_t size = UnserializeUInt64(stream);
		if(size != _accumulators.size())
			throw std::runtime_error("Unserializing flux spectrum accumulator with different number of channels: some of your input data do not match");
		for(size_t i=0; i!=size; ++i)
			_accumulators[i]->Unserialize(stream);
	}
	
	void AccumulateFromStream(std::istream& stream)
	{
		size_t size = UnserializeUInt64(stream);
		if(size != _accumulators.size())
			throw std::runtime_error("Unserializing flux spectrum accumulator with different number of channels: some of your input data do not match");
		for(size_t i=0; i!=size; ++i)
			_accumulators[i]->AccumulateFromStream(stream);
	}
	
	void SaveSpectrumAsTextFile(const char* filename) const
	{
		std::ofstream stream(filename);
		for(size_t ch=0; ch!=_bandData->ChannelCount(); ++ch)
		{
			stream << _bandData->ChannelFrequency(ch);
			double stokesValues[4];
			_accumulators[ch]->GetFlux(stokesValues);
			for(size_t p=0; p!=4; ++p)
				stream << '\t' << stokesValues[p];
		}
	}
	
	void GetSpectrum(ModelComponent& component) const
	{
		component = _component;
		SpectralEnergyDistribution sed;
		for(size_t ch=0; ch!=_bandData->ChannelCount(); ++ch)
		{
			Measurement m;
			m.SetFrequencyHz(_bandData->ChannelFrequency(ch));
			double stokesValues[4];
			_accumulators[ch]->GetFlux(stokesValues);
			for(size_t p=0; p!=4; ++p)
				m.SetFluxDensity(Polarization::IndexToStokes(p), stokesValues[p]);
			sed.AddMeasurement(m);
		}
		component.SetSED(sed);
	}
private:
	std::vector<FluxAccumulator*> _accumulators;
	ModelComponent _component;
	const BandData* _bandData;
	size_t _channelBlocks;
};

#endif
