#ifndef SPECTRUM_MAKER_H
#define SPECTRUM_MAKER_H

#include "lane.h"
#include "model.h"
#include "banddata.h"
#include "matrix2x2.h"

#include <complex>
#include <vector>
#include <memory>

class BeamEvaluator;

class SpectrumMaker
{
private:
	struct ThreadTaskInfo
	{
		size_t sourceIndex;
		std::complex<double>* flux;
		std::complex<double>* fluxWeights;
		double u, v, w;
		class Predicter *predicter;
		std::complex<double>* beamWeights;
		std::complex<double>* data;
		float* dataWeights;
		bool* flags;
		size_t a1, a2;
	};
	
	struct BeamEvalTaskInfo
	{
		ModelSource *source;
		std::complex<double> *weights;
	};

public:
	SpectrumMaker();
	
	~SpectrumMaker();
	
	void AddSource(const class ModelSource &source)
	{
		_sources.push_back(source);
	}
	
	void AddMeasurementSet(const std::string &filename)
	{
		_files.push_back(std::make_pair(filename, std::string()));
	}
	
	void AddMeasurementSet(const std::string &filename, const std::string &solutionsFile)
	{
		_files.push_back(std::make_pair(filename, solutionsFile));
	}
	
	void SetSubtractedModel(const class Model& model)
	{
		_subtractedModel = model;
	}
	
	void SetApplyBeam(bool applyBeam)
	{
		_applyBeam = applyBeam;
	}
	
	void Measure()
	{
		// Set all sources to flux 1 Jy
		_spectrumPerSource.resize(_sources.size());
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			_sources[s].MakeUnitFlux();
			_spectrumPerSource[s].Clear();
		}
		
		for(std::vector<std::pair<std::string,std::string>>::const_iterator i=_files.begin(); i!=_files.end(); ++i)
			measure(i->first, i->second);
		
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			_spectrumPerSource[s].Normalize();
		}
	}
	
	void FluxPerFrequency(std::map<double, double>& dest, size_t sourceIndex, size_t polIndex) const
	{ _spectrumPerSource[sourceIndex].ToMap(dest, polIndex); }
private:
	void measure(const std::string& filename, const std::string& solutionsFile);
	
	void measureThreadFunc(lane<SpectrumMaker::ThreadTaskInfo>* taskLane);
	
	void recalculateBeamWeights(size_t beamWeightIndex);
	
	void recalculateBeamWeightsThreadFunc(lane<BeamEvalTaskInfo>* taskLane);
	
	struct Measurement
	{
		std::complex<double> flux[4];
		std::complex<double> weights[4];
		
		Measurement() {
			for(size_t p=0; p!=4; ++p)
			{
				flux[p] = 0.0;
				weights[p] = 0.0;
			}
		}
		
		void Normalize()
		{
			// Calculate: W*FW sum(W*W)^-1
			// Where sum(W*W) is in variable weights, and W*FW in flux.
			if(Matrix2x2::Invert(weights))
			{
				std::complex<double> temp[4];
				Matrix2x2::ATimesB(temp, weights, flux);
				//Matrix2x2::ATimesHermB(flux, temp, weights);
				Matrix2x2::Assign(flux, temp);
				weights[0] = 0.0; weights[1] = 0.0;
				weights[2] = 0.0; weights[3] = 0.0;
			}
			else {
				for(size_t p=0; p!=4; ++p)
					flux[p] = std::numeric_limits<double>::quiet_NaN();
			}
		}
	};
	
	class Spectrum
	{
		public:
			typedef std::map<double, Measurement>::iterator iterator;
			typedef std::map<double, Measurement>::const_iterator const_iterator;
			iterator begin() { return _measurements.begin(); }
			iterator end() { return _measurements.end(); }
			const_iterator begin() const { return _measurements.begin(); }
			const_iterator end() const { return _measurements.end(); }
			
			void Normalize() {
				for(iterator i=begin(); i!=end(); ++i)
					i->second.Normalize();
			}
			void Clear() { _measurements.clear(); }
			
			void AddMeasurement(double frequency, const std::complex<double>* values, const std::complex<double>* weights)
			{
				iterator i = _measurements.find(frequency);
				if(i == end())
				{
					Measurement m;
					for(size_t p=0; p!=4; ++p)
					{
						m.flux[p] = values[p];
						m.weights[p] = weights[p];
					}
					_measurements.insert(std::pair<double, Measurement>(frequency, m));
				}
				else {
					for(size_t p=0; p!=4; ++p)
					{
						i->second.flux[p] += values[p];
						i->second.weights[p] += weights[p];
					}
				}
			}
			
			void ToMap(std::map<double, double>& dest, size_t polIndex) const
			{
				for(const_iterator i=begin(); i!=end(); ++i)
				{
					dest.insert(std::make_pair(i->first, i->second.flux[polIndex].real()));
				}
			}
		private:
			std::map<double, Measurement> _measurements;
	};
	
	std::vector<std::complex<double>> _beamWeights[2];
	std::vector<ModelSource> _sources;
	std::vector<std::pair<std::string,std::string>> _files;
	std::vector<Spectrum> _spectrumPerSource;
	Model _subtractedModel;
	std::unique_ptr<BeamEvaluator> _beamEvaluator;
	BandData _bandData;
	bool _applyBeam;
};

#endif
