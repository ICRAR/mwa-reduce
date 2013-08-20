#ifndef SPECTRUM_MAKER_H
#define SPECTRUM_MAKER_H

#include "lane.h"
#include "model.h"

#include <complex>
#include <vector>
#include <memory>

class SpectrumMaker
{
private:
	struct ThreadTaskInfo
	{
		size_t sourceIndex;
		long double *flux;
		long double *fluxWeights;
		double u, v, w;
		class Predicter *predicter;
		std::complex<double>* data;
		bool* flags;
		size_t a1, a2;
	};

public:
	SpectrumMaker() { }
	
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
	
	void Measure()
	{
		// Set all sources to flux 1 Jy
		_spectrumPerSource.resize(_sources.size());
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			_sources[s].SetSED(SpectralEnergyDistribution(1.0, 1.0));
			_spectrumPerSource[s].Clear();
		}
		
		for(std::vector<std::pair<std::string,std::string>>::const_iterator i=_files.begin(); i!=_files.end(); ++i)
			measure(i->first, i->second);
		
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			_spectrumPerSource[s].Normalize();
		}
	}
	
	void FluxPerFrequency(std::map<double, long double>& dest, size_t sourceIndex, size_t polIndex) const
	{ _spectrumPerSource[sourceIndex].ToMap(dest, polIndex); }
private:
	void measure(const std::string& filename, const std::string& solutionsFile);
	
	void measureThreadFunc(const class BandData* bandData, lane<ThreadTaskInfo>* taskLane, size_t channelCount, size_t polarizationCount);
	
	struct Measurement
	{
		double flux[4];
		double weights[4];
		
		Measurement() {
			for(size_t p=0; p!=4; ++p)
			{
				flux[p] = 0.0;
				weights[p] = 0;
			}
		}
		
		void Normalize()
		{
			for(size_t p=0; p!=4; ++p)
			{
				if(weights[p] != 0)
					flux[p] /= weights[p];
				else
					flux[p] = std::numeric_limits<double>::quiet_NaN();
				
				weights[p] = 0;
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
			
			void AddMeasurement(double frequency, const long double *values, const long double *weights)
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
			
			void ToMap(std::map<double, long double>& dest, size_t polIndex) const
			{
				for(const_iterator i=begin(); i!=end(); ++i)
				{
					dest.insert(std::pair<double, long double>(i->first, i->second.flux[polIndex]));
				}
			}
		private:
			std::map<double, Measurement> _measurements;
	};
	
	std::vector<ModelSource> _sources;
	std::vector<std::pair<std::string,std::string>> _files;
	std::vector<Spectrum> _spectrumPerSource;
	Model _subtractedModel;
};

#endif
