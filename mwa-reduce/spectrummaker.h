#ifndef SPECTRUM_MAKER_H
#define SPECTRUM_MAKER_H

#include "modelsource.h"
#include "spectralenergydistribution.h"
#include "predicter.h"
#include "model.h"
#include "banddata.h"
#include "beamevaluator.h"
#include "progressbar.h"
#include "lane.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <boost/thread/thread.hpp>

#include <vector>
#include <memory>

class SpectrumMaker
{
private:
	struct ThreadTaskInfo
	{
		size_t sourceIndex;
		long double *flux;
		unsigned long *count;
		double u, v, w;
		Predicter *predicter;
		casa::Complex* data;
		bool* flags;
	};

public:
	SpectrumMaker() { }
	
	void AddSource(const ModelSource &source)
	{
		_sources.push_back(source);
	}
	
	void AddMeasurementSet(const std::string &filename)
	{
		_filenames.push_back(filename);
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
		
		for(std::vector<std::string>::const_iterator i=_filenames.begin(); i!=_filenames.end(); ++i)
			measure(*i);
		
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			_spectrumPerSource[s].Normalize();
		}
	}
	
	void FluxPerFrequency(std::map<double, long double>& dest, size_t sourceIndex, size_t polIndex) const
	{ _spectrumPerSource[sourceIndex].ToMap(dest, polIndex); }
private:
	void measure(const std::string &filename)
	{
		ProgressBar progress(std::string("Measure spectra in ") + filename);
		casa::MeasurementSet ms(filename);
		
		/**
		 * Read some meta data from the measurement set
		 */
		BandData bandData(ms.spectralWindow());
		size_t channelCount = bandData.ChannelCount();
		
		casa::MSField fieldTable = ms.field();
		casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
		if(refDirColumn.nrow() != 1)
			throw std::runtime_error("Field table nrow != 1");
		casa::Array<double> refDir = refDirColumn(0);
		casa::Array<double>::const_iterator refDirIter = refDir.begin();
		long double phaseCentreRA = *refDirIter; ++refDirIter;
		long double phaseCentreDec = *refDirIter;
		
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
		casa::ROArrayColumn<casa::Complex> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
		casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
		casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));
		
		casa::IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		
		std::vector<long double> measFlux(channelCount * _sources.size() * 4);
		std::vector<unsigned long> measCount(channelCount * _sources.size() * 4);
		
		std::vector<std::unique_ptr<Predicter>> predicters;
		for(std::vector<ModelSource>::iterator sourceIter=_sources.begin();
				sourceIter!=_sources.end(); ++sourceIter)
		{
			predicters.push_back(std::unique_ptr<Predicter>(
				new Predicter(phaseCentreRA, phaseCentreDec, bandData.LowestFrequency(), bandData.HighestFrequency(), channelCount)));
			(*predicters.rbegin())->Initialize(*sourceIter);
		}
		
		const size_t BUFFER_COUNT = 16;
		
		size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
		boost::thread_group threadGroup;
		std::vector<lane<ThreadTaskInfo>*> taskLanes(cpuCount);
		for(size_t i=0; i!=cpuCount; ++i)
		{
			// The task lane must issue an "wait" when its size is
			// BUFFER_COUNT-2, because the pushing thread can start
			// overwritting the n-1 th buffer while the popping thread
			// can just get the n-1 th buffer out.
			taskLanes[i] = new lane<ThreadTaskInfo>(BUFFER_COUNT-2);
			threadGroup.add_thread(new boost::thread(&SpectrumMaker::measureThreadFunc, this, &bandData, &*taskLanes[i], channelCount, polarizationCount));
		}
		
		/**
		 * Calculate spectra
		 */
		std::vector<casa::Array<casa::Complex>*> dataBuffers(BUFFER_COUNT);
		std::vector<casa::Array<bool>*> flagBuffers(BUFFER_COUNT);
		for(size_t i=0; i!=BUFFER_COUNT; ++i)
		{
			dataBuffers[i] = new casa::Array<casa::Complex>(dataShape);
			flagBuffers[i] = new casa::Array<bool>(dataShape);
		}
		
		size_t bufferIndex = 0;
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			progress.SetProgress(rowIndex, ms.nrow());
			// Cross correlation?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 != a2)
			{
				casa::Array<casa::Complex> &data = *dataBuffers[bufferIndex];
				casa::Array<bool> &flags = *flagBuffers[bufferIndex];
				
				dataColumn.get(rowIndex, data);
				flagColumn.get(rowIndex, flags);
				casa::Array<double> uvwArray = uvwColumn(rowIndex);
				casa::Array<double>::const_iterator i = uvwArray.begin();
				double u = *i; ++i;
				double v = *i; ++i;
				double w = *i;
				
				for(size_t s=0; s!=_sources.size(); ++s)
				{
					size_t thread = s % cpuCount;
					ThreadTaskInfo task;
					task.flux = &measFlux[s * channelCount * polarizationCount];
					task.count = &measCount[s * channelCount * polarizationCount];
					task.data = data.cbegin();
					task.flags = flags.cbegin();
					task.u = u;
					task.v = v;
					task.w = w;
					task.sourceIndex = s;
					task.predicter = &*predicters[s];
					taskLanes[thread]->write(task);
				}
				
				bufferIndex = (bufferIndex + 1) % BUFFER_COUNT;
			}
		}
		for(size_t i=0; i!=cpuCount; ++i)
			taskLanes[i]->write_end();
		threadGroup.join_all();
		
		for(size_t i=0; i!=cpuCount; ++i)
			delete taskLanes[i];
		
		for(size_t i=0; i!=BUFFER_COUNT; ++i)
		{
			delete dataBuffers[i];
			delete flagBuffers[i];
		}
		
		progress.SetProgress(ms.nrow(), ms.nrow());
		
		// Apply beam
		BeamEvaluator beamEval(ms);
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				const ModelSource& src = _sources[s];
				beamEval.ApparentToAbs(src.PosRA(), src.PosDec(), bandData.ChannelFrequency(ch), &measCount[(s * channelCount + ch) * 4]);
			}
		}
		
		// Add to the total values
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			Spectrum &spectrum = _spectrumPerSource[s];
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				double freq = bandData.ChannelFrequency(ch);
				spectrum.AddMeasurement(freq, &measFlux[(s * channelCount + ch) * 4], &measCount[(s * channelCount + ch) * 4]);
			}
		}
	}
	
	void measureThreadFunc(const BandData *bandData, lane<ThreadTaskInfo> *taskLane, size_t channelCount, size_t polarizationCount)
	{
		ThreadTaskInfo taskInfo;
		while(taskLane->read(taskInfo))
		{
			casa::Complex* dataPtr = taskInfo.data;
			bool* flagPtr = taskInfo.flags;
			long double* measFluxIter = taskInfo.flux;
			long unsigned* measCountIter = taskInfo.count;
			Predicter& predicter = *taskInfo.predicter;
			double
				u = taskInfo.u,
				v = taskInfo.v,
				w = taskInfo.w;
			size_t s = taskInfo.sourceIndex;
			
			for(size_t ch=0; ch!=channelCount; ++ch)
			{
				double lambda = bandData->ChannelWavelength(ch);
				Predicter::CNumType predicted = predicter.Predict(_sources[s], u/lambda, v/lambda, w/lambda, ch, 0);
				for(size_t p=0; p!=polarizationCount; ++p)
				{
					float real = dataPtr->real(), imag = dataPtr->imag();
					if(!(*flagPtr) && std::isfinite(real) && std::isfinite(imag))
					{
						// add real(data * conj(predicted))
						(*measFluxIter) += real * predicted.real() + imag * predicted.imag();
						(*measCountIter) ++;
					}
					++measFluxIter;
					++measCountIter;
					++dataPtr;
					++flagPtr;
				}
			}
		}
	}
	
	struct Measurement
	{
		long double flux[4];
		long unsigned count[4];
		
		Measurement() {
			for(size_t p=0; p!=4; ++p)
			{
				flux[p] = 0.0;
				count[p] = 0;
			}
		}
		
		void Normalize()
		{
			for(size_t p=0; p!=4; ++p)
			{
				if(count[p] != 0)
				{
					flux[p] /= count[p];
					count[p] = 0;
				}
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
			
			void AddMeasurement(double frequency, const long double *values, const long unsigned *counts)
			{
				iterator i = _measurements.find(frequency);
				if(i == end())
				{
					Measurement m;
					for(size_t p=0; p!=4; ++p)
					{
						m.flux[p] = values[p];
						m.count[p] = counts[p];
					}
					_measurements.insert(std::pair<double, Measurement>(frequency, m));
				}
				else {
					for(size_t p=0; p!=4; ++p)
					{
						i->second.flux[p] += values[p];
						i->second.count[p] += counts[p];
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
	std::vector<std::string> _filenames;
	std::vector<Spectrum> _spectrumPerSource;
};

#endif
