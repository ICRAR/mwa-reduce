#ifndef ION_SPECTRUM_MAKER_H
#define ION_SPECTRUM_MAKER_H

#include <ms/MeasurementSets/MeasurementSet.h>
#include <boost/thread/thread.hpp>

#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include "banddata.h"
#include "fluxaccumulator.h"
#include "fluxspectrumaccumulator.h"
#include "lane.h"
#include "ionsolutionfile.h"
#include "uvector.h"

class IonSpectrumMaker
{
public:
	IonSpectrumMaker() { }
	
	void Accumulate(const char *msFilename, const char *ionFilename, const char *modelFilename)
	{
		casa::MeasurementSet ms(msFilename);
		std::string dataColumnName("DATA");
		IonSolutionFile ionSolutionFile;
		ionSolutionFile.OpenForReading(ionFilename);
		_model = Model(modelFilename);
		
		/**
			* Read some meta data from the measurement set
			*/
		_bandData = BandData(ms.spectralWindow());
		size_t channelCount = _bandData.ChannelCount();
		
		if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		casa::ROArrayColumn<casa::Complex> dataColumn(ms, dataColumnName);
		casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
		casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
		casa::MEpoch::ROScalarColumn timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
		casa::ROScalarColumn<double> timeAsDoubleColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
		casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
		casa::ROArrayColumn<double> uvwColumn(ms, ms.columnName(casa::MSMainEnums::UVW));

		casa::MSField fieldTable = ms.field();
		casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
		if(refDirColumn.nrow() != 1)
			throw std::runtime_error("Field table nrow != 1");
		casa::Array<double> refDir = refDirColumn(0);
		casa::Array<double>::const_iterator refDirIter = refDir.begin();
		long double phaseCentreRA = *refDirIter; ++refDirIter;
		long double phaseCentreDec = *refDirIter;
	
		casa::IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		if(polarizationCount != 4) throw std::runtime_error("Need 4 polarizations");
		
		double time = -1.0;
		size_t timestepCount = 0;
		for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
		{
			if(timeAsDoubleColumn(rowIndex) != time)
			{
				++timestepCount;
				time = timeAsDoubleColumn(rowIndex);
			}
		}
		
		_accumulatorPerSourceComponent.resize(_model.ComponentCount());
		size_t compIndex=0;
		for(Model::const_iterator s=_model.begin(); s!=_model.end(); ++s)
		{
			for(ModelSource::const_iterator c=s->begin(); c!=s->end(); ++c)
			{
				const ModelComponent &component = *c;
				_accumulatorPerSourceComponent[compIndex] = new FluxSpectrumAccumulator(component, &_bandData, ionSolutionFile.ChannelBlockCount(), phaseCentreRA, phaseCentreDec);
				++compIndex;
			}
		}
		
		_beamEvaluator = new BeamEvaluator(ms);
		
		ao::lane<RowData> lane(256);
		casa::Array<casa::Complex> dataArray(dataShape);
		casa::Array<float> weightArray(dataShape);
		
		_gSolutions.resize(ionSolutionFile.ChannelBlockCount()),
		_dlSolutions.resize(ionSolutionFile.ChannelBlockCount()),
		_dmSolutions.resize(ionSolutionFile.ChannelBlockCount());
		
		size_t timeIndex = 0, interval = 0;
		updateBeam(timeColumn(0), ionSolutionFile, interval);
		std::unique_ptr<boost::thread> processRowThread(new boost::thread(&IonSpectrumMaker::processRows, this, &lane));
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			size_t
				a1 = ant1Column(rowIndex),
				a2 = ant2Column(rowIndex);
			casa::MEpoch time = timeColumn(rowIndex);
			if(a1 != a2)
			{
				casa::Array<double> uvwArray = uvwColumn(rowIndex);
				casa::Array<double>::const_contiter uvwI = uvwArray.cbegin();
				double u = *uvwI; ++uvwI;
				double v = *uvwI; ++uvwI;
				double w = *uvwI;
				
				if(_beamEvaluator->Time().getValue() != time.getValue())
				{
					// Stop all threads, then update beam values, then restart threads.
					lane.write_end();
					processRowThread->join();
					
					lane.clear();
					updateBeam(time, ionSolutionFile, interval);
					processRowThread.reset(new boost::thread(&IonSpectrumMaker::processRows, this, &lane));
					
					++timeIndex;
					size_t nextIntervalStartTimestep = ((interval+1)*timestepCount) / ionSolutionFile.IntervalCount();
					if(timeIndex == nextIntervalStartTimestep)
					{
						++interval;
					}
					std::cout << "Starting timestep " << timeIndex << ',' << timestepCount << " of interval " << interval << '/' << ionSolutionFile.IntervalCount() << "...\n";
				}
				
				dataColumn.get(rowIndex, dataArray);
				weightColumn.get(rowIndex, weightArray);
				
				RowData rowData;
				rowData.uInM = u;
				rowData.vInM = v;
				rowData.wInM = w;
				rowData.data = new std::complex<float>[channelCount*4];
				rowData.weights = new float[channelCount*4];
				memcpy(rowData.data, dataArray.cbegin(), sizeof(std::complex<float>)*channelCount*4);
				memcpy(rowData.weights, weightArray.cbegin(), sizeof(float)*channelCount*4);
				lane.write(rowData);
			}
		}
		
		delete _beamEvaluator;
	}
private:
	struct RowData
	{
		RowData() : uInM(0.0), vInM(0.0), wInM(0.0), data(0), weights(0) { }
		double uInM, vInM, wInM;
		std::complex<float> *data;
		float *weights;
	};

	struct SampleData
	{
		size_t channelIndex, componentIndex;
		double u, v, w;
		std::complex<double> data[4];
		double weight;
	};

	void updateBeam(const casa::MEpoch& time, IonSolutionFile& ionSolutionFile, size_t interval)
	{
		std::cout << "Calculating beam gains...\n";
		_beamEvaluator->SetTime(time);
		size_t compIndex=0;
		for(size_t sIndex=0; sIndex!=_model.SourceCount(); ++sIndex)
		{
			for(size_t cb=0; cb!=ionSolutionFile.ChannelBlockCount(); ++cb)
			{
				IonSolutionFile::Solution solution;
				ionSolutionFile.ReadSolution(solution, interval, 0, sIndex, cb);
				_gSolutions[cb] = solution.gain;
				_dlSolutions[cb] = solution.dl;
				_dmSolutions[cb] = solution.dm;
			}
			for(ModelSource::const_iterator c=_model.Source(sIndex).begin(); c!=_model.Source(sIndex).end(); ++c)
			{
				_accumulatorPerSourceComponent[compIndex]->UpdateBeam(*_beamEvaluator, _gSolutions.data(), _dlSolutions.data(), _dmSolutions.data());
				++compIndex;
			}
		}
	}		
	
	void processRows(ao::lane<RowData>* lane)
	{
		size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
		std::vector<ao::lane<SampleData> > outLanes(cpuCount);
		boost::thread_group threadGroup;
		for(size_t c=0; c!=cpuCount; ++c)
		{
			outLanes[c].resize(256);
			threadGroup.add_thread(new boost::thread(&IonSpectrumMaker::processSamples, this, &outLanes[c]));
		}
		RowData rowData;
		ao::uvector<double> wavelengths(_bandData.ChannelCount());
		for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
			wavelengths[ch] = _bandData.ChannelWavelength(ch);
		
		while(lane->read(rowData))
		{
			size_t compIndex=0;
			for(std::vector<FluxSpectrumAccumulator*>::const_iterator c=_accumulatorPerSourceComponent.begin(); c!=_accumulatorPerSourceComponent.end(); ++c)
			{
				SampleData sample;
				sample.componentIndex = compIndex;
				for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
				{
					sample.weight = 0.0;
					for(size_t p=0; p!=4; ++p)
					{
						sample.data[p] = rowData.data[ch*4 + p];
						sample.weight += rowData.weights[ch*4 + p];
					}
					const double lambda = wavelengths[ch];
					sample.channelIndex = ch;
					sample.u = rowData.uInM/lambda;
					sample.v = rowData.vInM/lambda;
					sample.w = rowData.wInM/lambda;
					outLanes[ch%cpuCount].write(sample);
				}
				++compIndex;
			}
			delete[] rowData.data;
		}
		
		for(size_t c=0; c!=cpuCount; ++c)
			outLanes[c].write_end();
	}
	
	void processSamples(ao::lane<SampleData>* lane)
	{
		SampleData sample;
		while(lane->read(sample))
		{
			_accumulatorPerSourceComponent[sample.componentIndex]->Accumulate(sample.data, sample.weight, sample.channelIndex, sample.u, sample.v, sample.w);
		}
	}

	Model _model;
	BandData _bandData;
	BeamEvaluator *_beamEvaluator;
	std::vector<FluxSpectrumAccumulator*> _accumulatorPerSourceComponent;
	ao::uvector<double> _gSolutions, _dlSolutions, _dmSolutions;
};

#endif
