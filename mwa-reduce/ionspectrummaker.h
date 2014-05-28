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
#include "buffered_lane.h"

#define ION_SPECTRUM_ROW_LANE_SIZE 1024
#define ION_SPECTRUM_ROW_BUFFER_SIZE 512
#define ION_SPECTRUM_SAMPLE_LANE_SIZE 65536
#define ION_SPECTRUM_SAMPLE_BUFFER_SIZE 32768

class IonSpectrumMaker
{
public:
	IonSpectrumMaker(const char *msFilename, const char *ionFilename, const char *modelFilename) :
		_ms(msFilename),
		_model(modelFilename),
		_bandData(_ms.spectralWindow())
	{
		_ionSolutionFile.OpenForReading(ionFilename);
		
		casa::MSField fieldTable = _ms.field();
		casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
		if(refDirColumn.nrow() != 1)
			throw std::runtime_error("Field table nrow != 1");
		casa::Array<double> refDir = refDirColumn(0);
		casa::Array<double>::const_iterator refDirIter = refDir.begin();
		long double phaseCentreRA = *refDirIter; ++refDirIter;
		long double phaseCentreDec = *refDirIter;
		
		_accumulatorPerSourceComponent.resize(_model.ComponentCount());
		size_t compIndex=0;
		for(Model::const_iterator s=_model.begin(); s!=_model.end(); ++s)
		{
			for(ModelSource::const_iterator c=s->begin(); c!=s->end(); ++c)
			{
				const ModelComponent &component = *c;
				_accumulatorPerSourceComponent[compIndex] = new FluxSpectrumAccumulator(component, &_bandData, _ionSolutionFile.ChannelBlockCount(), phaseCentreRA, phaseCentreDec);
				++compIndex;
			}
		}
	}
	
	~IonSpectrumMaker()
	{
		for(std::vector<FluxSpectrumAccumulator*>::iterator i=_accumulatorPerSourceComponent.begin(); i!=_accumulatorPerSourceComponent.end(); ++i)
			delete *i;
	}
	
	void Accumulate()
	{
		std::string dataColumnName("DATA");
		
		/**
			* Read some meta data from the measurement set
			*/
		size_t channelCount = _bandData.ChannelCount();
		
		if(_ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
		
		casa::ROArrayColumn<casa::Complex> dataColumn(_ms, dataColumnName);
		casa::ROArrayColumn<float> weightColumn(_ms, _ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
		casa::ROArrayColumn<bool> flagColumn(_ms, _ms.columnName(casa::MSMainEnums::FLAG));
		casa::MEpoch::ROScalarColumn timeColumn(_ms, _ms.columnName(casa::MSMainEnums::TIME));
		casa::ROScalarColumn<double> timeAsDoubleColumn(_ms, _ms.columnName(casa::MSMainEnums::TIME));
		casa::ROScalarColumn<int> ant1Column(_ms, _ms.columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> ant2Column(_ms, _ms.columnName(casa::MSMainEnums::ANTENNA2));
		casa::ROArrayColumn<double> uvwColumn(_ms, _ms.columnName(casa::MSMainEnums::UVW));
	
		casa::IPosition dataShape = dataColumn.shape(0);
		unsigned polarizationCount = dataShape[0];
		if(polarizationCount != 4) throw std::runtime_error("Need 4 polarizations");
		
		double time = -1.0;
		size_t timestepCount = 0;
		for(size_t rowIndex=0;rowIndex!=_ms.nrow();++rowIndex)
		{
			if(timeAsDoubleColumn(rowIndex) != time)
			{
				++timestepCount;
				time = timeAsDoubleColumn(rowIndex);
			}
		}
		
		_beamEvaluator.reset(new BeamEvaluator(_ms));
		
		ao::lane<RowData> internalLane(ION_SPECTRUM_ROW_LANE_SIZE);
		lane_write_buffer<RowData> bufferedLane(&internalLane, ION_SPECTRUM_ROW_BUFFER_SIZE);
		
		casa::Array<casa::Complex> dataArray(dataShape);
		casa::Array<float> weightArray(dataShape);
		casa::Array<bool> flagArray(dataShape);
		
		_gSolutions.resize(_ionSolutionFile.ChannelBlockCount()),
		_dlSolutions.resize(_ionSolutionFile.ChannelBlockCount()),
		_dmSolutions.resize(_ionSolutionFile.ChannelBlockCount());
		
		size_t timeIndex = 0, interval = 0;
		updateBeam(timeColumn(0), _ionSolutionFile, interval);
		std::cout << "Starting timestep " << timeIndex << '/' << timestepCount << " of interval " << interval << '/' << _ionSolutionFile.IntervalCount() << "...\n";
		std::unique_ptr<boost::thread> processRowThread(new boost::thread(&IonSpectrumMaker::processRows, this, &internalLane));
		for(size_t rowIndex=0; rowIndex!=_ms.nrow(); ++rowIndex)
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
					bufferedLane.write_end();
					processRowThread->join();
					
					bufferedLane.clear();
					updateBeam(time, _ionSolutionFile, interval);
					processRowThread.reset(new boost::thread(&IonSpectrumMaker::processRows, this, &internalLane));
					
					++timeIndex;
					size_t nextIntervalStartTimestep = ((interval+1)*timestepCount) / _ionSolutionFile.IntervalCount();
					if(timeIndex == nextIntervalStartTimestep)
					{
						++interval;
					}
					std::cout << "Starting timestep " << timeIndex << '/' << timestepCount << " of interval " << interval << '/' << _ionSolutionFile.IntervalCount() << "...\n";
				}
				
				dataColumn.get(rowIndex, dataArray);
				weightColumn.get(rowIndex, weightArray);
				flagColumn.get(rowIndex, flagArray);
				
				RowData rowData;
				rowData.uInM = u;
				rowData.vInM = v;
				rowData.wInM = w;
				rowData.data = new std::complex<float>[channelCount*4];
				rowData.weights = new float[channelCount*4];
				rowData.flags = new bool[channelCount*4];
				memcpy(rowData.data, dataArray.cbegin(), sizeof(std::complex<float>)*channelCount*4);
				memcpy(rowData.weights, weightArray.cbegin(), sizeof(float)*channelCount*4);
				memcpy(rowData.flags, flagArray.cbegin(), sizeof(bool)*channelCount*4);
				bufferedLane.write(rowData);
			}
		}
		
		bufferedLane.write_end();
		processRowThread->join();
		
		_beamEvaluator.reset();
	}
	
	void Save(const char* filename)
	{
		std::ofstream stream(filename);
		for(std::vector<FluxSpectrumAccumulator*>::iterator i=_accumulatorPerSourceComponent.begin(); i!=_accumulatorPerSourceComponent.end(); ++i)
		{
			(*i)->Finish();
			(*i)->Serialize(stream);
		}
	}
	
	void AccumulateFile(const char* filename)
	{
		std::ifstream stream(filename);
		for(std::vector<FluxSpectrumAccumulator*>::iterator i=_accumulatorPerSourceComponent.begin(); i!=_accumulatorPerSourceComponent.end(); ++i)
			(*i)->AccumulateFromStream(stream);
	}
	
	void GetModel(Model& model) const
	{
		size_t compIndex = 0;
		for(size_t s=0; s!=_model.SourceCount(); ++s)
		{
			const ModelSource& source = _model.Source(s);
			ModelSource newSource(source);
			newSource.ClearComponents();
			for(size_t c=0; c!=source.ComponentCount(); ++c)
			{
				ModelComponent newComponent;
				_accumulatorPerSourceComponent[compIndex]->GetSpectrum(newComponent);
				newSource.AddComponent(newComponent);
				
				++compIndex;
			}
			model.AddSource(newSource);
		}
	}
private:
	struct RowData
	{
		RowData() : uInM(0.0), vInM(0.0), wInM(0.0), data(0), weights(0) { }
		double uInM, vInM, wInM;
		std::complex<float> *data;
		float *weights;
		bool *flags;
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
				ionSolutionFile.ReadSolution(solution, interval, cb, 0, sIndex);
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
		lane_read_buffer<RowData> bufferedInputLane(lane, ION_SPECTRUM_ROW_BUFFER_SIZE);
		size_t cpuCount = (size_t) sysconf(_SC_NPROCESSORS_ONLN);
		std::vector<ao::lane<SampleData>*> outLanesInternal(cpuCount);
		std::vector<lane_write_buffer<SampleData>> bufferedOutLanes(cpuCount);
		boost::thread_group threadGroup;
		for(size_t c=0; c!=cpuCount; ++c)
		{
			outLanesInternal[c] = new ao::lane<SampleData>(ION_SPECTRUM_SAMPLE_LANE_SIZE);
			bufferedOutLanes[c].reset(outLanesInternal[c], ION_SPECTRUM_SAMPLE_BUFFER_SIZE);
			threadGroup.add_thread(new boost::thread(&IonSpectrumMaker::processSamples, this, outLanesInternal[c]));
		}
		RowData rowData;
		ao::uvector<double> wavelengths(_bandData.ChannelCount());
		for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
			wavelengths[ch] = _bandData.ChannelWavelength(ch);
		
		while(bufferedInputLane.read(rowData))
		{
			size_t compIndex=0;
			for(std::vector<FluxSpectrumAccumulator*>::const_iterator c=_accumulatorPerSourceComponent.begin(); c!=_accumulatorPerSourceComponent.end(); ++c)
			{
				SampleData sample;
				sample.componentIndex = compIndex;
				for(size_t ch=0; ch!=_bandData.ChannelCount(); ++ch)
				{
					sample.weight = 0.0;
					bool flagged = false;
					for(size_t p=0; p!=4; ++p)
					{
						sample.data[p] = rowData.data[ch*4 + p];
						sample.weight += rowData.weights[ch*4 + p];
						flagged = flagged || rowData.flags[ch*4 + p];
					}
					if(!flagged)
					{
						const double lambda = wavelengths[ch];
						sample.channelIndex = ch;
						sample.u = rowData.uInM/lambda;
						sample.v = rowData.vInM/lambda;
						sample.w = rowData.wInM/lambda;
						bufferedOutLanes[ch%cpuCount].write(sample);
					}
				}
				++compIndex;
			}
			delete[] rowData.data;
			delete[] rowData.weights;
			delete[] rowData.flags;
		}
		
		for(size_t c=0; c!=cpuCount; ++c)
			bufferedOutLanes[c].write_end();
		threadGroup.join_all();
		for(size_t c=0; c!=cpuCount; ++c)
			delete outLanesInternal[c];
	}
	
	void processSamples(ao::lane<SampleData>* lane)
	{
		lane_read_buffer<SampleData> bufferedLane(lane, ION_SPECTRUM_SAMPLE_BUFFER_SIZE);
		
		SampleData sample;
		while(bufferedLane.read(sample))
		{
			_accumulatorPerSourceComponent[sample.componentIndex]->Accumulate(sample.data, sample.weight, sample.channelIndex, sample.u, sample.v, sample.w);
		}
	}

	casa::MeasurementSet _ms;
	IonSolutionFile _ionSolutionFile;
	Model _model;
	BandData _bandData;
	std::unique_ptr<BeamEvaluator> _beamEvaluator;
	std::vector<FluxSpectrumAccumulator*> _accumulatorPerSourceComponent;
	ao::uvector<double> _gSolutions, _dlSolutions, _dmSolutions;
};

#endif
