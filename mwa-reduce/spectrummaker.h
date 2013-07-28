#ifndef SPECTRUM_MAKER_H
#define SPECTRUM_MAKER_H

#include "modelsource.h"
#include "spectralenergydistribution.h"
#include "predicter.h"
#include "model.h"
#include "banddata.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <vector>
#include <memory>

class SpectrumMaker
{
public:
	SpectrumMaker();
	
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
		_measuredFlux.resize(_sources.size());
		_measureCount.resize(_sources.size());
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			_sources[s].SetSED(SpectralEnergyDistribution(1.0, 1.0));
			_measuredFlux[s].clear();
			_measureCount[s].clear();
		}
		
		for(std::vector<std::string>::const_iterator i=_filenames.begin(); i!=_filenames.end(); ++i)
			measure(*i);
		
		for(size_t s=0; s!=_sources.size(); ++s)
		{
			std::map<double, long unsigned>::const_iterator countPtr = _measureCount[s].begin();
			for(std::map<double, long double>::iterator fluxPtr = _measuredFlux[s].begin(); fluxPtr != _measuredFlux[s].end(); ++fluxPtr)
			{
				if(*countPtr != 0)
					fluxPtr->second /= countPtr->second;
				++countPtr;
			}
			_measureCount[s].clear();
		}
	}
	
	const std::map<double, long double>& FluxPerFrequency() const { return _measuredFlux; }
private:
	void measure(const std::string &filename)
	{
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
		
		typedef float num_t;
		typedef std::complex<num_t> complex_t;
		casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
		casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
		casa::ROArrayColumn<complex_t> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
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
		
		/**
		 * Calculate spectra
		 */
		casa::Array<complex_t> data(dataShape);
		casa::Array<bool> flags(dataShape);
		
		for(size_t rowIndex=0; rowIndex!=ms.nrow(); ++rowIndex)
		{
			// Cross correlation?
			size_t a1 = ant1Column.get(rowIndex), a2 = ant2Column.get(rowIndex);
			if(a1 != a2)
			{
				dataColumn.get(rowIndex, data);
				flagColumn.get(rowIndex, flags);
				casa::Array<double> uvwArray = uvwColumn(rowIndex);
				casa::Array<double>::const_iterator i = uvwArray.begin();
				double u = *i; ++i;
				double v = *i; ++i;
				double w = *i;
				
				std::vector<long double>::iterator measFluxIter = measFlux.begin();
				std::vector<unsigned long>::iterator measCountIter = measCount.begin();
				for(size_t s=0; s!=_sources.size(); ++s)
				{
					casa::Array<complex_t>::contiter dataPtr = data.cbegin();
					casa::Array<bool>::contiter flagPtr = flags.cbegin();
					for(size_t ch=0; ch!=channelCount; ++ch)
					{
						double lambda = bandData.ChannelWavelength(ch);
						for(size_t p=0; p!=polarizationCount; ++p)
						{
							float real = dataPtr->real(), imag = dataPtr->imag();
							if(!(*flagPtr) && std::isfinite(real) && std::isfinite(imag))
							{
								Predicter::CNumType predicted = predicters[s]->Predict(_sources[s], u/lambda, v/lambda, w/lambda, ch, p);
								// add real(data * conj(predicted))
								(*measFluxIter) += real * predicted.real() + imag * predicted.imag();
								(*measCountIter) ++;
								++measFluxIter;
								++measCountIter;
							}
						}
						++dataPtr;
						++flagPtr;
					}
				}
			}
		}
		
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			double freq = bandData.ChannelFrequency(ch);
			std::map<double, long double>::iterator i = _measuredFlux.find(freq);
			if(i == _measuredFlux.end())
			{
				_measuredFlux.insert(std::pair<double, long double>(freq, measFlux[ch]));
				_measureCount.insert(std::pair<double, long unsigned>(freq, measCount[ch]));
			}
			else {
				i->second += measFlux[ch];
				_measureCount.find(freq)->second += measCount[ch];
			}
		}
	}
	
	std::vector<ModelSource> _sources;
	std::vector<std::string> _filenames;
	std::vector<std::map<double, long double>> _measuredFlux;
	std::vector<std::map<double, long unsigned>> _measureCount;
};

#endif
