#ifndef BANDDATA_H
#define BANDDATA_H

#include <stdexcept>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>

/**
 * Contains information about a single band ("spectral window").
 * A band consists of a sequence of contiguous channels.
 */
class BandData
{
	public:
		typedef std::reverse_iterator<double*> reverse_iterator;
		typedef std::reverse_iterator<const double*> const_reverse_iterator;
		
		/**
		 * Construct an empty instance.
		 */
		BandData() : _channelCount(0), _channelFrequencies(0), _frequencyStep(0.0)
		{
		}
		
		/**
		 * Construct an instance from a spectral window table. The spectral window table
		 * can only have a single entry, otherwise an exception is thrown.
		 * @param spwTable The CASA Measurement Set spectral window table.
		 */
		BandData(casacore::MSSpectralWindow& spwTable)
		{
			if(spwTable.nrow() != 1) throw std::runtime_error("Set should have exactly one spectral window");
			
			initFromTable(spwTable, 0);
		}
		
		/**
		 * Construct an instance from a specified entry of a spectral window table.
		 * @param spwTable The CASA Measurement Set spectral window table.
		 * @param bandIndex The entry index of the spectral window table.
		 */
		BandData(casacore::MSSpectralWindow& spwTable, size_t bandIndex)
		{
			initFromTable(spwTable, bandIndex);
		}
		
		/**
		 * Copy constructor.
		 * @param source Copied to the new banddata.
		 */
		BandData(const BandData& source) : _channelCount(source._channelCount), _frequencyStep(source._frequencyStep)
		{
			_channelFrequencies = new double[_channelCount];
			for(size_t index = 0; index != _channelCount; ++index)
			{
				_channelFrequencies[index] = source._channelFrequencies[index];
			}
		}
		
		/**
		 * Construct a new instance from a part of another band.
		 * @param source Instance that is partially copied.
		 * @param startChannel Start of range of channels that are copied.
		 * @param endChannel End of range, exclusive.
		 */
		BandData(const BandData &source, size_t startChannel, size_t endChannel) :
			_channelCount(endChannel - startChannel), _frequencyStep(source._frequencyStep)
		{
			if(_channelCount == 0) throw std::runtime_error("No channels in set");
			
			_channelFrequencies = new double[_channelCount];
			for(size_t index = 0; index != _channelCount; ++index)
			{
				_channelFrequencies[index] = source._channelFrequencies[index + startChannel];
			}
		}
		
		/**
		 * Construct a new BandData class and initialize it with an array of frequencies.
		 * @param channelCount Number of channels in the new instance.
		 * @param frequencies Array of @p channelCount doubles containing the channel frequencies.
		 */
		BandData(size_t channelCount, const double* frequencies) :
			_channelCount(channelCount)
		{
			_channelFrequencies = new double[channelCount];
			memcpy(_channelFrequencies, frequencies, sizeof(double)*channelCount);
			if(channelCount >= 2)
				_frequencyStep = _channelFrequencies[1] - _channelFrequencies[0];
			else
				_frequencyStep = 0.0;
		}
		
		/** Destructor. */
		~BandData()
		{
			delete[] _channelFrequencies;
		}
		
		/** Assignment operator */
		void operator=(const BandData& source)
		{
			_channelCount = source._channelCount;
			_frequencyStep = source._frequencyStep;
			if(_channelCount != 0)
			{
				_channelFrequencies = new double[_channelCount];
				for(size_t index = 0; index != _channelCount; ++index)
				{
					_channelFrequencies[index] = source._channelFrequencies[index];
				}
			}
			else {
				_channelFrequencies = 0;
			}
		}
		
		/** Iterator over frequencies, pointing to first channel */
		double* begin()
		{ return _channelFrequencies; }
		/** Iterator over frequencies, pointing past last channel */
		double* end()
		{ return _channelFrequencies+_channelCount; }
		/** Constant iterator over frequencies, pointing to first channel */
		const double* begin() const
		{ return _channelFrequencies; }
		/** Constant iterator over frequencies, pointing to last channel */
		const double* end() const
		{ return _channelFrequencies+_channelCount; }
		
		std::reverse_iterator<double*> rbegin()
		{ return std::reverse_iterator<double*>(end()); }
		
		std::reverse_iterator<double*> rend()
		{ return std::reverse_iterator<double*>(begin()); }
	
		std::reverse_iterator<const double*> rbegin() const
		{ return std::reverse_iterator<const double*>(end()); }
		
		std::reverse_iterator<const double*> rend() const
		{ return std::reverse_iterator<const double*>(begin()); }
	
		/**
		 * Assign new frequencies to this instance.
		 * @param channelCount Number of channels.
		 * @param frequencies Array of @p channelCount doubles containing the channel frequencies.
		 */
		void Set(size_t channelCount, const double* frequencies)
		{
			_channelCount = channelCount;
			delete[] _channelFrequencies;
			_channelFrequencies = new double[channelCount];
			memcpy(_channelFrequencies, frequencies, sizeof(double)*channelCount);
		}
		
		/** Retrieve number of channels in this band.
		 * @returns Number of channels.
		 */
		size_t ChannelCount() const { return _channelCount; }
		
		/** Get the frequency in Hz of a specified channel.
		 * @param channelIndex Zero-indexed channel index.
		 */
		double ChannelFrequency(size_t channelIndex) const
		{
			return _channelFrequencies[channelIndex];
		}
		
		/** Get the wavelength in m of a specified channel.
		 * @param channelIndex Zero-indexed channel index.
		 */
		double ChannelWavelength(size_t channelIndex) const
		{
			return 299792458.0L / _channelFrequencies[channelIndex];
		}
		
		/** Get the frequency of the last channel.
		 * @returns Highest frequency.
		 */
		double HighestFrequency() const
		{
			return _channelFrequencies[_channelCount-1];
		}
		
		/** Get the frequency of the first channel.
		 * @returns Lowest frequency.
		 */
		double LowestFrequency() const
		{
			return _channelFrequencies[0];
		}
		
		/** Get the centre frequency.
		 * @returns 0.5 * (HighestFrequency + LowestFrequency)
		 */
		double CentreFrequency() const
		{
			return (HighestFrequency() + LowestFrequency()) * 0.5;
		}
		
		/** Convert a frequency to a wavelength.
		 * @param frequencyHz Frequency in Hz.
		 * @returns Wavelength in m.
		 */
		static double FrequencyToLambda(double frequencyHz)
		{
			return 299792458.0L / frequencyHz;
		}
		
		/** Get the wavelength of the central channel.
		 * @returns Central channel wavelength.
		 */
		double CentreWavelength() const
		{
			return 299792458.0L / ((HighestFrequency() + LowestFrequency()) * 0.5);
		}
		
		/** Get the distance between channels in Hz.
		 * @returns Distance between channels.
		 */
		double FrequencyStep() const
		{
			return _frequencyStep;
		}
		
		/** Get the wavelength of the first channel.
		 * @returns longest wavelength. */
		double LongestWavelength() const
		{
			return ChannelWavelength(0);
		}
		
		/** Get the wavelength of the last channel.
		 * @returns smallest wavelength. */
		double SmallestWavelength() const
		{
			return ChannelWavelength(_channelCount-1);
		}
		
		/** Get the start of the frequency range covered by this band.
		 * @returns Start of the band in Hz. */
		double BandStart() const
		{
			return LowestFrequency() - FrequencyStep()*0.5;
		}
		/** Get the end of the frequency range covered by this band.
		 * @returns End of the band in Hz. */
		double BandEnd() const
		{
			return HighestFrequency() + FrequencyStep()*0.5;
		}
		
		/** Get the total bandwidth covered by this band.
		 * @returns Bandwidth in Hz. */
		double Bandwidth() const
		{
			return HighestFrequency() - LowestFrequency() + FrequencyStep();
		}
		
	private:
		void initFromTable(casacore::MSSpectralWindow& spwTable, size_t bandIndex)
		{
			casacore::ROScalarColumn<int> numChanCol(spwTable, casacore::MSSpectralWindow::columnName(casacore::MSSpectralWindowEnums::NUM_CHAN));
			int temp;
			numChanCol.get(bandIndex, temp);
			_channelCount = temp;
			if(_channelCount == 0) throw std::runtime_error("No channels in set");
			
			casacore::ROArrayColumn<double> chanFreqCol(spwTable, casacore::MSSpectralWindow::columnName(casacore::MSSpectralWindowEnums::CHAN_FREQ));
			casacore::Array<double> channelFrequencies;
			chanFreqCol.get(bandIndex, channelFrequencies, true);
			
			_channelFrequencies = new double[_channelCount];
			size_t index = 0;
			for(casacore::Array<double>::const_iterator i=channelFrequencies.begin();
					i != channelFrequencies.end(); ++i)
			{
				_channelFrequencies[index] = *i;
				++index;
			}
		  if(_channelCount > 1)
		    _frequencyStep = _channelFrequencies[1] - _channelFrequencies[0];
		  else
		    _frequencyStep = 0.0;
		}
		
		size_t _channelCount;
		double *_channelFrequencies;
		double _frequencyStep;
};

#endif
