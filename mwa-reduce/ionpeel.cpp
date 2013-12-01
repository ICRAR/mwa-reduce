#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <measures/Measures/MEpoch.h>
#include <measures/TableMeasures/ScalarMeasColumn.h>

#include "model.h"
#include "banddata.h"
#include "beamevaluator.h"
#include "predicter.h"
#include "visibilityarray.h"

/**
 * Approach:
 * - Read a few timesteps (the solution interval) from the measurement set
 * - For each source in the model (multi-thread this over channels)
 *   + Predict its visibility (using Predicter)
 *   + Find the global ionospheric phase term
 *   + Subtract the source with i-term from the visibilities
 * - Write the residual flux back to the measurement set
 */

class IonPeeler
{
public:
	void Peel(const char* msName, const char* modelName);
private:
	void processChannel(size_t channelIndex);
	
	Model _model;
	std::vector<std::unique_ptr<VisibilityArray<std::complex<double>, 2>>> _dataArrays;
	std::vector<std::unique_ptr<VisibilityArray<double, 2>>> _weightArrays;
	std::vector<std::unique_ptr<Predicter>> _predicters;
};

void IonPeeler::Peel(const char* msName, const char* modelName)
{
	casa::MeasurementSet ms(msName);
	_model = Model(modelName);
	
	bool applyBeam = true;
	size_t startRow = 0;
	std::string solutionFile;
	size_t solutionInterval = 1;
	
	if(ms.nrow() == 0) throw std::runtime_error("Table has no rows (no data)");
	
	casa::ROArrayColumn<casa::Complex> dataColumn(ms, ms.columnName(casa::MSMainEnums::DATA));
	casa::IPosition dataShape = dataColumn.shape(0);
	unsigned polarizationCount = dataShape[0];
	if(polarizationCount != 4)
		throw std::runtime_error("Expecting MS with 4 polarizations");
	
	BandData bandData(ms.spectralWindow());
	const size_t channelCount = bandData.ChannelCount();
	const size_t antennaCount = ms.antenna().nrow();
	
	casa::MSField fieldTable = ms.field();
	casa::ROArrayColumn<double> refDirColumn(fieldTable, fieldTable.columnName(casa::MSFieldEnums::REFERENCE_DIR));
	if(refDirColumn.nrow() != 1)
		throw std::runtime_error("Field table nrow != 1");
	casa::Array<double> refDir = refDirColumn(0);
	casa::Array<double>::const_iterator refDirIter = refDir.begin();
	long double phaseCentreRA = *refDirIter; ++refDirIter;
	long double phaseCentreDec = *refDirIter;
	// By setting the time beforehand, we don't waste time calculating a time step we don't need.
	casa::MEpoch::ROScalarColumn timeMColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	casa::MEpoch startTime = timeMColumn(startRow);
	BeamEvaluator beamEvaluator(ms);
	beamEvaluator.SetTime(startTime);
	
	std::cout << "Initializing model predicter for " << _model.SourceCount() << " sources...\n";
	std::vector<Model> predictionModels;
	for(size_t s=0; s!=_model.SourceCount(); ++s)
	{
		predictionModels.push_back(Model());
		predictionModels.back().AddSource(_model.Source(s));
		
		_predicters[s].reset(new Predicter(phaseCentreRA, phaseCentreDec, bandData.LowestFrequency(), bandData.HighestFrequency(), channelCount));
		if(applyBeam)
			_predicters[s]->Initialize(predictionModels.back(), solutionFile, &beamEvaluator);
		else
			_predicters[s]->Initialize(predictionModels.back(), solutionFile);
	}
	
	std::cout << "Counting timesteps...\n";
	double time = -1.0;
	ao::uvector<size_t> timestepRows;
	casa::ROScalarColumn<double> timeColumn(ms, ms.columnName(casa::MSMainEnums::TIME));
	for(size_t rowIndex=0;rowIndex!=ms.nrow();++rowIndex)
	{
		if(timeColumn(rowIndex) != time)
		{
			timestepRows.push_back(rowIndex);
			time = timeColumn(rowIndex);
		}
	}
	size_t timestepCount = timestepRows.size();
	timestepRows.push_back(ms.nrow());
	
	casa::ROArrayColumn<float> weightColumn(ms, ms.columnName(casa::MSMainEnums::WEIGHT_SPECTRUM));
	casa::ROArrayColumn<bool> flagColumn(ms, ms.columnName(casa::MSMainEnums::FLAG));
	casa::ROScalarColumn<int> ant1Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA1));
	casa::ROScalarColumn<int> ant2Column(ms, ms.columnName(casa::MSMainEnums::ANTENNA2));
	
	size_t passCount = (solutionInterval==0) ? 1 : (timestepCount + solutionInterval - 1) / solutionInterval;
	casa::Array<std::complex<float> > data(dataShape), modelData(dataShape);
	casa::Array<float> weights(dataShape);
	casa::Array<bool> flags(dataShape);
	_dataArrays.resize(channelCount);
	_weightArrays.resize(channelCount);
	for(size_t pass=0; pass!=passCount; ++pass)
	{
		const size_t
			startTimestep = timestepCount * pass / passCount,
			endTimestep = timestepCount * (pass+1) / passCount,
			timestepsInPass = endTimestep - startTimestep,
			startRow = timestepRows[startTimestep],
			endRow = timestepRows[endTimestep];
		
		std::cout << "Reading...\n";
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			_dataArrays[ch].reset(new VisibilityArray<std::complex<double>, 2>(1, antennaCount, timestepsInPass));
			_weightArrays[ch].reset(new VisibilityArray<double, 2>(1, antennaCount, timestepsInPass));
		}
		double time = timeColumn(startRow);
		size_t timeIndex = startTimestep;
		for(size_t rowIndex=startRow; rowIndex!=endRow; ++rowIndex)
		{
			size_t
				a1 = ant1Column(rowIndex),
				a2 = ant2Column(rowIndex);
			if(a1 != a2)
			{
				if(timeColumn(rowIndex) != time)
				{
					++timeIndex;
					time = timeColumn(rowIndex);
				}
				dataColumn.get(rowIndex, data);
				weightColumn.get(rowIndex, weights);
				flagColumn.get(rowIndex, flags);
				std::complex<float> *dataPtr = data.cbegin();
				float *weightsReadPtr = weights.cbegin();
				bool *flagPtr = flags.cbegin();
				
				for(size_t ch = 0; ch!=channelCount; ++ch)
				{
					size_t chIndex = ch * 4;
					for(size_t p=0; p!=4; ++p)
					{
						if(flagPtr[chIndex+p]) weightsReadPtr[chIndex+p] = 0.0;
					}
					
					std::complex<double> *arrPtr = _dataArrays[ch]->ValuePtr(a1, a2, timeIndex - startTimestep);
					arrPtr[0] = dataPtr[0];
					arrPtr[1] = dataPtr[3];
					double* weightsWritePtr = _weightArrays[ch]->ValuePtr(a1, a2, timeIndex - startTimestep);
					weightsWritePtr[0] = weightsReadPtr[0];
					weightsWritePtr[1] = weightsReadPtr[3];
				}
			}
		}
		
		std::cout << "Predicting...\n";
		// Data is in memory: predict visibilities
	}
}

void IonPeeler::processChannel(size_t channelIndex)
{
	for(size_t sourceIndex=0; sourceIndex!=_model.SourceCount(); ++sourceIndex)
	{
		// Predict visibility
		
		// Calculate phase term
		
		// Subtract source from visibilities in mem
	}
}

int main(int argc, char* argv[])
{
	if(argc != 3) {
		std::cout << "syntax: <ms> <model>";
		return -1;
	}
	IonPeeler peeler;
	peeler.Peel(argv[1], argv[2]);
}
