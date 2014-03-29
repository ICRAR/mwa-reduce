#ifndef IONPEELER_H
#define IONPEELER_H

#include "banddata.h"
#include "model.h"
#include "visibilityarray.h"
#include "weightmode.h"

#include <mutex>
#include <thread>

#ifdef HAVE_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#endif

class ImageWeights;
class IonSolutionFile;
class Predicter;
class ProgressBar;

class IonPeeler
{
private:
	typedef std::size_t size_t;
public:
	IonPeeler();
	~IonPeeler();
	
	void SetApplyBeam(bool applyBeam) { _applyBeam = applyBeam; }
	void SetSolutionInterval(size_t solutionInterval) { _solutionInterval = solutionInterval; }
	void SetDataColumnName(const std::string& dataColumnName) { _dataColumnName = dataColumnName; }
	void SetWeighting(WeightMode mode, size_t gridSize, double pixelScale)
	{
		_weightMode = mode;
		_weightGridSize = gridSize;
		_weightPixelScale = pixelScale;
	}
	void SetChannelBlockSize(size_t newSize) { _channelBlockSize = newSize; }
	void SetClusterFluxLimit(double limit) { _clusterFluxLimit = limit; }
	void SetDistanceLimit(double limit) { _distanceLimit = limit; }
	
	void Peel(const char* msName, const char* modelName, const char* solutionFilename);
	
	const Model& GetUsedModel() const { return _model; }
private:
	struct RowData
	{
		double u, v, w;
		size_t a1, a2, timeIndex;
	};
	struct FittingInfo
	{
		IonPeeler* ionPeeler;
		ao::uvector<std::complex<double>>* modelData;
		size_t channelBlockIndex, startChannel, endChannel;
		double lambda;
	};
	struct PeelingStats
	{
		size_t lsFits;
		size_t lsFittingIterations;
		double gSum, dlSum, dmSum;
		double gSumSq, dlSumSq, dmSumSq;
		PeelingStats() :
			lsFits(0),
			lsFittingIterations(0),
			gSum(0.0),
			dlSum(0.0),
			dmSum(0.0),
			gSumSq(0.0),
			dlSumSq(0.0),
			dmSumSq(0.0)
		{ }
		void operator+=(const PeelingStats& rhs)
		{
			lsFits += rhs.lsFits;
			lsFittingIterations += rhs.lsFittingIterations;
			gSum += rhs.gSum;
			dlSum += rhs.dlSum;
			dmSum += rhs.dmSum;
			gSumSq += rhs.gSumSq;
			dlSumSq += rhs.dlSumSq;
			dmSumSq += rhs.dmSumSq;
		}
	};

	void processChannel(size_t channelIndex, PeelingStats& stats);
	void processingThreadFunction(std::mutex* mutex, std::vector<size_t>* tasks);
	static bool isfinite(const std::complex<double>& val) { return std::isfinite(val.real()) && std::isfinite(val.imag()); }
	void initWeighting(casa::MeasurementSet& ms);
	
	void positionFitter(size_t channelIndex, PeelingStats& stats);
	
	void outputStats(const PeelingStats& stats);
	std::string radToString(double r);
	
	static int posMinimizationFunc(const gsl_vector *xvec, void *data, gsl_vector *f);
	static int posMinimizationFuncDeriv(const gsl_vector *xvec, void *data, gsl_matrix *J);
	static int posMinimizationFuncBoth(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J);
	
	size_t _solutionInterval;
	size_t _fitIterationCount;
	bool _applyBeam;
	std::string _dataColumnName;
	
	Model _model;
	std::vector<VisibilityArray<std::complex<double>, 4>*> _dataArrays;
	std::vector<VisibilityArray<double, 2>*> _weightArrays;
	std::vector<Predicter*> _predicters;
	std::vector<Model> _predictionModels;
	std::vector<RowData> _rowData;
	size_t _pass, _passCount;
	size_t _curStartRow, _curEndRow;
	size_t _startTimestep, _endTimestep;
	BandData _bandData;
	size_t _antennaCount, _channelBlockSize, _channelBlockCount;
	size_t _cpuCount;
	ao::uvector<size_t> _failedConvergencesPerSource, _failedConvergencesPerChannelGroup;
	struct PeelingStats _stats;
	
	WeightMode _weightMode;
	size_t _weightGridSize;
	double _weightPixelScale, _clusterFluxLimit, _distanceLimit;
	std::unique_ptr<ImageWeights> _imageWeights;
	std::unique_ptr<ProgressBar> _progressBar;
	std::unique_ptr<IonSolutionFile> _solutionFile;
};

#endif
