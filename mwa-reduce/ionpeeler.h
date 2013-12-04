#ifndef IONPEELER_H
#define IONPEELER_H

#include "banddata.h"
#include "model.h"
#include "visibilityarray.h"
#include "weightmode.h"

#include <mutex>
#include <thread>

class ImageWeights;
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
	
	void Peel(const char* msName, const char* modelName);
private:
	struct RowData
	{
		double u, v, w;
		size_t a1, a2, timeIndex;
	};

	void processChannel(size_t channelIndex);
	void processingThreadFunction(std::mutex* mutex, std::vector<size_t>* tasks);
	static bool isfinite(const std::complex<double>& val) { return std::isfinite(val.real()) && std::isfinite(val.imag()); }
	void initWeighting(casa::MeasurementSet& ms);
	
	size_t _solutionInterval;
	bool _applyBeam;
	std::string _dataColumnName;
	
	Model _model;
	std::vector<VisibilityArray<std::complex<double>, 4>*> _dataArrays;
	std::vector<VisibilityArray<double, 2>*> _weightArrays;
	std::vector<Predicter*> _predicters;
	std::vector<Model> _predictionModels;
	std::vector<RowData> _rowData;
	size_t _curStartRow, _curEndRow;
	BandData _bandData;
	size_t _cpuCount;
	
	WeightMode _weightMode;
	size_t _weightGridSize;
	double _weightPixelScale;
	std::unique_ptr<ImageWeights> _imageWeights;
	std::unique_ptr<ProgressBar> _progressBar;
};

#endif
