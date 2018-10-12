#ifndef CALIBRATOR_H
#define CALIBRATOR_H

#include "solutionfile.h"
#include "msselection.h"

#include "model/model.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <memory>
#include <mutex>
#include <queue>
#include <vector>

class Calibrator
{
public:
	Calibrator(casacore::MeasurementSet& ms, size_t threadCount);
	
	void Perform();

	void SetApplyBeam(bool applyBeam)
	{
		_applyBeam = applyBeam;
	}
	
	void SetOnlyScalar(bool onlyScalar)
	{
		_onlyScalar = onlyScalar;
	}
	
	void SetOnlyDiag(bool onlyDiag)
	{
		_onlyDiag = onlyDiag;
	}
	
	void SetOnlyRotation(bool onlyRotation)
	{
		_onlyRotation = onlyRotation;
	}
	
	void SetModelFilename(const std::string& modelFilename)
	{
		_modelFilename = modelFilename;
	}
	
	void SetSolutionOutputFilename(const std::string& solutionOutputFilename)
	{
		_solutionFilename = solutionOutputFilename;
	}
	
	void SetModel(const Model& model)
	{
		_model = model;
	}
	
	void SetDataColumnName(const std::string& dataColumnName)
	{
		_dataColumnName = dataColumnName;
	}
	
	void SetNIter(size_t nIter)
	{
		_nIter = nIter;
	}
	
	void SetAccuracy(double minAccuracy, double stoppingAccuracy)
	{
		_minAccuracy = minAccuracy;
		_stoppingAccuracy = stoppingAccuracy;
	}
	
	void SetMinUVW(double minUVW)
	{
		_minUVW = minUVW;
	}
	
	void SetMaxUVW(double maxUVW)
	{
		_maxUVW = maxUVW;
	}
	
	void SetSolutionInterval(size_t solutionInterval)
	{
		_solutionInterval = solutionInterval;
	}
	
	void SetSolutionChannels(size_t nChannels)
	{
		_solutionChannels = nChannels;
	}
	
	void SetVerbose(bool verbose)
	{
		_verbose = verbose;
	}
	
	void SetSavePlotFiles(bool savePlotFiles)
	{
		_savePlotFiles = savePlotFiles;
	}
	
	void SetPlotFilenames(const std::string& phasePlotFilename, const std::string& gainPlotFilename)
	{
		_phasePlotFilename = phasePlotFilename;
		_gainPlotFilename = gainPlotFilename;
	}
	
	void SetAbsMem(double absmem)
	{
	  _absmem = absmem;
	}

	void SetInterval(size_t startTimestep, size_t endTimestep)
	{
		_selection.SetInterval(startTimestep, endTimestep);
	}

	void SetFollowAntenna(size_t followAntenna)
	{
		_followAntenna = followAntenna;
	}
	void SetMWAPath(const std::string& path)
	{
		_mwaPath = path;
	}	
	SolutionFile& GetSolutionFile() {
		return _solutionFile;
	}
private:
	struct ThreadData
	{
		ThreadData(size_t nChBlocks) : lastSuccesfulChBlock(nChBlocks) { }
		size_t lastSuccesfulChBlock;
	};
	std::vector<ThreadData> _threadData;
	std::vector<std::unique_ptr<class CalibrationMethod>> _calMethods;
	
	void calibrateChannelBlock(size_t channelBlockIndex, size_t threadIndex);

	casacore::MeasurementSet _ms;
	SolutionFile _solutionFile;
	std::string _modelFilename, _solutionFilename, _dataColumnName;
	Model _model;
	double _minAccuracy, _stoppingAccuracy;
	size_t _nIter, _solutionInterval, _solutionChannels, _threadCount;
	bool _onlyScalar, _onlyDiag, _onlyRotation;
	bool _applyBeam;
	double _minUVW, _maxUVW, _absmem;
	bool _savePlotFiles, _saveFaradayPlotFiles, _saveCrossTermsPlotFile, _verbose;
	std::string _phasePlotFilename, _gainPlotFilename, _faradayPlotFilename, _crossTermsPlotFilename, _mwaPath;
	MSSelection _selection;
	size_t _followAntenna;
};

#endif
