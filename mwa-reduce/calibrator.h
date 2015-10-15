#ifndef CALIBRATOR_H
#define CALIBRATOR_H

#include "model.h"
#include "solutionfile.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <boost/thread/mutex.hpp>

#include <queue>

class Calibrator
{
public:
	Calibrator(casacore::MeasurementSet& ms, size_t threadCount);
	
	void Perform();

	void SetBeamOnSource(bool beamOnSource)
	{
		_beamOnSource = beamOnSource;
	}
	
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
	
	void SetRHSSolutionFile(const std::string& rhsSolutionFile)
	{
		_rhsSolutionFilename = rhsSolutionFile;
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
	
	SolutionFile& GetSolutionFile() {
		return _solutionFile;
	}
private:
	struct ThreadData
	{
		ThreadData() { }
		
		boost::mutex *mutex;
		std::queue<size_t> *tasks;
		std::vector<class CalibrationMethod*> *calMethods;
	};
	
	void threadFunction(ThreadData data);

	casacore::MeasurementSet _ms;
	SolutionFile _solutionFile;
	std::string _modelFilename, _solutionFilename, _rhsSolutionFilename, _dataColumnName;
	Model _model;
	double _minAccuracy, _stoppingAccuracy;
	size_t _nIter, _solutionInterval, _threadCount;
	bool _onlyScalar, _onlyDiag, _onlyRotation;
	bool _beamOnSource, _applyBeam;
	double _minUVW, _maxUVW;
	bool _savePlotFiles, _saveFaradayPlotFiles, _saveCrossTermsPlotFile, _verbose;
	std::string _phasePlotFilename, _gainPlotFilename, _faradayPlotFilename, _crossTermsPlotFilename;
};

#endif
