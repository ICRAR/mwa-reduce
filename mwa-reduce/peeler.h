#ifndef PEELER_H
#define PEELER_H

#include <ms/MeasurementSets/MeasurementSet.h>

#include <boost/thread/thread.hpp>

#include <queue>
#include "banddata.h"
#include "lane.h"
#include "model.h"

class Peeler
{
public:
	Peeler(casa::MeasurementSet& ms);
	
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
	
	void SetModel(const Model& model)
	{
		_model = model;
	}
	
	void SetRHSSolutionFile(const std::string& rhsSolutionFile)
	{
		_rhsSolutionFile = rhsSolutionFile;
	}
	
	void SetDataColumName(const std::string& dataColumnName)
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
	
	void SetSolutionInterval(size_t solutionInterval)
	{
		_solutionInterval = solutionInterval;
	}
	
private:
	struct ThreadData
	{
		ThreadData() { }
		
		boost::mutex *mutex;
		std::queue<size_t> *tasks;
		std::vector<class CalibrationMethod*> *calMethods;
		double limit;
		size_t nIter;
	};
	struct SubtractThreadInfo
	{
		size_t rowIndex;
		casa::Array<casa::Complex> *data;
		bool readyForWrite;
		double u, v, w;
		size_t a1, a2;
	};

	void calibrateThreadFunction(ThreadData data);
	
	casa::MeasurementSet& _ms;
	BandData _bandData;
	
	bool
		_beamOnSource,
		_applyBeam,
		_onlyScalar,
		_onlyDiag,
		_onlyRotation;
	std::string
		_modelFilename,
		_rhsSolutionFile,
		_dataColumnName;
	size_t _nIter;
	double _minAccuracy, _stoppingAccuracy, _minUVW;
	size_t _solutionInterval;
	Model _model;
};

#endif
