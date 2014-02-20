#ifndef WSCLEAN_H
#define WSCLEAN_H

#include "msprovider/msprovider.h"
#include "msprovider/partitionedms.h"

#include "layeredimager.h"
#include "msselection.h"
#include "polarizationenum.h"
#include "weightmode.h"
#include "imagebufferallocator.h"
#include "stopwatch.h"

class WSClean
{
public:
	WSClean();
	~WSClean();
	
	void SetImageSize(size_t width, size_t height) { _imgWidth = width; _imgHeight = height; }
	void SetPixelScale(double pixelScale) { _pixelScaleX = pixelScale; _pixelScaleY = pixelScale; }
	void SetNWlayers(size_t nWLayers) { _nWLayers = nWLayers; }
	void SetCleanGain(double gain) { _gain = gain; }
	void SetCleanMGain(double mGain) { _mGain = mGain; }
	void SetNIter(size_t nIter) { _nIter = nIter; }
	void SetThreshold(double threshold) { _threshold = threshold; }
	void SetColumnName(const std::string& columnName) { _columnName = columnName; }
	void SetPolarizations(const std::vector<PolarizationEnum>& polarizations) { _polarizations = polarizations; }
	void SetImaginaryPart(bool imagPart) { _imaginaryPart = imagPart; }
	void SetAllowNegative(bool allowNegative) { _allowNegative = allowNegative; }
	void SetStopOnNegative(bool stopOnNegative) { _stopOnNegative = stopOnNegative; }
	void SetMakePSF(bool makePSF) { _makePSF = makePSF; }
	void SetAddModelFilename(const std::string& modelFilename) { _addModelFilename = modelFilename; }
	void SetAddAppModel(bool addAppModel) { _addApparentModel = addAppModel; }
	void SetSaveModelFilename(const std::string& modelFilename) { _saveModelFilename = modelFilename; }
	void SetCleanAreasFilename(const std::string& filename) { _cleanAreasFilename = filename; }
	void SetPrefixName(const std::string& prefixName) { _prefixName = prefixName; }
	void SetGridMode(LayeredImager::GridModeEnum gridMode) { _gridMode = gridMode; }
	void SetSmallPSF(bool smallPSF) { _smallPSF = smallPSF; }
	void SetIntervalSelection(size_t startTimestep, size_t endTimestep) {
		_globalSelection.SetInterval(startTimestep, endTimestep);
	}
	void SetChannelSelection(size_t startChannel, size_t endChannel) {
		_globalSelection.SetChannelRange(startChannel, endChannel);
	}
	void SetFieldSelection(size_t fieldId) {
		_globalSelection.SetFieldId(fieldId);
	}
	void SetChannelsOut(size_t channelsOut) { _channelsOut = channelsOut; }
	void SetWeightMode(enum WeightMode::WeightingEnum weighting) {
		_weightMode.SetMode(WeightMode(weighting));
	}
	void SetBriggsWeighting(double robustness) {
		_weightMode.SetMode(WeightMode::Briggs(robustness));
	}
	void SetSuperWeight(double superWeight) { _weightMode.SetSuperWeight(superWeight); }
	void SetBeamSize(double beamSize) { _beamSize = beamSize; }
	void SetAntialiasingKernelSize(size_t kernelSize) { _antialiasingKernelSize = kernelSize; }
	void SetOversamplingFactor(size_t oversampling) { _overSamplingFactor = oversampling; }
	void SetForceReorder(bool forceReorder) { _forceReorder = forceReorder; }
	void SetForceNoReorder(bool forceNoReorder) { _forceNoReorder = forceNoReorder; }
	void SetMemFraction(double memFraction) { _memFraction = memFraction; }
	void SetWLimit(double wLimit) { _wLimit = wLimit; }
	void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }
	
	void AddInputMS(const std::string& msPath) { _filenames.push_back(msPath); }
	
	void Run();
private:
	void runChannel(size_t outChannelIndex);
	void runPolarization(size_t outChannelIndex, PolarizationEnum polarization);
	void prepareInversionAlgorithm(PolarizationEnum polarization);
	
	void initFitsWriter(class FitsWriter& writer);
	void setCleanParameters(class FitsWriter& writer, const class CleanAlgorithm& clean);
	void updateCleanParameters(class FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr);
	void initializeImageWeights(const MSSelection& partSelection);
	
	void imagePSF();
	void imageGridding();
	void imageMainFirst();
	void imageMainNonFirst();
	void predict(const double* modelImage);
	
	size_t _imgWidth, _imgHeight, _channelsOut;
	double _pixelScaleX, _pixelScaleY, _threshold, _gain, _mGain, _beamSize, _memFraction, _wLimit;
	size_t _nWLayers, _nIter, _antialiasingKernelSize, _overSamplingFactor;
	MSSelection _globalSelection, _currentPartSelection;
	std::string _columnName, _addModelFilename, _saveModelFilename, _cleanAreasFilename;
	std::vector<PolarizationEnum> _polarizations;
	WeightMode _weightMode;
	std::string _prefixName;
	bool _allowNegative, _smallPSF, _addApparentModel, _stopOnNegative, _imaginaryPart, _makePSF;
	bool _forceReorder, _forceNoReorder;
	enum LayeredImager::GridModeEnum _gridMode;
	std::vector<std::string> _filenames;
	std::string _commandLine;
	
	std::unique_ptr<class InversionAlgorithm> _inversionAlgorithm;
	std::unique_ptr<class ImageWeights> _imageWeights;
	ImageBufferAllocator<double> _imageAllocator;
	Stopwatch _inversionWatch, _predictingWatch, _cleaningWatch;
	bool _isFirstInversion, _doReorder;
	double *_psfImage;
	std::vector<PartitionedMS::Handle> _partitionedMSHandles;
};

#endif
