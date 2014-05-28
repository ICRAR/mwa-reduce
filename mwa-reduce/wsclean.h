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
#include "cachedimageset.h"

#include <set>

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
	void SetCleanBorderRatio(double borderRatio) { _cleanBorderRatio = borderRatio; }
	void SetThreshold(double threshold) { _threshold = threshold; }
	void SetColumnName(const std::string& columnName) { _columnName = columnName; }
	void SetPolarizations(const std::set<PolarizationEnum>& polarizations) { _polarizations = polarizations; }
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
	void SetSmallInversion(bool smallInversion) { _smallInversion = smallInversion; }
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
	void SetJoinPolarizations(bool joinPolarizations) { _joinedPolarizationCleaning = joinPolarizations; }
	void SetJoinChannels(bool joinChannels) { _joinedFrequencyCleaning = joinChannels; }
	void SetMultiscale(bool multiscale) { _multiscale = multiscale; }
	void SetMultiscaleThresholdBias(double thresholdBias)
	{ _multiscaleThresholdBias = thresholdBias; }
	void SetMultiscaleScaleBias(double scaleBias)
	{ _multiscaleScaleBias = scaleBias; }
	void SetMFSWeighting(bool mfsWeighting) { _mfsWeighting = mfsWeighting; }
	void SetWeightMode(enum WeightMode::WeightingEnum weighting) {
		_weightMode.SetMode(WeightMode(weighting));
	}
	void SetBriggsWeighting(double robustness) {
		_weightMode.SetMode(WeightMode::Briggs(robustness));
	}
	void SetSuperWeight(double superWeight) { _weightMode.SetSuperWeight(superWeight); }
	void SetBeamSize(double beamSize) { _manualBeamSize = beamSize; }
	void SetAntialiasingKernelSize(size_t kernelSize) { _antialiasingKernelSize = kernelSize; }
	void SetOversamplingFactor(size_t oversampling) { _overSamplingFactor = oversampling; }
	void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }
	void SetForceReorder(bool forceReorder) { _forceReorder = forceReorder; }
	void SetForceNoReorder(bool forceNoReorder) { _forceNoReorder = forceNoReorder; }
	void SetMemFraction(double memFraction) { _memFraction = memFraction; }
	void SetMemAbsLimit(double absMemLimit) { _absMemLimit = absMemLimit; }
	void SetMinUVW(double minUVW) { _globalSelection.SetMinUVW(minUVW); }
	void SetMaxUVW(double maxUVW) { _globalSelection.SetMaxUVW(maxUVW); }
	void SetWLimit(double wLimit) { _wLimit = wLimit; }
	void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }
	
	void AddInputMS(const std::string& msPath) { _filenames.push_back(msPath); }
	
	void RunClean();
	
	void RunPredict();
	
	bool JoinedFrequencyCleaning() const { return _joinedFrequencyCleaning; }
private:
	void runIndependentChannel(size_t outChannelIndex);
	void predictChannel(size_t outChannelIndex);
	
	void runFirstInversion(size_t outChannelIndex, PolarizationEnum polarization, size_t joinedChannelIndex);
	void performClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr);
	void performSimpleClean(class CleanAlgorithm& cleanAlgorithm, size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr, PolarizationEnum polarization);
	template<size_t PolCount>
	void performJoinedPolClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr);
	template<size_t PolCount>
	void performJoinedPolFreqClean(bool& reachedMajorThreshold, size_t majorIterationNr);
	void prepareInversionAlgorithm(PolarizationEnum polarization);
	
	void checkPolarizations();
	void performReordering();
	
	void initFitsWriter(class FitsWriter& writer);
	void setCleanParameters(class FitsWriter& writer, const class CleanAlgorithm& clean);
	void updateCleanParameters(class FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr);
	void initializeImageWeights(const MSSelection& partSelection);
	void initializeMFSImageWeights();
	void initializeCleanAlgorithm();
	void freeCleanAlgorithms();
	void initializeCurMSProviders(size_t currentChannelIndex, PolarizationEnum polarization);
	void clearCurMSProviders();
	void storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image);
	void selectChannels(MSSelection& selection, size_t outChannelIndex, size_t channelsOut);
	
	void imagePSF(size_t currentChannelIndex, size_t joinedChannelIndex);
	void imageGridding();
	void imageMainFirst(PolarizationEnum polarization, size_t joinedChannelIndex);
	void imageMainNonFirst(PolarizationEnum polarization, size_t joinedChannelIndex);
	void predict(PolarizationEnum polarization, size_t joinedChannelIndex);
	
	void makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary);
	void writeFits(const string& suffix, const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary);
	
	std::string getPSFPrefix(size_t channelIndex) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_channelsOut != 1)
		{
			partPrefixNameStr << '-';
			if(channelIndex < 1000) partPrefixNameStr << '0';
			if(channelIndex < 100) partPrefixNameStr << '0';
			if(channelIndex < 10) partPrefixNameStr << '0';
			partPrefixNameStr << channelIndex;
		}
		return partPrefixNameStr.str();
	}
	
	std::string getPrefix(PolarizationEnum polarization, size_t channelIndex, bool isImaginary) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_channelsOut != 1)
		{
			partPrefixNameStr << '-';
			if(channelIndex < 1000) partPrefixNameStr << '0';
			if(channelIndex < 100) partPrefixNameStr << '0';
			if(channelIndex < 10) partPrefixNameStr << '0';
			partPrefixNameStr << channelIndex;
		}
		if(_polarizations.size() != 1)
		{
			partPrefixNameStr << '-' << Polarization::TypeToShortString(polarization);
			if(isImaginary)
				partPrefixNameStr << 'i';
		}
		return partPrefixNameStr.str();
	}
	
	std::string getMFSPrefix(PolarizationEnum polarization, bool isImaginary) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_channelsOut != 1)
			partPrefixNameStr << "-MFS";
		if(_polarizations.size() != 1)
		{
			partPrefixNameStr << '-' << Polarization::TypeToShortString(polarization);
			if(isImaginary)
				partPrefixNameStr << 'i';
		}
		return partPrefixNameStr.str();
	}
	
	bool preferReordering() const
	{
		return ((_channelsOut != 1) || (_polarizations.size()>=4) || _forceReorder) && !_forceNoReorder;
	}
	
	size_t _imgWidth, _imgHeight, _channelsOut;
	double _pixelScaleX, _pixelScaleY, _threshold, _gain, _mGain, _cleanBorderRatio, _manualBeamSize, _memFraction, _absMemLimit, _wLimit, _multiscaleThresholdBias, _multiscaleScaleBias;
	size_t _nWLayers, _nIter, _antialiasingKernelSize, _overSamplingFactor, _threadCount;
	MSSelection _globalSelection, _currentPartSelection;
	std::string _columnName, _addModelFilename, _saveModelFilename, _cleanAreasFilename;
	std::set<PolarizationEnum> _polarizations;
	WeightMode _weightMode;
	std::string _prefixName;
	bool _allowNegative, _smallPSF, _smallInversion, _addApparentModel, _stopOnNegative, _makePSF;
	bool _forceReorder, _forceNoReorder, _joinedPolarizationCleaning, _joinedFrequencyCleaning, _mfsWeighting, _multiscale;
	enum LayeredImager::GridModeEnum _gridMode;
	std::vector<std::string> _filenames;
	std::string _commandLine;
	std::vector<size_t> _weightPerChannel;
	
	std::unique_ptr<class InversionAlgorithm> _inversionAlgorithm;
	std::unique_ptr<class ImageWeights> _imageWeights;
	std::vector<class CleanAlgorithm*> _cleanAlgorithms;
	std::unique_ptr<class AreaSet> _cleanAreas;
	ImageBufferAllocator<double> _imageAllocator;
	Stopwatch _inversionWatch, _predictingWatch, _cleaningWatch;
	bool _isFirstInversion, _doReorder;
	size_t _majorIterationNr;
	CachedImageSet _psfImages, _modelImages, _residualImages;
	std::vector<PartitionedMS::Handle> _partitionedMSHandles;
	FitsWriter _fitsWriter;
	std::vector<MSProvider*> _currentPolMSes;
	MultiBandData _firstMSBand;
};

#endif
