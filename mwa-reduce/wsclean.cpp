#include "wsclean.h"

#include "areaset.h"
#include "beamevaluator.h"
#include "fitswriter.h"
#include "imageweights.h"
#include "inversionalgorithm.h"
#include "modelrenderer.h"
#include "model.h"
#include "msselection.h"
#include "wsinversion.h"

#include "cleanalgorithms/cleanalgorithm.h"
#include "cleanalgorithms/joinedpolclean.h"
#include "cleanalgorithms/simpleclean.h"

#include "msprovider/contiguousms.h"

#include "parser/areaparser.h"

#include <iostream>
#include <memory>

std::string commandLine;

WSClean::WSClean() :
	_imgWidth(2048), _imgHeight(2048), _channelsOut(1),
	_pixelScaleX(0.01 * M_PI / 180.0), _pixelScaleY(0.01 * M_PI / 180.0),
	_threshold(0.0), _gain(0.1), _mGain(1.0), _manualBeamSize(0.0), _memFraction(1.0), _wLimit(0.0),
	_nWLayers(0), _nIter(0), _antialiasingKernelSize(7), _overSamplingFactor(63),
	_globalSelection(),
	_columnName(), _addModelFilename(), _saveModelFilename(), _cleanAreasFilename(),
	_polarizations(),
	_weightMode(WeightMode::UniformWeighted),
	_prefixName("wsclean"),
	_allowNegative(true), _smallPSF(false), _smallInversion(false), _stopOnNegative(false), _makePSF(false),
	_forceReorder(false), _forceNoReorder(false), _joinedPolarizationCleaning(false),
	_gridMode(LayeredImager::KaiserBessel),
	_filenames(),
	_commandLine(),
	_inversionWatch(false), _predictingWatch(false), _cleaningWatch(false),
	_isFirstInversion(true)
{
	_polarizations.insert(Polarization::StokesI);
}

WSClean::~WSClean()
{ }

void WSClean::initFitsWriter(FitsWriter& writer)
{
	double
		ra = _inversionAlgorithm->PhaseCentreRA(),
		dec = _inversionAlgorithm->PhaseCentreDec(),
		pixelScaleX = _inversionAlgorithm->PixelSizeX(),
		pixelScaleY = _inversionAlgorithm->PixelSizeY(),
		freqHigh = _inversionAlgorithm->HighestFrequencyChannel(),
		freqLow = _inversionAlgorithm->LowestFrequencyChannel(),
		freqCentre = (freqHigh + freqLow) * 0.5,
		bandwidth = _inversionAlgorithm->BandEnd() - _inversionAlgorithm->BandStart(),
		beamSize = _inversionAlgorithm->BeamSize(),
		dateObs = _inversionAlgorithm->StartTime();
		
	writer.SetImageDimensions(_inversionAlgorithm->ImageWidth(), _inversionAlgorithm->ImageHeight(), ra, dec, pixelScaleX, pixelScaleY);
	writer.SetFrequency(freqCentre, bandwidth);
	writer.SetDate(dateObs);
	writer.SetPolarization(_inversionAlgorithm->Polarization());
	writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
	writer.AddHistory(commandLine);
	if(_manualBeamSize != 0.0) {
		double beamSizeRad = _manualBeamSize * (M_PI / 60.0 / 180.0);
		writer.SetBeamInfo(beamSizeRad, beamSizeRad, 0.0);
	}
	else {
		writer.SetBeamInfo(beamSize);
	}
	if(_inversionAlgorithm->HasDenormalPhaseCentre())
		writer.SetPhaseCentreShift(_inversionAlgorithm->PhaseCentreDL(), _inversionAlgorithm->PhaseCentreDM());
	
	writer.SetExtraKeyword("WSCNWLAY", _inversionAlgorithm->WGridSize());
	writer.SetExtraKeyword("WSCDATAC", _inversionAlgorithm->DataColumnName());
	writer.SetExtraKeyword("WSCWEIGH", _inversionAlgorithm->Weighting().ToString());
	writer.SetExtraKeyword("WSCGKRNL", _inversionAlgorithm->AntialiasingKernelSize());
	if(_inversionAlgorithm->Selection().HasChannelRange())
	{
		writer.SetExtraKeyword("WSCCHANS", _inversionAlgorithm->Selection().ChannelRangeStart());
		writer.SetExtraKeyword("WSCCHANE", _inversionAlgorithm->Selection().ChannelRangeEnd());
	}
	if(_inversionAlgorithm->Selection().HasInterval())
	{
		writer.SetExtraKeyword("WSCTIMES", _inversionAlgorithm->Selection().IntervalStart());
		writer.SetExtraKeyword("WSCTIMEE", _inversionAlgorithm->Selection().IntervalEnd());
	}
	writer.SetExtraKeyword("WSCFIELD", _inversionAlgorithm->Selection().FieldId());
}

void WSClean::setCleanParameters(FitsWriter& writer, const CleanAlgorithm& clean)
{
	writer.SetExtraKeyword("WSCNITER", clean.MaxNIter());
	writer.SetExtraKeyword("WSCTHRES", clean.Threshold());
	writer.SetExtraKeyword("WSCGAIN", clean.SubtractionGain());
	writer.SetExtraKeyword("WSCMGAIN", clean.StopGain());
	writer.SetExtraKeyword("WSCNEGCM", clean.AllowNegativeComponents());
	writer.SetExtraKeyword("WSCNEGST", clean.StopOnNegativeComponents());
	writer.SetExtraKeyword("WSCSMPSF", clean.ResizePSF());
}

void WSClean::updateCleanParameters(FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr)
{
	writer.SetExtraKeyword("WSCMINOR", minorIterationNr);
	writer.SetExtraKeyword("WSCMAJOR", majorIterationNr);
}

void WSClean::imagePSF()
{
	std::cout << std::flush << " == Constructing PSF ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoImagePSF(true);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
		
	CleanAlgorithm::RemoveNaNsInPSF(_inversionAlgorithm->ImageRealResult(), _imgWidth, _imgHeight);
	initFitsWriter(_fitsWriter);
	_psfImages.Initialize(_fitsWriter, _polarizations.size(), _prefixName + "-psf", _imageAllocator);
	_psfImages.Store(_inversionAlgorithm->ImageRealResult(), *_polarizations.begin(), false);
	_inversionWatch.Pause();
	
	_isFirstInversion = false;
	std::cout << "Beam size is " << _inversionAlgorithm->BeamSize()*(180.0*60.0/M_PI) << " arcmin.\n";
	
	std::cout << "Writing psf image... " << std::flush;
	_fitsWriter.Write(_prefixName + "-psf.fits", _inversionAlgorithm->ImageRealResult());
	std::cout << "DONE\n";
}

void WSClean::imageGridding()
{
	std::cout << "Writing gridding correction image... " << std::flush;
	double* gridding = _imageAllocator.Allocate(_imgWidth * _imgHeight);
	_inversionAlgorithm->GetGriddingCorrectionImage(&gridding[0]);
	FitsWriter fitsWriter;
	initFitsWriter(fitsWriter);
	fitsWriter.Write(_prefixName + "-gridding.fits", &gridding[0]);
	_imageAllocator.Free(gridding);
	std::cout << "DONE\n";
}

void WSClean::imageMainFirst(PolarizationEnum polarization)
{
	std::cout << std::flush << " == Constructing image ==\n";
	_inversionWatch.Start();
	if(_nWLayers != 0)
		_inversionAlgorithm->SetWGridSize(_nWLayers);
	else
		_inversionAlgorithm->SetNoWGridSize();
	_inversionAlgorithm->SetDoImagePSF(false);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	_inversionAlgorithm->SetVerbose(false);
	
	storeAndCombineXYandYX(_residualImages, polarization, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		storeAndCombineXYandYX(_residualImages, polarization, true, _inversionAlgorithm->ImageImaginaryResult());
}

void WSClean::imageMainNonFirst(PolarizationEnum polarization)
{
	std::cout << std::flush << " == Constructing image ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoSubtractModel(true);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	
	storeAndCombineXYandYX(_residualImages, polarization, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		storeAndCombineXYandYX(_residualImages, polarization, true, _inversionAlgorithm->ImageImaginaryResult());
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, bool isImaginary, const double* image)
{
	if(polarization == Polarization::YX && _polarizations.count(Polarization::XY)!=0)
	{
		std::cout << "Adding XY and YX together...\n";
		double
			*xyImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
		dest.Load(xyImage, Polarization::XY, isImaginary);
		size_t count = _imgWidth*_imgHeight;
		if(isImaginary)
		{
			for(size_t i=0; i!=count; ++i)
				xyImage[i] = (xyImage[i]-image[i])*0.5;
		}
		else {
			for(size_t i=0; i!=count; ++i)
				xyImage[i] = (xyImage[i]+image[i])*0.5;
		}
		dest.Store(xyImage, Polarization::XY, isImaginary);
		_imageAllocator.Free(xyImage);
	}
	else {
		dest.Store(image, polarization, isImaginary);
	}
}

void WSClean::predict(PolarizationEnum polarization)
{
	std::cout << std::flush << " == Converting model image to visibilities ==\n";
	const size_t size = _imgWidth*_imgHeight;
	double
		*modelImageReal = _imageAllocator.Allocate(size),
		*modelImageImaginary = 0;
		
	if(polarization == Polarization::YX)
	{
		_modelImages.Load(modelImageReal, Polarization::XY, false);
		modelImageImaginary = _imageAllocator.Allocate(size);
		_modelImages.Load(modelImageImaginary, Polarization::XY, true);
		for(size_t i=0; i!=size; ++i)
			modelImageImaginary[i] = -modelImageImaginary[i];
	}
	else {
		_modelImages.Load(modelImageReal, polarization, false);
		if(Polarization::IsComplex(polarization))
		{
			modelImageImaginary = _imageAllocator.Allocate(size);
			_modelImages.Load(modelImageImaginary, polarization, true);
		}
	}
	
	_predictingWatch.Start();
	_inversionAlgorithm->SetAddToModel(false);
	if(Polarization::IsComplex(polarization))
		_inversionAlgorithm->Predict(modelImageReal, modelImageImaginary);
	else
		_inversionAlgorithm->Predict(modelImageReal);
	_predictingWatch.Pause();
	_imageAllocator.Free(modelImageReal);
	_imageAllocator.Free(modelImageImaginary);
}

void WSClean::initializeImageWeights(const MSSelection& partSelection)
{
	if(_weightMode.RequiresGridding())
	{
		std::cout << "Precalculating weights for " << _weightMode.ToString() << " weighting... " << std::flush;
		_imageWeights.reset(new ImageWeights(_imgWidth, _imgHeight, _pixelScaleX, _pixelScaleY, _weightMode.SuperWeight()));
		for(size_t i=0; i!=_inversionAlgorithm->MeasurementSetCount(); ++i)
		{
			_imageWeights->Grid(_inversionAlgorithm->MeasurementSet(i), _weightMode, partSelection);
			if(_inversionAlgorithm->MeasurementSetCount() > 1)
				std::cout << i << ' ' << std::flush;
		}
		_inversionAlgorithm->SetPrecalculatedWeightInfo(_imageWeights.get());
		std::cout << "DONE\n";
	}
}

void WSClean::initializeCleanAlgorithm()
{
	if(_joinedPolarizationCleaning)
		_cleanAlgorithm.reset(new JoinedPolClean());
	else
		_cleanAlgorithm.reset(new SimpleClean());
	_cleanAlgorithm->SetMaxNIter(_nIter);
	_cleanAlgorithm->SetThreshold(_threshold);
	_cleanAlgorithm->SetSubtractionGain(_gain);
	_cleanAlgorithm->SetStopGain(_mGain);
	_cleanAlgorithm->SetAllowNegativeComponents(_allowNegative);
	_cleanAlgorithm->SetStopOnNegativeComponents(_stopOnNegative);
	_cleanAlgorithm->SetResizePSF(_smallPSF);
	if(!_cleanAreasFilename.empty())
	{
		_cleanAreas.reset(new AreaSet());
		AreaParser parser;
		std::ifstream caFile(_cleanAreasFilename.c_str());
		parser.Parse(*_cleanAreas, caFile);
		_cleanAreas->SetImageProperties(_pixelScaleX, _pixelScaleY, _inversionAlgorithm->PhaseCentreRA(), _inversionAlgorithm->PhaseCentreDec(), _imgWidth, _imgHeight);
		_cleanAlgorithm->SetCleanAreas(*_cleanAreas);
	}
}

void WSClean::prepareInversionAlgorithm(PolarizationEnum polarization)
{
	_inversionAlgorithm->SetImageWidth(_imgWidth);
	_inversionAlgorithm->SetImageHeight(_imgHeight);
	_inversionAlgorithm->SetPixelSizeX(_pixelScaleX);
	_inversionAlgorithm->SetPixelSizeY(_pixelScaleY);
	if(_nWLayers != 0)
		_inversionAlgorithm->SetWGridSize(_nWLayers);
	else
		_inversionAlgorithm->SetNoWGridSize();
	_inversionAlgorithm->SetAntialiasingKernelSize(_antialiasingKernelSize);
	_inversionAlgorithm->SetOverSamplingFactor(_overSamplingFactor);
	_inversionAlgorithm->SetPolarization(polarization);
	_inversionAlgorithm->SetIsComplex(polarization == Polarization::XY || polarization == Polarization::YX);
	_inversionAlgorithm->SetDataColumnName(_columnName);
	_inversionAlgorithm->SetWeighting(_weightMode);
	_inversionAlgorithm->SetSelection(_currentPartSelection);
	_inversionAlgorithm->SetWLimit(_wLimit/100.0);
	_inversionAlgorithm->SetSmallInversion(_smallInversion);
}

void WSClean::Run()
{
	// If no column specified, determine column to use
	if(_columnName.empty())
	{
		casa::MeasurementSet ms(_filenames.front());
		bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
		if(hasCorrected) {
			std::cout << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
			_columnName = "CORRECTED_DATA";
		} else {
			std::cout << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
			_columnName= "DATA";
		}
	}

	_joinedPolarizationCleaning = Polarization::Has4Polarizations(_polarizations);
	
	_doReorder = ((_channelsOut != 1) || _forceReorder) && !_forceNoReorder;
	
	if(_doReorder)
	{
		for(std::vector<std::string>::const_iterator i=_filenames.begin(); i != _filenames.end(); ++i)
		{
			_partitionedMSHandles.push_back(PartitionedMS::Partition(*i, _channelsOut, _globalSelection, _columnName, true, _mGain != 1.0, _polarizations));
		}
	}
	
	_firstMSBand = BandData(casa::MeasurementSet(_filenames[0]).spectralWindow());
	
	for(size_t outChannelIndex=0; outChannelIndex!=_channelsOut; ++outChannelIndex)
	{
		runChannel(outChannelIndex);
	}
}

void WSClean::runChannel(size_t outChannelIndex)
{
	MSSelection partSelection = _globalSelection;
	_inversionAlgorithm.reset(new WSInversion(&_imageAllocator, _memFraction));
	static_cast<WSInversion&>(*_inversionAlgorithm).SetGridMode(_gridMode);
	
	if(_doReorder)
	{
		size_t startCh, endCh;
		if(_globalSelection.HasChannelRange())
		{
			startCh = _globalSelection.ChannelRangeStart();
			endCh = _globalSelection.ChannelRangeEnd();
		}
		else {
			startCh = 0;
			endCh = _firstMSBand.ChannelCount();
		}
		size_t newStart = startCh + (endCh - startCh) * outChannelIndex / _channelsOut;
		size_t newEnd = startCh + (endCh - startCh) * (outChannelIndex+1) / _channelsOut;
		partSelection.SetChannelRange(newStart, newEnd);
	}
		
	const std::string rootPrefix = _prefixName;
	if(_channelsOut != 1)
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << rootPrefix << '-';
		if(outChannelIndex < 1000) partPrefixNameStr << '0';
		if(outChannelIndex < 100) partPrefixNameStr << '0';
		if(outChannelIndex < 10) partPrefixNameStr << '0';
		partPrefixNameStr << outChannelIndex;
		_prefixName = partPrefixNameStr.str();
	}
		
	for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
	{
		runPolarizationStart(outChannelIndex, *curPol);
	}
	
	initializeCleanAlgorithm();

	initFitsWriter(_fitsWriter);
	setCleanParameters(_fitsWriter, *_cleanAlgorithm);
	updateCleanParameters(_fitsWriter, 0, 0);
		
	if(_nIter > 0)
	{
		// Start major cleaning loop
		size_t majorIterationNr = 1;
		bool reachedMajorThreshold = false;
		do {
			performClean(reachedMajorThreshold, majorIterationNr);
			
			if(_mGain != 1.0)
			{
				for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
				{
					prepareInversionAlgorithm(*curPol);
					
					// TODO handle imaginary
					initializeCurMSProviders(outChannelIndex, *curPol);
					predict(*curPol);
					
					imageMainNonFirst(*curPol);
					clearCurMSProviders();
				
					if(!reachedMajorThreshold)
					{
						// This was the final major iteration: save results
						if(!(*curPol == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
						{
							double* residualImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
							_residualImages.Load(residualImage, *curPol, false);
							std::cout << "Writing residual image... " << std::flush;
							_fitsWriter.Write(polPrefix(*curPol, false) + "-residual.fits", residualImage);
							_imageAllocator.Free(residualImage);
							std::cout << "DONE\n";
						}
					}
				}
				
				++majorIterationNr;
			}
			
		} while(reachedMajorThreshold);
		
		std::cout << majorIterationNr << " major iterations were performed.\n";
	}
	
	for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
	{
		if(!(*curPol == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
		{
			std::vector<bool> imaginaryStates(1, false);
			if(*curPol == Polarization::XY || *curPol == Polarization::YX)
				imaginaryStates.push_back(true);
			
			for(std::vector<bool>::const_iterator isImaginary=imaginaryStates.begin(); isImaginary!=imaginaryStates.end(); ++isImaginary)
			{
				Model model;
				if(!_addModelFilename.empty())
				{
					std::cout << "Reading model from " << _addModelFilename << "... " << std::flush;
					model = Model(_addModelFilename.c_str());
				}
				double* modelImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
				_modelImages.Load(modelImage, *curPol, *isImaginary);
				// A model cannot hold instrumental pols (xx/xy/yx/yy), hence always use Stokes I here
				CleanAlgorithm::GetModelFromImage(model, modelImage, _imgWidth, _imgHeight, _fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, 0.0, _fitsWriter.Frequency(), Polarization::StokesI);
				_imageAllocator.Free(modelImage);
				
				if(!_saveModelFilename.empty())
				{
					std::cout << "Saving model to " << _saveModelFilename << "... " << std::flush;
					model.Save(_saveModelFilename.c_str());
				}
			
				ModelRenderer renderer(_fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY);
				double
					freqLow = _fitsWriter.Frequency() - _fitsWriter.Bandwidth()*0.5,
					freqHigh = _fitsWriter.Frequency() + _fitsWriter.Bandwidth()*0.5;
				double* restoredImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
				_residualImages.Load(restoredImage, *curPol, *isImaginary);
				std::cout << "Rendering " << model.SourceCount() << " sources to restored image... " << std::flush;
				renderer.Restore(restoredImage, _imgWidth, _imgHeight, model, _fitsWriter.BeamSizeMajorAxis(), freqLow, freqHigh, Polarization::StokesI);
				std::cout << "DONE\n";
				
				std::cout << "Writing restored image... " << std::flush;
				_fitsWriter.Write(polPrefix(*curPol, *isImaginary) + "-image.fits", restoredImage);
				std::cout << "DONE\n";
				_imageAllocator.Free(restoredImage);
			}
		}
	}
	
	_imageAllocator.ReportStatistics();
	std::cout << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", cleaning: " << _cleaningWatch.ToString() << '\n';
	
	_prefixName = rootPrefix;
	
	// Needs to be desctructed before image allocator, or image allocator will report error caused by leaked memory
	_inversionAlgorithm.reset();
}

void WSClean::initializeCurMSProviders(size_t outChannelIndex, PolarizationEnum polarization)
{
	_inversionAlgorithm->ClearMeasurementSetList();
	for(size_t i=0; i != _filenames.size(); ++i)
	{
		MSProvider* msProvider;
		if(_doReorder)
			msProvider = new PartitionedMS(_partitionedMSHandles[i], outChannelIndex, polarization);
		else
			msProvider = new ContiguousMS(_filenames[i], _columnName, _currentPartSelection, polarization, _mGain != 1.0);
		_inversionAlgorithm->AddMeasurementSet(msProvider);
		_currentPolMSes.push_back(msProvider);
	}
}

void WSClean::clearCurMSProviders()
{
	for(std::vector<MSProvider*>::iterator i=_currentPolMSes.begin(); i != _currentPolMSes.end(); ++i)
		delete *i;
	_currentPolMSes.clear();
}

void WSClean::runPolarizationStart(size_t outChannelIndex, PolarizationEnum polarization)
{
	initializeCurMSProviders(outChannelIndex, polarization);
	
	if(polarization == *_polarizations.begin())
		initializeImageWeights(_currentPartSelection);
	
	prepareInversionAlgorithm(polarization);
	
	const bool firstBeforePSF = _isFirstInversion;

	if((_nIter > 0 || _makePSF) && polarization == *_polarizations.begin())
		imagePSF();
	
	initFitsWriter(_fitsWriter);
	_modelImages.Initialize(_fitsWriter, _polarizations.size(), _prefixName + "-model", _imageAllocator);
	_residualImages.Initialize(_fitsWriter, _polarizations.size(), _prefixName + "-residual", _imageAllocator);
	
	imageMainFirst(polarization);
	
	if(firstBeforePSF && _inversionAlgorithm->HasGriddingCorrectionImage())
		imageGridding();
	
	_isFirstInversion = false;
	
	if(!(polarization == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
	{
		double* modelImage = _imageAllocator.Allocate(_imgWidth * _imgHeight);
		memset(modelImage, 0, _imgWidth * _imgHeight * sizeof(double));
		_modelImages.Store(modelImage, polarization, false);
		if(Polarization::IsComplex(polarization))
			_modelImages.Store(modelImage, polarization, true);
		_imageAllocator.Free(modelImage);
	}
	
	std::cout << "Writing dirty image... " << std::flush;
	_fitsWriter.Write(polPrefix(polarization, false) + "-dirty.fits", _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		_fitsWriter.Write(polPrefix(polarization, true) + "-dirty.fits", _inversionAlgorithm->ImageImaginaryResult());
	std::cout << "DONE\n";
	
	clearCurMSProviders();
}

void WSClean::performClean(bool& reachedMajorThreshold, size_t majorIterationNr)
{
	std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
	
	if(!_joinedPolarizationCleaning)
		performSimpleClean(reachedMajorThreshold, majorIterationNr);
	else
		performJoinedPolClean(reachedMajorThreshold, majorIterationNr);
}

void WSClean::performSimpleClean(bool& reachedMajorThreshold, size_t majorIterationNr)
{
	PolarizationEnum polarization = *_polarizations.begin();
	double
		*residualImage = _imageAllocator.Allocate(_imgWidth*_imgHeight),
		*modelImage = _imageAllocator.Allocate(_imgWidth*_imgHeight),
		*psfImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
	_residualImages.Load(residualImage, polarization, false);
	_modelImages.Load(modelImage, polarization, false);
	_psfImages.Load(psfImage, polarization, false);
	
	_cleaningWatch.Start();
	static_cast<SimpleClean&>(*_cleanAlgorithm).ExecuteMajorIteration(residualImage, modelImage, psfImage, _imgWidth, _imgHeight, reachedMajorThreshold);
	_cleaningWatch.Pause();
	
	_modelImages.Store(modelImage, polarization, false);
	
	updateCleanParameters(_fitsWriter, _cleanAlgorithm->IterationNumber(), majorIterationNr);
	
	if(majorIterationNr == 1)
	{
		if(_mGain == 1.0)
		{
			std::cout << "Writing residual image... " << std::flush;
			_fitsWriter.Write(polPrefix(polarization, false) + "-residual.fits", residualImage);
		}
		else {
			std::cout << "Writing first iteration image... " << std::flush;
			_fitsWriter.Write(polPrefix(polarization, false) + "-first-residual.fits", residualImage);
		}
		std::cout << "DONE\n";
	}
	if(!reachedMajorThreshold)
	{
		std::cout << "Writing model image... " << std::flush;
		_fitsWriter.Write(polPrefix(polarization, false) + "-model.fits", modelImage);
		std::cout << "DONE\n";
	}
	
	_imageAllocator.Free(residualImage);
	_imageAllocator.Free(modelImage);
	_imageAllocator.Free(psfImage);
}

void WSClean::performJoinedPolClean(bool& reachedMajorThreshold, size_t majorIterationNr)
{
	JoinedPolClean::ImageSet
		modelSet(_imgWidth*_imgHeight, _imageAllocator),
		residualSet(_imgWidth*_imgHeight, _imageAllocator);
		
	double* psfImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
	_psfImages.Load(psfImage, *_polarizations.begin(), false);
	modelSet.Load(_modelImages);
	residualSet.Load(_residualImages);

	_cleaningWatch.Start();
	static_cast<JoinedPolClean&>(*_cleanAlgorithm).ExecuteMajorIteration(residualSet, modelSet, psfImage, _imgWidth, _imgHeight, reachedMajorThreshold);
	_cleaningWatch.Pause();
	
	_imageAllocator.Free(psfImage);
	modelSet.Store(_modelImages);
	residualSet.Store(_residualImages);
	
	updateCleanParameters(_fitsWriter, _cleanAlgorithm->IterationNumber(), majorIterationNr);
	
	PolarizationEnum pols[4] = { Polarization::XX, Polarization::XY, Polarization::XY, Polarization::YY };
	for(size_t i=0; i!=4; ++i)
	{
		PolarizationEnum polarization = pols[i];
		bool isImaginary = (i==2);
		
		if(majorIterationNr == 1)
		{
			if(_mGain == 1.0)
			{
				std::cout << "Writing residual image... " << std::flush;
				_fitsWriter.Write(polPrefix(polarization, isImaginary) + "-residual.fits", residualSet.Get(i));
			}
			else {
				std::cout << "Writing first iteration image... " << std::flush;
				_fitsWriter.Write(polPrefix(polarization, isImaginary) + "-first-residual.fits", residualSet.Get(i));
			}
			std::cout << "DONE\n";
		}
		if(!reachedMajorThreshold)
		{
			std::cout << "Writing model image... " << std::flush;
			_fitsWriter.Write(polPrefix(polarization, isImaginary) + "-model.fits", modelSet.Get(i));
			std::cout << "DONE\n";
		}
	}
}
