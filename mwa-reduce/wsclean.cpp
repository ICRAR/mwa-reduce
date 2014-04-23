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
#include "cleanalgorithms/joinedclean.h"
#include "cleanalgorithms/simpleclean.h"

#include "msprovider/contiguousms.h"

#include "parser/areaparser.h"

#include "uvector.h"

#include <iostream>
#include <memory>

std::string commandLine;

WSClean::WSClean() :
	_imgWidth(2048), _imgHeight(2048), _channelsOut(1),
	_pixelScaleX(0.01 * M_PI / 180.0), _pixelScaleY(0.01 * M_PI / 180.0),
	_threshold(0.0), _gain(0.1), _mGain(1.0), _manualBeamSize(0.0), _memFraction(1.0), _absMemLimit(0.0), _wLimit(0.0),
	_nWLayers(0), _nIter(0), _antialiasingKernelSize(7), _overSamplingFactor(63),
	_globalSelection(),
	_columnName(), _addModelFilename(), _saveModelFilename(), _cleanAreasFilename(),
	_polarizations(),
	_weightMode(WeightMode::UniformWeighted),
	_prefixName("wsclean"),
	_allowNegative(true), _smallPSF(false), _smallInversion(false), _stopOnNegative(false), _makePSF(false),
	_forceReorder(false), _forceNoReorder(false), _joinedPolarizationCleaning(false), _joinedFrequencyCleaning(false),
	_gridMode(LayeredImager::KaiserBessel),
	_filenames(),
	_commandLine(),
	_inversionWatch(false), _predictingWatch(false), _cleaningWatch(false),
	_isFirstInversion(true), _doReorder(false),
	_majorIterationNr(0)
{
	_polarizations.insert(Polarization::StokesI);
}

WSClean::~WSClean()
{
	freeCleanAlgorithms();
}

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
	
	writer.SetExtraKeyword("WSCIMGWG", _inversionAlgorithm->ImageWeight());
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

void WSClean::imagePSF(size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Constructing PSF ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoImagePSF(true);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
		
	CleanAlgorithm::RemoveNaNsInPSF(_inversionAlgorithm->ImageRealResult(), _imgWidth, _imgHeight);
	initFitsWriter(_fitsWriter);
	_psfImages.SetFitsWriter(_fitsWriter);
	_psfImages.Store(_inversionAlgorithm->ImageRealResult(), *_polarizations.begin(), joinedChannelIndex, false);
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
	fitsWriter.SetImageDimensions(_inversionAlgorithm->ActualInversionWidth(), _inversionAlgorithm->ActualInversionHeight());
	fitsWriter.Write(_prefixName + "-gridding.fits", &gridding[0]);
	_imageAllocator.Free(gridding);
	std::cout << "DONE\n";
}

void WSClean::imageMainFirst(PolarizationEnum polarization, size_t joinedChannelIndex)
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
	
	storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, true, _inversionAlgorithm->ImageImaginaryResult());
}

void WSClean::imageMainNonFirst(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Constructing image ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoSubtractModel(true);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	
	storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, true, _inversionAlgorithm->ImageImaginaryResult());
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image)
{
	if(polarization == Polarization::YX && _polarizations.count(Polarization::XY)!=0)
	{
		std::cout << "Adding XY and YX together...\n";
		double
			*xyImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
		dest.Load(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
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
		dest.Store(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
		_imageAllocator.Free(xyImage);
	}
	else {
		dest.Store(image, polarization, joinedChannelIndex, isImaginary);
	}
}

void WSClean::predict(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Converting model image to visibilities ==\n";
	const size_t size = _imgWidth*_imgHeight;
	double
		*modelImageReal = _imageAllocator.Allocate(size),
		*modelImageImaginary = 0;
		
	if(polarization == Polarization::YX)
	{
		_modelImages.Load(modelImageReal, Polarization::XY, joinedChannelIndex, false);
		modelImageImaginary = _imageAllocator.Allocate(size);
		_modelImages.Load(modelImageImaginary, Polarization::XY, joinedChannelIndex, true);
		for(size_t i=0; i!=size; ++i)
			modelImageImaginary[i] = -modelImageImaginary[i];
	}
	else {
		_modelImages.Load(modelImageReal, polarization, joinedChannelIndex, false);
		if(Polarization::IsComplex(polarization))
		{
			modelImageImaginary = _imageAllocator.Allocate(size);
			_modelImages.Load(modelImageImaginary, polarization, joinedChannelIndex, true);
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
		if(_mfsWeighting)
		{
			std::cout << "Reusing MFS weights for " << _weightMode.ToString() << " weighting.\n";
			_inversionAlgorithm->SetPrecalculatedWeightInfo(_imageWeights.get());
		}
		else {
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
}

void WSClean::initializeMFSImageWeights()
{
	if(_weightMode.RequiresGridding())
	{
		std::cout << "Precalculating MFS weights for " << _weightMode.ToString() << " weighting...\n";
		_imageWeights.reset(new ImageWeights(_imgWidth, _imgHeight, _pixelScaleX, _pixelScaleY, _weightMode.SuperWeight()));
		for(size_t i=0; i!=_partitionedMSHandles.size(); ++i)
		{
			if(_doReorder)
			{
				for(size_t ch=0; ch!=_channelsOut; ++ch)
				{
					MSSelection partSelection(_globalSelection);
					selectChannels(partSelection, ch, _channelsOut);
					PartitionedMS msProvider(_partitionedMSHandles[i], ch, *_polarizations.begin());
					_imageWeights->Grid(msProvider, _weightMode, partSelection);
				}
			}
			else {
				ContiguousMS msProvider(_filenames[i], _columnName, _globalSelection, *_polarizations.begin(), _mGain != 1.0);
				_imageWeights->Grid(msProvider, _weightMode, _globalSelection);
				std::cout << '.' << std::flush;
			}
		}
	}
}

void WSClean::freeCleanAlgorithms()
{
	for(std::vector<CleanAlgorithm*>::iterator caPtr=_cleanAlgorithms.begin(); caPtr!=_cleanAlgorithms.end(); ++caPtr)
		delete *caPtr;
	_cleanAlgorithms.clear();
}

void WSClean::initializeCleanAlgorithm()
{
	freeCleanAlgorithms();
	size_t count = 0;
	if(_joinedPolarizationCleaning)
	{
		_cleanAlgorithms.resize(1);
		_cleanAlgorithms[0] = new JoinedClean<>();
		count = 1;
	}
	else {
		_cleanAlgorithms.resize(_polarizations.size());
		for(size_t p=0; p!=_polarizations.size(); ++p)
			_cleanAlgorithms[p] = new SimpleClean();
		count = _polarizations.size();
	}
	for(size_t p=0; p!=count; ++p)
	{
		_cleanAlgorithms[p]->SetMaxNIter(_nIter);
		_cleanAlgorithms[p]->SetThreshold(_threshold);
		_cleanAlgorithms[p]->SetSubtractionGain(_gain);
		_cleanAlgorithms[p]->SetStopGain(_mGain);
		_cleanAlgorithms[p]->SetAllowNegativeComponents(_allowNegative);
		_cleanAlgorithms[p]->SetStopOnNegativeComponents(_stopOnNegative);
		_cleanAlgorithms[p]->SetResizePSF(_smallPSF);
		if(!_cleanAreasFilename.empty())
		{
			_cleanAreas.reset(new AreaSet());
			AreaParser parser;
			std::ifstream caFile(_cleanAreasFilename.c_str());
			parser.Parse(*_cleanAreas, caFile);
			_cleanAreas->SetImageProperties(_pixelScaleX, _pixelScaleY, _inversionAlgorithm->PhaseCentreRA(), _inversionAlgorithm->PhaseCentreDec(), _imgWidth, _imgHeight);
			_cleanAlgorithms[p]->SetCleanAreas(*_cleanAreas);
		}
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

	if(_joinedPolarizationCleaning && !Polarization::HasFullPolarization(_polarizations))
		throw std::runtime_error("Joined polarization cleaning requested, but not all 4 polarizations are imaged");
	
	_doReorder = ((_channelsOut != 1) || (_polarizations.size()>=4) || _forceReorder) && !_forceNoReorder;
	
	if(_doReorder)
	{
		for(std::vector<std::string>::const_iterator i=_filenames.begin(); i != _filenames.end(); ++i)
		{
			_partitionedMSHandles.push_back(PartitionedMS::Partition(*i, _channelsOut, _globalSelection, _columnName, true, _mGain != 1.0, _polarizations));
		}
	}
	
	_firstMSBand = BandData(casa::MeasurementSet(_filenames[0]).spectralWindow());
	_weightPerChannel.assign(_channelsOut, 0);
	
	if(_mfsWeighting)
		initializeMFSImageWeights();
	
	if(_joinedFrequencyCleaning)
	{
		// run all frequencies
		runIndependentChannel(0);
	}
	else {
		for(size_t outChannelIndex=0; outChannelIndex!=_channelsOut; ++outChannelIndex)
		{
			runIndependentChannel(outChannelIndex);
		}
	}
	
	if(_channelsOut > 1)
	{
		for(std::set<PolarizationEnum>::const_iterator pol=_polarizations.begin(); pol!=_polarizations.end(); ++pol)
		{
			if(!(*pol == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
			{
				makeMFSImage("image.fits", *pol, false);
				if(_nIter > 0)
				{
					makeMFSImage("residual.fits", *pol, false);
					makeMFSImage("model.fits", *pol, false);
				}
				if(Polarization::IsComplex(*pol))
				{
					makeMFSImage("image.fits", *pol, true);
					if(_nIter > 0)
					{
					  makeMFSImage("residual.fits", *pol, true);
					  makeMFSImage("model.fits", *pol, true);
					}
				}
			}
		}
	}
}

void WSClean::selectChannels(MSSelection& selection, size_t outChannelIndex, size_t channelsOut)
{
	size_t startCh, endCh;
	if(selection.HasChannelRange())
	{
		startCh = selection.ChannelRangeStart();
		endCh = selection.ChannelRangeEnd();
	}
	else {
		startCh = 0;
		endCh = _firstMSBand.ChannelCount();
	}
	size_t newStart = startCh + (endCh - startCh) * outChannelIndex / channelsOut;
	size_t newEnd = startCh + (endCh - startCh) * (outChannelIndex+1) / channelsOut;
	selection.SetChannelRange(newStart, newEnd);
}

void WSClean::runIndependentChannel(size_t outChannelIndex)
{
	_inversionAlgorithm.reset(new WSInversion(&_imageAllocator, _memFraction, _absMemLimit));
	static_cast<WSInversion&>(*_inversionAlgorithm).SetGridMode(_gridMode);
	
	size_t  joinedChannelsOut;
	if(_joinedFrequencyCleaning)
		joinedChannelsOut = _channelsOut;
	else
		joinedChannelsOut = 1;
	
	_currentPartSelection = _globalSelection;
	if(!_joinedFrequencyCleaning)
		selectChannels(_currentPartSelection, outChannelIndex, _channelsOut);
	
	_modelImages.Initialize(_fitsWriter, _polarizations.size(), joinedChannelsOut, _prefixName + "-model", _imageAllocator);
	_residualImages.Initialize(_fitsWriter, _polarizations.size(), joinedChannelsOut, _prefixName + "-residual", _imageAllocator);
	_psfImages.Initialize(_fitsWriter, 1, joinedChannelsOut, _prefixName + "-psf", _imageAllocator);
	
	const std::string rootPrefix = _prefixName;
		
	for(size_t ch=0; ch!=joinedChannelsOut; ++ch)
	{
		const MSSelection selectionBefore = _currentPartSelection;
		size_t currentChannelIndex;
		if(_joinedFrequencyCleaning)
		{
			selectChannels(_currentPartSelection, ch, _channelsOut);
			currentChannelIndex = ch;
		}
		else {
			currentChannelIndex = outChannelIndex;
		}
		for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
		{
			runFirstInversion(currentChannelIndex, *curPol, ch);
		}
		_currentPartSelection = selectionBefore;
	}
	
	initializeCleanAlgorithm();

	initFitsWriter(_fitsWriter);
	setCleanParameters(_fitsWriter, *_cleanAlgorithms[0]);
	updateCleanParameters(_fitsWriter, 0, 0);
		
	if(_nIter > 0)
	{
		// Start major cleaning loop
		_majorIterationNr = 1;
		bool reachedMajorThreshold = false;
		do {
			performClean(outChannelIndex, reachedMajorThreshold, _majorIterationNr);
			
			if(_mGain != 1.0)
			{
				for(size_t ch=0; ch!=joinedChannelsOut; ++ch)
				{
					const MSSelection selectionBefore = _currentPartSelection;
					size_t currentChannelIndex;
					if(_joinedFrequencyCleaning)
					{
						selectChannels(_currentPartSelection, ch, _channelsOut);
						currentChannelIndex = ch;
					}
					else {
						currentChannelIndex = outChannelIndex;
					}
					for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
					{
						prepareInversionAlgorithm(*curPol);
						
						initializeCurMSProviders(currentChannelIndex, *curPol);
						if(*curPol == *_polarizations.begin())
							initializeImageWeights(_currentPartSelection);
	
						predict(*curPol, ch);
						
						imageMainNonFirst(*curPol, ch);
						clearCurMSProviders();
					
						if(!reachedMajorThreshold)
						{
							// This was the final major iteration: save results
							if(!(*curPol == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
							{
								double* residualImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
								_residualImages.Load(residualImage, *curPol, ch, false);
								std::cout << "Writing residual image... " << std::flush;
								writeFits("residual.fits", residualImage, *curPol, currentChannelIndex, false);
								if(Polarization::IsComplex(*curPol))
									writeFits("residual.fits", residualImage, *curPol, currentChannelIndex, true);
								_imageAllocator.Free(residualImage);
								std::cout << "DONE\n";
							}
						}
					} // end of polarization loop
					_currentPartSelection = selectionBefore;
				} // end of joined channels loop
				
				++_majorIterationNr;
			}
			
		} while(reachedMajorThreshold);
		
		std::cout << _majorIterationNr << " major iterations were performed.\n";
	}
	
	
	// Restore model to residuals and save all images
	for(size_t ch=0; ch!=joinedChannelsOut; ++ch)
	{
		size_t currentChannelIndex = _joinedFrequencyCleaning ? ch : outChannelIndex;
		MSSelection curSelection(_globalSelection);
		selectChannels(curSelection, currentChannelIndex, _channelsOut);
		const BandData curBand(_firstMSBand, curSelection.ChannelRangeStart(), curSelection.ChannelRangeEnd());
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
					_modelImages.Load(modelImage, *curPol, currentChannelIndex, *isImaginary);
					// A model cannot hold instrumental pols (xx/xy/yx/yy), hence always use Stokes I here
					CleanAlgorithm::GetModelFromImage(model, modelImage, _imgWidth, _imgHeight, _fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, 0.0, _fitsWriter.Frequency(), Polarization::StokesI);
					_imageAllocator.Free(modelImage);
					
					if(!_saveModelFilename.empty())
					{
						std::cout << "Saving model to " << _saveModelFilename << "... " << std::flush;
						model.Save(_saveModelFilename.c_str());
					}
				
					ModelRenderer renderer(_fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY);
					double freqLow = curBand.LowestFrequency(), freqHigh = curBand.HighestFrequency();
					double* restoredImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
					_residualImages.Load(restoredImage, *curPol, currentChannelIndex, *isImaginary);
					std::cout << "Rendering " << model.SourceCount() << " sources to restored image... " << std::flush;
					renderer.Restore(restoredImage, _imgWidth, _imgHeight, model, _fitsWriter.BeamSizeMajorAxis(), freqLow, freqHigh, Polarization::StokesI);
					std::cout << "DONE\n";
					
					std::cout << "Writing restored image... " << std::flush;
					writeFits("image.fits", restoredImage, *curPol, currentChannelIndex, *isImaginary);
					std::cout << "DONE\n";
					_imageAllocator.Free(restoredImage);
				}
			}
		}
	}
	
	_imageAllocator.ReportStatistics();
	std::cout << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", cleaning: " << _cleaningWatch.ToString() << '\n';
	
	_prefixName = rootPrefix;
	
	// Needs to be desctructed before image allocator, or image allocator will report error caused by leaked memory
	_inversionAlgorithm.reset();
}

void WSClean::initializeCurMSProviders(size_t currentChannelIndex, PolarizationEnum polarization)
{
	_inversionAlgorithm->ClearMeasurementSetList();
	for(size_t i=0; i != _filenames.size(); ++i)
	{
		MSProvider* msProvider;
		if(_doReorder)
			msProvider = new PartitionedMS(_partitionedMSHandles[i], currentChannelIndex, polarization);
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

void WSClean::runFirstInversion(size_t currentChannelIndex, PolarizationEnum polarization, size_t joinedChannelIndex)
{
	initializeCurMSProviders(currentChannelIndex, polarization);
	
	if(polarization == *_polarizations.begin())
		initializeImageWeights(_currentPartSelection);
	
	prepareInversionAlgorithm(polarization);
	
	const bool firstBeforePSF = _isFirstInversion;

	if((_nIter > 0 || _makePSF) && polarization == *_polarizations.begin())
		imagePSF(joinedChannelIndex);
	
	initFitsWriter(_fitsWriter);
	_modelImages.SetFitsWriter(_fitsWriter);
	_residualImages.SetFitsWriter(_fitsWriter);
	
	imageMainFirst(polarization, joinedChannelIndex);
	
	_weightPerChannel[currentChannelIndex] = _inversionAlgorithm->ImageWeight();
	
	if(firstBeforePSF && _inversionAlgorithm->HasGriddingCorrectionImage())
		imageGridding();
	
	_isFirstInversion = false;
	
	if(!(polarization == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
	{
		double* modelImage = _imageAllocator.Allocate(_imgWidth * _imgHeight);
		memset(modelImage, 0, _imgWidth * _imgHeight * sizeof(double));
		_modelImages.Store(modelImage, polarization, joinedChannelIndex, false);
		if(Polarization::IsComplex(polarization))
			_modelImages.Store(modelImage, polarization, joinedChannelIndex, true);
		_imageAllocator.Free(modelImage);
	}
	
	std::cout << "Writing dirty image... " << std::flush;
	writeFits("dirty.fits", _inversionAlgorithm->ImageRealResult(), polarization, currentChannelIndex, false);
	if(Polarization::IsComplex(polarization))
		writeFits("dirty.fits", _inversionAlgorithm->ImageImaginaryResult(), polarization, currentChannelIndex, true);
	std::cout << "DONE\n";
	
	clearCurMSProviders();
}

void WSClean::performClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr)
{
	std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
	
	if(_joinedFrequencyCleaning)
	{
		if(!_joinedPolarizationCleaning)
			throw std::runtime_error("Can only joinedly clean frequencies when cleaning all polarizations simultaneously");
		performJoinedPolFreqClean(reachedMajorThreshold, majorIterationNr);
	}
	else if(_joinedPolarizationCleaning)
		performJoinedPolClean(currentChannelIndex, reachedMajorThreshold, majorIterationNr);
	else {
		std::set<PolarizationEnum>::const_iterator pol = _polarizations.begin();
		for(size_t pIndex=0; pIndex!=_polarizations.size(); ++pIndex)
		{
			performSimpleClean(*_cleanAlgorithms[pIndex], currentChannelIndex, reachedMajorThreshold, majorIterationNr, *pol);
			++pol;
		}
	}
}

void WSClean::performSimpleClean(CleanAlgorithm& cleanAlgorithm, size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr, PolarizationEnum polarization)
{
	double
		*residualImage = _imageAllocator.Allocate(_imgWidth*_imgHeight),
		*modelImage = _imageAllocator.Allocate(_imgWidth*_imgHeight),
		*psfImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
	_residualImages.Load(residualImage, polarization, 0, false);
	_modelImages.Load(modelImage, polarization, 0, false);
	_psfImages.Load(psfImage, polarization, 0, false);
	
	_cleaningWatch.Start();
	static_cast<SimpleClean&>(cleanAlgorithm).ExecuteMajorIteration(residualImage, modelImage, psfImage, _imgWidth, _imgHeight, reachedMajorThreshold);
	_cleaningWatch.Pause();
	
	_modelImages.Store(modelImage, polarization, 0, false);
	
	updateCleanParameters(_fitsWriter, cleanAlgorithm.IterationNumber(), majorIterationNr);
	
	if(majorIterationNr == 1)
	{
		if(_mGain == 1.0)
		{
			std::cout << "Writing residual image... " << std::flush;
			writeFits("residual.fits", residualImage, polarization, currentChannelIndex, false);
		}
		else {
			std::cout << "Writing first iteration image... " << std::flush;
			writeFits("first-residual.fits", residualImage, polarization, currentChannelIndex, false);
		}
		std::cout << "DONE\n";
	}
	if(!reachedMajorThreshold)
	{
		std::cout << "Writing model image... " << std::flush;
		writeFits("model.fits", modelImage, polarization, currentChannelIndex, false);
		std::cout << "DONE\n";
	}
	
	_imageAllocator.Free(residualImage);
	_imageAllocator.Free(modelImage);
	_imageAllocator.Free(psfImage);
}

void WSClean::performJoinedPolClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr)
{
	JoinedClean<>::ImageSet
		modelSet(_imgWidth*_imgHeight, _imageAllocator),
		residualSet(_imgWidth*_imgHeight, _imageAllocator);
	
	std::vector<double*> psfImages;
	double* psfImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
	_psfImages.Load(psfImage, *_polarizations.begin(), 0, false);
	psfImages.push_back(psfImage);
	
	bool hasStokesPols = Polarization::HasFullStokesPolarization(_polarizations);
	if(hasStokesPols)
	{
		modelSet.LoadStokes(_modelImages, 0);
		residualSet.LoadStokes(_residualImages, 0);
	}
	else {
		modelSet.LoadLinear(_modelImages, 0);
		residualSet.LoadLinear(_residualImages, 0);
	}

	_cleaningWatch.Start();
	static_cast<JoinedClean<>&>(*_cleanAlgorithms[0]).ExecuteMajorIteration(residualSet, modelSet, psfImages, _imgWidth, _imgHeight, reachedMajorThreshold);
	_cleaningWatch.Pause();
	
	_imageAllocator.Free(psfImage);
	if(hasStokesPols)
	{
		modelSet.StoreStokes(_modelImages, 0);
		residualSet.StoreStokes(_residualImages, 0);
	}
	else {
		modelSet.StoreLinear(_modelImages, 0);
		residualSet.StoreLinear(_residualImages, 0);
	}
	
	updateCleanParameters(_fitsWriter, _cleanAlgorithms[0]->IterationNumber(), majorIterationNr);
	
	const PolarizationEnum
		linPols[4] = { Polarization::XX, Polarization::XY, Polarization::XY, Polarization::YY },
		stokesPols[4] = { Polarization::StokesI, Polarization::StokesQ, Polarization::StokesU, Polarization::StokesV },
		*pols = hasStokesPols ? stokesPols : linPols;
	for(size_t i=0; i!=4; ++i)
	{
		PolarizationEnum polarization = pols[i];
		bool isImaginary = (i==2) && !hasStokesPols;
		
		if(majorIterationNr == 1)
		{
			if(_mGain == 1.0)
			{
				std::cout << "Writing residual image... " << std::flush;
				writeFits("residual.fits", residualSet.GetImage(i), polarization, currentChannelIndex, isImaginary);
			}
			else {
				std::cout << "Writing first iteration image... " << std::flush;
				writeFits("first-residual.fits", residualSet.GetImage(i), polarization, currentChannelIndex, isImaginary);
			}
			std::cout << "DONE\n";
		}
		if(!reachedMajorThreshold)
		{
			std::cout << "Writing model image... " << std::flush;
			writeFits("model.fits", modelSet.GetImage(i), polarization, currentChannelIndex, isImaginary);
			std::cout << "DONE\n";
		}
	}
}

void WSClean::performJoinedPolFreqClean(bool& reachedMajorThreshold, size_t majorIterationNr)
{
	JoinedClean<joined_pol_clean::MultiImageSet>::ImageSet
		modelSet(_imgWidth*_imgHeight, _channelsOut, _imageAllocator),
		residualSet(_imgWidth*_imgHeight, _channelsOut, _imageAllocator);
	
	bool hasStokesPols = Polarization::HasFullStokesPolarization(_polarizations);
	
	std::vector<double*> psfImages;
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		double* psfImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
		_psfImages.Load(psfImage, *_polarizations.begin(), ch, false);
		psfImages.push_back(psfImage);
		
		if(hasStokesPols)
		{
			modelSet.LoadStokes(_modelImages, ch);
			residualSet.LoadStokes(_residualImages, ch);
		}
		else {
			modelSet.LoadLinear(_modelImages, ch);
			residualSet.LoadLinear(_residualImages, ch);
		}
	}

	_cleaningWatch.Start();
	static_cast<JoinedClean<joined_pol_clean::MultiImageSet>&>(*_cleanAlgorithms[0]).ExecuteMajorIteration(residualSet, modelSet, psfImages, _imgWidth, _imgHeight, reachedMajorThreshold);
	_cleaningWatch.Pause();
	
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		_imageAllocator.Free(psfImages[ch]);
		if(hasStokesPols)
		{
			modelSet.StoreStokes(_modelImages, ch);
			residualSet.StoreStokes(_residualImages, ch);
		}
		else {
			modelSet.StoreLinear(_modelImages, ch);
			residualSet.StoreLinear(_residualImages, ch);
		}
	}
	
	updateCleanParameters(_fitsWriter, _cleanAlgorithms[0]->IterationNumber(), majorIterationNr);
	
	const PolarizationEnum
		linPols[4] = { Polarization::XX, Polarization::XY, Polarization::XY, Polarization::YY },
		stokesPols[4] = { Polarization::StokesI, Polarization::StokesQ, Polarization::StokesU, Polarization::StokesV },
		*pols = hasStokesPols ? stokesPols : linPols;
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		for(size_t i=0; i!=4; ++i)
		{
			PolarizationEnum polarization = pols[i];
			bool isImaginary = (i==2) && !hasStokesPols;
			
			if(majorIterationNr == 1)
			{
				if(_mGain == 1.0)
				{
					std::cout << "Writing residual image... " << std::flush;
					writeFits("residual.fits", residualSet.GetImage(i,ch), polarization, ch, isImaginary);
				}
				else {
					std::cout << "Writing first iteration image... " << std::flush;
					writeFits("first-residual.fits", residualSet.GetImage(i,ch), polarization, ch, isImaginary);
				}
				std::cout << "DONE\n";
			}
			if(!reachedMajorThreshold)
			{
				std::cout << "Writing model image... " << std::flush;
				writeFits("model.fits", modelSet.GetImage(i,ch), polarization, ch, isImaginary);
				std::cout << "DONE\n";
			}
		}
	}
}

void WSClean::makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary)
{
	double lowestFreq = 0.0, highestFreq = 0.0;
	const size_t size = _imgWidth * _imgHeight;
	ao::uvector<double> mfsImage(size, 0.0), addedImage(size), weightImage(size, 0.0);
	double weightSum = 0.0;
	FitsWriter writer;
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		const std::string name(getPrefix(pol, ch, isImaginary) + '-' + suffix);
		FitsReader reader(name);
		if(ch == 0)
		{
			writer = FitsWriter(reader);
			lowestFreq = reader.Frequency() - reader.Bandwidth()*0.5;
			highestFreq = reader.Frequency() + reader.Bandwidth()*0.5;
		}
		else {
			lowestFreq = std::min(lowestFreq, reader.Frequency() - reader.Bandwidth()*0.5);
			highestFreq = std::max(highestFreq, reader.Frequency() + reader.Bandwidth()*0.5);
		}
		double weight;
		if(!reader.ReadDoubleKeyIfExists("WSCIMGWG", weight))
			weight = 0.0;
		weightSum += weight;
		reader.Read(addedImage.data());
		for(size_t i=0; i!=size; ++i)
		{
			if(std::isfinite(addedImage[i])) {
				mfsImage[i] += addedImage[i] * weight;
				weightImage[i] += weight;
			}
		}
	}
	for(size_t i=0; i!=size; ++i)
		mfsImage[i] /= weightImage[i];
	writer.SetFrequency((lowestFreq+highestFreq)*0.5, highestFreq-lowestFreq);
	writer.SetExtraKeyword("WSCIMGWG", weightSum);
	writer.Write(getMFSPrefix(pol, isImaginary) + '-' + suffix, mfsImage.data());
}

void WSClean::writeFits(const string& suffix, const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary)
{
	MSSelection selection(_globalSelection);
	selectChannels(selection, channelIndex, _channelsOut);
	BandData band(_firstMSBand, selection.ChannelRangeStart(), selection.ChannelRangeEnd());
	const std::string name(getPrefix(pol, channelIndex, isImaginary) + '-' + suffix);
	initFitsWriter(_fitsWriter);
	_fitsWriter.SetPolarization(pol);
	_fitsWriter.SetFrequency(band.CentreFrequency(), band.Bandwidth());
	_fitsWriter.SetExtraKeyword("WSCIMGWG", _weightPerChannel[channelIndex]);
	size_t polIndex;
	if(_joinedPolarizationCleaning)
		polIndex = 0;
	else
		Polarization::TypeToIndex(pol, _polarizations, polIndex);
	if(_cleanAlgorithms.size()>polIndex && _cleanAlgorithms[polIndex] != 0)
	{
		setCleanParameters(_fitsWriter, *_cleanAlgorithms[polIndex]);
		updateCleanParameters(_fitsWriter, _cleanAlgorithms[polIndex]->IterationNumber(), _majorIterationNr);
	}
	_fitsWriter.Write(name, image);
}
