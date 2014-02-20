#include "wsclean.h"

#include "areaset.h"
#include "beamevaluator.h"
#include "cleanalgorithm.h"
#include "fitswriter.h"
#include "imageweights.h"
#include "inversionalgorithm.h"
#include "modelrenderer.h"
#include "model.h"
#include "msselection.h"
#include "wsinversion.h"

#include "msprovider/contiguousms.h"

#include "parser/areaparser.h"

#include <iostream>
#include <memory>

std::string commandLine;

WSClean::WSClean() :
	_imgWidth(2048), _imgHeight(2048), _channelsOut(1),
	_pixelScaleX(0.01 * M_PI / 180.0), _pixelScaleY(0.01 * M_PI / 180.0),
	_threshold(0.0), _gain(0.1), _mGain(1.0), _beamSize(0.0), _memFraction(1.0), _wLimit(0.0),
	_nWLayers(0), _nIter(0), _antialiasingKernelSize(7), _overSamplingFactor(63),
	_globalSelection(),
	_columnName(), _addModelFilename(), _saveModelFilename(), _cleanAreasFilename(),
	_polarizations(1, Polarization::StokesI),
	_weightMode(WeightMode::UniformWeighted),
	_prefixName("wsclean"),
	_allowNegative(true), _smallPSF(false), _addApparentModel(false), _stopOnNegative(false), _imaginaryPart(false), _makePSF(false),
	_forceReorder(false), _forceNoReorder(false),
	_gridMode(LayeredImager::KaiserBessel),
	_filenames(),
	_commandLine(),
	_inversionWatch(false), _predictingWatch(false), _cleaningWatch(false),
	_isFirstInversion(true)
{
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
	if(_beamSize != 0.0) {
		double beamSizeRad = _beamSize * (M_PI / 60.0 / 180.0);
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
		
	_psfImage = _imageAllocator.Allocate(_imgWidth * _imgHeight);
	memcpy(_psfImage, _inversionAlgorithm->ImageResult(), _imgWidth * _imgHeight * sizeof(double));
	_inversionWatch.Pause();
	
	_isFirstInversion = false;
	std::cout << "Beam size is " << _inversionAlgorithm->BeamSize()*(180.0*60.0/M_PI) << " arcmin.\n";
	
	std::cout << "Writing psf image... " << std::flush;
	FitsWriter fitsWriter;
	initFitsWriter(fitsWriter);
	fitsWriter.Write(_prefixName + "-psf.fits", _psfImage);
	std::cout << "DONE\n";
	
	CleanAlgorithm::RemoveNaNsInPSF(_psfImage, _imgWidth, _imgHeight);
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

void WSClean::imageMainFirst()
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
	_inversionAlgorithm->SetDataColumnName(_columnName);
	_inversionAlgorithm->SetWeighting(_weightMode);
	_inversionAlgorithm->SetImaginaryPart(_imaginaryPart);
	_inversionAlgorithm->SetSelection(_currentPartSelection);
	_inversionAlgorithm->SetWLimit(_wLimit/100.0);
}

void WSClean::imageMainNonFirst()
{
	std::cout << std::flush << " == Constructing image ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoSubtractModel(true);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
}

void WSClean::predict(const double* modelImage)
{
	std::cout << std::flush << " == Converting model image to visibilities ==\n";
	_predictingWatch.Start();
	_inversionAlgorithm->SetAddToModel(false);
	_inversionAlgorithm->InvertToVisibilities(modelImage);
	_predictingWatch.Pause();
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

	_doReorder = ((_channelsOut != 1) || _forceReorder) && !_forceNoReorder;
	
	std::vector<PartitionedMS::Handle> partitionedMSHandles;
	if(_doReorder)
	{
		for(std::vector<std::string>::const_iterator i=_filenames.begin(); i != _filenames.end(); ++i)
		{
			partitionedMSHandles.push_back(PartitionedMS::Partition(*i, _channelsOut, _globalSelection, _columnName, true, _mGain != 1.0, _polarizations));
		}
	}
	
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
			BandData band(_inversionAlgorithm->MeasurementSet(0).MS().spectralWindow());
			startCh = 0;
			endCh = band.ChannelCount();
		}
		size_t newStart = startCh + (endCh - startCh) * outChannelIndex / _channelsOut;
		size_t newEnd = startCh + (endCh - startCh) * (outChannelIndex+1) / _channelsOut;
		partSelection.SetChannelRange(newStart, newEnd);
	}
		
	std::string rootPrefix = _prefixName;
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
		
	for(std::vector<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
	{
		runPolarization(outChannelIndex, *curPol);
	}
	
	_prefixName = rootPrefix;
}

void WSClean::runPolarization(size_t outChannelIndex, PolarizationEnum polarization)
{
	std::vector<MSProvider*> msProviders;
	for(size_t i=0; i != _filenames.size(); ++i)
	{
		MSProvider* msProvider;
		if(_doReorder)
			msProvider = new PartitionedMS(_partitionedMSHandles[i], outChannelIndex, polarization);
		else
			msProvider = new ContiguousMS(_filenames[i], _columnName, _currentPartSelection, polarization, _mGain != 1.0);
		_inversionAlgorithm->AddMeasurementSet(msProvider);
		msProviders.push_back(msProvider);
	}
	
	initializeImageWeights(_currentPartSelection);
	
	prepareInversionAlgorithm(polarization);
	
	const bool firstBeforePSF = _isFirstInversion;
	_psfImage = 0;
	if(_nIter > 0 || _makePSF)
		imagePSF();
	
	imageMainFirst();
	
	if(firstBeforePSF && _inversionAlgorithm->HasGriddingCorrectionImage())
		imageGridding();
	
	_isFirstInversion = false;
		
	double
		*modelImage = _imageAllocator.Allocate(_imgWidth * _imgHeight),
		*residual = _imageAllocator.Allocate(_imgWidth * _imgHeight);
	memcpy(residual, _inversionAlgorithm->ImageResult(), _imgWidth * _imgHeight * sizeof(double));
	memset(modelImage, 0, _imgWidth * _imgHeight * sizeof(double));
	
	std::cout << "Writing dirty image... " << std::flush;
	FitsWriter fitsWriter;
	initFitsWriter(fitsWriter);
	fitsWriter.Write(std::string(_prefixName) + "-dirty.fits", residual);
	std::cout << "DONE\n";
	
	if(_mGain == 1.0)
		_inversionAlgorithm.reset();
	
	CleanAlgorithm cleanAlgorithm;
	cleanAlgorithm.SetMaxNIter(_nIter);
	cleanAlgorithm.SetThreshold(_threshold);
	cleanAlgorithm.SetSubtractionGain(_gain);
	cleanAlgorithm.SetStopGain(_mGain);
	cleanAlgorithm.SetAllowNegativeComponents(_allowNegative);
	cleanAlgorithm.SetStopOnNegativeComponents(_stopOnNegative);
	cleanAlgorithm.SetResizePSF(_smallPSF);

	setCleanParameters(fitsWriter, cleanAlgorithm);
	updateCleanParameters(fitsWriter, 0, 0);
		
	if(_nIter > 0)
	{
		std::unique_ptr<AreaSet> cleanAreas;
		if(!_cleanAreasFilename.empty())
		{
			cleanAreas.reset(new AreaSet());
			AreaParser parser;
			std::ifstream caFile(_cleanAreasFilename.c_str());
			parser.Parse(*cleanAreas, caFile);
			cleanAreas->SetImageProperties(_pixelScaleX, _pixelScaleY, _inversionAlgorithm->PhaseCentreRA(), _inversionAlgorithm->PhaseCentreDec(), _imgWidth, _imgHeight);
			cleanAlgorithm.SetCleanAreas(*cleanAreas);
		}
		
		// Start major cleaning loop
		size_t majorIterationNr = 1;
		bool reachedMajorThreshold = false;
		do {
			std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
			_cleaningWatch.Start();
			cleanAlgorithm.ExecuteMajorIteration(residual, modelImage, _psfImage, _imgWidth, _imgHeight, reachedMajorThreshold);
			_cleaningWatch.Pause();
			
			updateCleanParameters(fitsWriter, cleanAlgorithm.IterationNumber(), majorIterationNr);
			
			if(majorIterationNr == 1)
			{
				if(_mGain == 1.0)
				{
					std::cout << "Writing residual image... " << std::flush;
					fitsWriter.Write(_prefixName + "-residual.fits", &residual[0]);
				}
				else {
					std::cout << "Writing first iteration image... " << std::flush;
					fitsWriter.Write(_prefixName + "-first-residual.fits", &residual[0]);
				}
				std::cout << "DONE\n";
			}
			
			if(!reachedMajorThreshold)
			{
				std::cout << "Writing model image... " << std::flush;
				fitsWriter.Write(_prefixName + "-model.fits", &modelImage[0]);
				std::cout << "DONE\n";
			}
			
			if(_mGain != 1.0)
			{
				predict(&modelImage[0]);
				
				imageMainNonFirst();
				
				memcpy(&residual[0], _inversionAlgorithm->ImageResult(), _imgWidth * _imgHeight * sizeof(double));
				
				if(!reachedMajorThreshold)
				{
					// This was the final major iteration: clean up & save results
					_inversionAlgorithm.reset();
					
					std::cout << "Writing residual image... " << std::flush;
					fitsWriter.Write(_prefixName + "-residual.fits", &residual[0]);
					std::cout << "DONE\n";
				}
				
				++majorIterationNr;
			}
			
		} while(reachedMajorThreshold);
		
		std::cout << majorIterationNr << " major iterations were performed.\n";
	}
	
	Model model;
	if(!_addModelFilename.empty())
	{
		std::cout << "Reading model from " << _addModelFilename << "... " << std::flush;
		model = Model(_addModelFilename.c_str());
		if(_addApparentModel)
		{
			casa::MeasurementSet ms(_filenames.front());
			BeamEvaluator beamEval(ms, false);
			beamEval.AbsToApparent(model);
		}
	}
	CleanAlgorithm::GetModelFromImage(model, &modelImage[0], _imgWidth, _imgHeight, fitsWriter.RA(), fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, 0.0, fitsWriter.Frequency(), polarization);
	if(!_saveModelFilename.empty())
	{
		std::cout << "Saving model to " << _saveModelFilename << "... " << std::flush;
		model.Save(_saveModelFilename.c_str());
	}
	
	std::cout << "Rendering " << model.SourceCount() << " sources to restored image... " << std::flush;
	ModelRenderer renderer(fitsWriter.RA(), fitsWriter.Dec(), _pixelScaleX, _pixelScaleY);
	size_t polarizationIndex;
	switch(polarization)
	{
		case Polarization::StokesI:
		case Polarization::XX: polarizationIndex = 0; break;
		case Polarization::XY: polarizationIndex = 1; break;
		case Polarization::YX: polarizationIndex = 2; break;
		case Polarization::YY: polarizationIndex = 3; break;
		default: throw std::runtime_error("Unsupported polarization");
	}
	double
		freqLow = fitsWriter.Frequency() - fitsWriter.Bandwidth()*0.5,
		freqHigh = fitsWriter.Frequency() + fitsWriter.Bandwidth()*0.5;
	renderer.Restore(&residual[0], _imgWidth, _imgHeight, model, fitsWriter.BeamSizeMajorAxis(), freqLow, freqHigh, polarizationIndex);
	std::cout << "DONE\n";
	
	std::cout << "Writing restored image... " << std::flush;
	fitsWriter.Write(std::string(_prefixName) + "-image.fits", &residual[0]);
	std::cout << "DONE\n";
	
	_imageAllocator.ReportStatistics();
	std::cout << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", cleaning: " << _cleaningWatch.ToString() << '\n';
	
	_imageAllocator.Free(residual);
	_imageAllocator.Free(modelImage);
	_imageAllocator.Free(_psfImage);
	
	for(std::vector<MSProvider*>::iterator i=msProviders.begin(); i != msProviders.end(); ++i)
		delete *i;
}

