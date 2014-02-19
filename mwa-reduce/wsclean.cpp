#include "wsclean.h"

#include "areaset.h"
#include "beamevaluator.h"
#include "cleanalgorithm.h"
#include "fitswriter.h"
#include "imagebufferallocator.h"
#include "imageweights.h"
#include "inversionalgorithm.h"
#include "modelrenderer.h"
#include "model.h"
#include "msselection.h"
#include "stopwatch.h"
#include "wsinversion.h"

#include "msprovider/contiguousms.h"
#include "msprovider/msprovider.h"
#include "msprovider/partitionedms.h"

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
	_commandLine()
{
}

void initFitsWriter(FitsWriter& writer, const InversionAlgorithm& inversionAlgorithm, double beamSizeArcmin)
{
	double
		ra = inversionAlgorithm.PhaseCentreRA(),
		dec = inversionAlgorithm.PhaseCentreDec(),
		pixelScaleX = inversionAlgorithm.PixelSizeX(),
		pixelScaleY = inversionAlgorithm.PixelSizeY(),
		freqHigh = inversionAlgorithm.HighestFrequencyChannel(),
		freqLow = inversionAlgorithm.LowestFrequencyChannel(),
		freqCentre = (freqHigh + freqLow) * 0.5,
		bandwidth = inversionAlgorithm.BandEnd() - inversionAlgorithm.BandStart(),
		beamSize = inversionAlgorithm.BeamSize(),
		dateObs = inversionAlgorithm.StartTime();
		
	writer.SetImageDimensions(inversionAlgorithm.ImageWidth(), inversionAlgorithm.ImageHeight(), ra, dec, pixelScaleX, pixelScaleY);
	writer.SetFrequency(freqCentre, bandwidth);
	writer.SetDate(dateObs);
	writer.SetPolarization(inversionAlgorithm.Polarization());
	writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
	writer.AddHistory(commandLine);
	if(beamSizeArcmin != 0.0) {
		double beamSizeRad = beamSizeArcmin * (M_PI / 60.0 / 180.0);
		writer.SetBeamInfo(beamSizeRad, beamSizeRad, 0.0);
	}
	else {
		writer.SetBeamInfo(beamSize);
	}
	if(inversionAlgorithm.HasDenormalPhaseCentre())
		writer.SetPhaseCentreShift(inversionAlgorithm.PhaseCentreDL(), inversionAlgorithm.PhaseCentreDM());
	
	writer.SetExtraKeyword("WSCNWLAY", inversionAlgorithm.WGridSize());
	writer.SetExtraKeyword("WSCDATAC", inversionAlgorithm.DataColumnName());
	writer.SetExtraKeyword("WSCWEIGH", inversionAlgorithm.Weighting().ToString());
	writer.SetExtraKeyword("WSCGKRNL", inversionAlgorithm.AntialiasingKernelSize());
	if(inversionAlgorithm.Selection().HasChannelRange())
	{
		writer.SetExtraKeyword("WSCCHANS", inversionAlgorithm.Selection().ChannelRangeStart());
		writer.SetExtraKeyword("WSCCHANE", inversionAlgorithm.Selection().ChannelRangeEnd());
	}
	if(inversionAlgorithm.Selection().HasInterval())
	{
		writer.SetExtraKeyword("WSCTIMES", inversionAlgorithm.Selection().IntervalStart());
		writer.SetExtraKeyword("WSCTIMEE", inversionAlgorithm.Selection().IntervalEnd());
	}
	writer.SetExtraKeyword("WSCFIELD", inversionAlgorithm.Selection().FieldId());
}

void SetCleanParameters(FitsWriter& writer, const CleanAlgorithm& clean)
{
	writer.SetExtraKeyword("WSCNITER", clean.MaxNIter());
	writer.SetExtraKeyword("WSCTHRES", clean.Threshold());
	writer.SetExtraKeyword("WSCGAIN", clean.SubtractionGain());
	writer.SetExtraKeyword("WSCMGAIN", clean.StopGain());
	writer.SetExtraKeyword("WSCNEGCM", clean.AllowNegativeComponents());
	writer.SetExtraKeyword("WSCNEGST", clean.StopOnNegativeComponents());
	writer.SetExtraKeyword("WSCSMPSF", clean.ResizePSF());
}

void UpdateCleanParameters(FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr)
{
	writer.SetExtraKeyword("WSCMINOR", minorIterationNr);
	writer.SetExtraKeyword("WSCMAJOR", majorIterationNr);
}

void WSClean::Run()
{
	ImageBufferAllocator<double> imageAllocator;
	
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

	bool doReorder = ((_channelsOut != 1) || _forceReorder) && !_forceNoReorder;
	
	std::vector<PartitionedMS::Handle> partitionedMSHandles;
	if(doReorder)
	{
		for(std::vector<std::string>::const_iterator i=_filenames.begin(); i != _filenames.end(); ++i)
		{
			partitionedMSHandles.push_back(PartitionedMS::Partition(*i, _channelsOut, _globalSelection, _columnName, true, _mGain != 1.0, _polarizations));
		}
	}
	
	for(size_t outChannelIndex=0; outChannelIndex!=_channelsOut; ++outChannelIndex)
	{
		MSSelection partSelection = _globalSelection;
		std::unique_ptr<InversionAlgorithm> inversionAlgorithm(new WSInversion(&imageAllocator, _memFraction));
		static_cast<WSInversion&>(*inversionAlgorithm).SetGridMode(_gridMode);
		
		for(std::vector<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
		{
			std::vector<MSProvider*> msProviders;
			for(size_t i=0; i != _filenames.size(); ++i) {
				MSProvider* msProvider;
				if(doReorder)
					msProvider = new PartitionedMS(partitionedMSHandles[i], outChannelIndex, *curPol);
				else
					msProvider = new ContiguousMS(_filenames[i], _columnName, partSelection, *curPol, _mGain != 1.0);
				inversionAlgorithm->AddMeasurementSet(msProvider);
				msProviders.push_back(msProvider);
			}
			if(doReorder)
			{
				size_t startCh, endCh;
				if(_globalSelection.HasChannelRange())
				{
					startCh = _globalSelection.ChannelRangeStart();
					endCh = _globalSelection.ChannelRangeEnd();
				}
				else {
					BandData band(inversionAlgorithm->MeasurementSet(0).MS().spectralWindow());
					startCh = 0;
					endCh = band.ChannelCount();
				}
				size_t newStart = startCh + (endCh - startCh) * outChannelIndex / _channelsOut;
				size_t newEnd = startCh + (endCh - startCh) * (outChannelIndex+1) / _channelsOut;
				partSelection.SetChannelRange(newStart, newEnd);
			}
			std::string partPrefixName;
			if(_channelsOut != 1)
			{
				std::ostringstream partPrefixNameStr;
				partPrefixNameStr << _prefixName << '-';
				if(outChannelIndex < 1000) partPrefixNameStr << '0';
				if(outChannelIndex < 100) partPrefixNameStr << '0';
				if(outChannelIndex < 10) partPrefixNameStr << '0';
				partPrefixNameStr << outChannelIndex;
				partPrefixName = partPrefixNameStr.str();
			}
			else {
				partPrefixName = _prefixName;
			}
			
			// Initialize weight grid if necessary.
			std::unique_ptr<ImageWeights> imageWeights;
			if(_weightMode.RequiresGridding())
			{
				std::cout << "Precalculating weights for " << _weightMode.ToString() << " weighting... " << std::flush;
				imageWeights.reset(new ImageWeights(_imgWidth, _imgHeight, _pixelScaleX, _pixelScaleY, _weightMode.SuperWeight()));
				for(size_t i=0; i!=inversionAlgorithm->MeasurementSetCount(); ++i)
				{
					imageWeights->Grid(inversionAlgorithm->MeasurementSet(i), _weightMode, partSelection);
					if(inversionAlgorithm->MeasurementSetCount() > 1)
						std::cout << i << ' ' << std::flush;
				}
				std::cout << "DONE\n";
				inversionAlgorithm->SetPrecalculatedWeightInfo(imageWeights.get());
			}
			
			inversionAlgorithm->SetImageWidth(_imgWidth);
			inversionAlgorithm->SetImageHeight(_imgHeight);
			inversionAlgorithm->SetPixelSizeX(_pixelScaleX);
			inversionAlgorithm->SetPixelSizeY(_pixelScaleY);
			if(_nWLayers != 0)
				inversionAlgorithm->SetWGridSize(_nWLayers);
			else
				inversionAlgorithm->SetNoWGridSize();
			inversionAlgorithm->SetAntialiasingKernelSize(_antialiasingKernelSize);
			inversionAlgorithm->SetOverSamplingFactor(_overSamplingFactor);
			inversionAlgorithm->SetPolarization(*curPol);
			inversionAlgorithm->SetDataColumnName(_columnName);
			inversionAlgorithm->SetWeighting(_weightMode);
			inversionAlgorithm->SetImaginaryPart(_imaginaryPart);
			inversionAlgorithm->SetSelection(partSelection);
			inversionAlgorithm->SetWLimit(_wLimit/100.0);
			
			double* psf = 0;
			bool isFirstInversion = true;
			
			Stopwatch inversionWatch(false), predictingWatch(false), cleaningWatch(false);
			
			if(_nIter > 0 || _makePSF)
			{
				std::cout << std::flush << " == Constructing PSF ==\n";
				inversionWatch.Start();
				inversionAlgorithm->SetDoImagePSF(true);
				inversionAlgorithm->SetVerbose(isFirstInversion);
				inversionAlgorithm->Invert();
					
				psf = imageAllocator.Allocate(_imgWidth * _imgHeight);
				memcpy(psf, inversionAlgorithm->ImageResult(), _imgWidth * _imgHeight * sizeof(double));
				inversionWatch.Pause();
				
				isFirstInversion = false;
				std::cout << "Beam size is " << inversionAlgorithm->BeamSize()*(180.0*60.0/M_PI) << " arcmin.\n";
				
				std::cout << "Writing psf image... " << std::flush;
				FitsWriter fitsWriter;
				initFitsWriter(fitsWriter, *inversionAlgorithm, _beamSize);
				fitsWriter.Write(std::string(partPrefixName) + "-psf.fits", psf);
				std::cout << "DONE\n";
				
				CleanAlgorithm::RemoveNaNsInPSF(psf, _imgWidth, _imgHeight);
			}
			
			std::cout << std::flush << " == Constructing image ==\n";
			inversionWatch.Start();
			if(_nWLayers != 0)
				inversionAlgorithm->SetWGridSize(_nWLayers);
			else
				inversionAlgorithm->SetNoWGridSize();
			inversionAlgorithm->SetDoImagePSF(false);
			inversionAlgorithm->SetVerbose(isFirstInversion);
			inversionAlgorithm->Invert();
			inversionWatch.Pause();
			inversionAlgorithm->SetVerbose(false);
			
			if(inversionAlgorithm->HasGriddingCorrectionImage())
			{
				std::cout << "Writing gridding correction image... " << std::flush;
				double* gridding = imageAllocator.Allocate(_imgWidth * _imgHeight);
				inversionAlgorithm->GetGriddingCorrectionImage(&gridding[0]);
				FitsWriter fitsWriter;
				initFitsWriter(fitsWriter, *inversionAlgorithm, _beamSize);
				fitsWriter.Write(std::string(partPrefixName) + "-gridding.fits", &gridding[0]);
				imageAllocator.Free(gridding);
				std::cout << "DONE\n";
			}
				
			double
				*modelImage = imageAllocator.Allocate(_imgWidth * _imgHeight),
				*residual = imageAllocator.Allocate(_imgWidth * _imgHeight);
			memcpy(residual, inversionAlgorithm->ImageResult(), _imgWidth * _imgHeight * sizeof(double));
			memset(modelImage, 0, _imgWidth * _imgHeight * sizeof(double));
			
			std::cout << "Writing dirty image... " << std::flush;
			FitsWriter fitsWriter;
			initFitsWriter(fitsWriter, *inversionAlgorithm, _beamSize);
			fitsWriter.Write(std::string(partPrefixName) + "-dirty.fits", residual);
			std::cout << "DONE\n";
			
			if(_mGain == 1.0)
				inversionAlgorithm.reset();
			
			CleanAlgorithm cleanAlgorithm;
			cleanAlgorithm.SetMaxNIter(_nIter);
			cleanAlgorithm.SetThreshold(_threshold);
			cleanAlgorithm.SetSubtractionGain(_gain);
			cleanAlgorithm.SetStopGain(_mGain);
			cleanAlgorithm.SetAllowNegativeComponents(_allowNegative);
			cleanAlgorithm.SetStopOnNegativeComponents(_stopOnNegative);
			cleanAlgorithm.SetResizePSF(_smallPSF);

			SetCleanParameters(fitsWriter, cleanAlgorithm);
			UpdateCleanParameters(fitsWriter, 0, 0);
				
			if(_nIter > 0)
			{
				std::unique_ptr<AreaSet> cleanAreas;
				if(!_cleanAreasFilename.empty())
				{
					cleanAreas.reset(new AreaSet());
					AreaParser parser;
					std::ifstream caFile(_cleanAreasFilename.c_str());
					parser.Parse(*cleanAreas, caFile);
					cleanAreas->SetImageProperties(_pixelScaleX, _pixelScaleY, inversionAlgorithm->PhaseCentreRA(), inversionAlgorithm->PhaseCentreDec(), _imgWidth, _imgHeight);
					cleanAlgorithm.SetCleanAreas(*cleanAreas);
				}
				
				// Start major cleaning loop
				size_t majorIterationNr = 1;
				bool reachedMajorThreshold = false;
				do {
					std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
					cleaningWatch.Start();
					cleanAlgorithm.ExecuteMajorIteration(residual, modelImage, psf, _imgWidth, _imgHeight, reachedMajorThreshold);
					cleaningWatch.Pause();
					
					UpdateCleanParameters(fitsWriter, cleanAlgorithm.IterationNumber(), majorIterationNr);
					
					if(majorIterationNr == 1)
					{
						if(_mGain == 1.0)
						{
							std::cout << "Writing residual image... " << std::flush;
							fitsWriter.Write(std::string(partPrefixName) + "-residual.fits", &residual[0]);
						}
						else {
							std::cout << "Writing first iteration image... " << std::flush;
							fitsWriter.Write(std::string(partPrefixName) + "-first-residual.fits", &residual[0]);
						}
						std::cout << "DONE\n";
					}
					
					if(!reachedMajorThreshold)
					{
						std::cout << "Writing model image... " << std::flush;
						fitsWriter.Write(std::string(partPrefixName) + "-model.fits", &modelImage[0]);
						std::cout << "DONE\n";
					}
					
					if(_mGain != 1.0)
					{
						std::cout << std::flush << " == Converting model image to visibilities ==\n";
						predictingWatch.Start();
						inversionAlgorithm->SetAddToModel(false);
						inversionAlgorithm->InvertToVisibilities(&modelImage[0]);
						predictingWatch.Pause();
						
						std::cout << std::flush << " == Constructing image ==\n";
						inversionWatch.Start();
						inversionAlgorithm->SetDoSubtractModel(true);
						inversionAlgorithm->Invert();
						inversionWatch.Pause();
						
						memcpy(&residual[0], inversionAlgorithm->ImageResult(), _imgWidth * _imgHeight * sizeof(double));
						
						if(!reachedMajorThreshold)
						{
							// This was the final major iteration: clean up & save results
							inversionAlgorithm.reset();
							
							std::cout << "Writing residual image... " << std::flush;
							fitsWriter.Write(std::string(partPrefixName) + "-residual.fits", &residual[0]);
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
			CleanAlgorithm::GetModelFromImage(model, &modelImage[0], _imgWidth, _imgHeight, fitsWriter.RA(), fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, 0.0, fitsWriter.Frequency(), *curPol);
			if(!_saveModelFilename.empty())
			{
				std::cout << "Saving model to " << _saveModelFilename << "... " << std::flush;
				model.Save(_saveModelFilename.c_str());
			}
			
			std::cout << "Rendering " << model.SourceCount() << " sources to restored image... " << std::flush;
			ModelRenderer renderer(fitsWriter.RA(), fitsWriter.Dec(), _pixelScaleX, _pixelScaleY);
			size_t polarizationIndex;
			switch(*curPol)
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
			fitsWriter.Write(std::string(partPrefixName) + "-image.fits", &residual[0]);
			std::cout << "DONE\n";
			
			imageAllocator.ReportStatistics();
			std::cout << "Inversion: " << inversionWatch.ToString() << ", prediction: " << predictingWatch.ToString() << ", cleaning: " << cleaningWatch.ToString() << '\n';
			
			imageAllocator.Free(residual);
			imageAllocator.Free(modelImage);
			imageAllocator.Free(psf);
			
			for(std::vector<MSProvider*>::iterator i=msProviders.begin(); i != msProviders.end(); ++i)
				delete *i;
		}
	}
}
