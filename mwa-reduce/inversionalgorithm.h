#ifndef INVERSION_ALGORITHM_H
#define INVERSION_ALGORITHM_H

#include "polarizationenum.h"
#include "msselection.h"
#include "weightmode.h"

#include <cmath>
#include <string>
#include <vector>

class InversionAlgorithm
{
	public:
		InversionAlgorithm() :
			_imageWidth(0),
			_imageHeight(0),
			_pixelSizeX((1.0 / 60.0) * M_PI / 180.0),
			_pixelSizeY((1.0 / 60.0) * M_PI / 180.0),
			_wGridSize(0),
			_measurementSets(),
			_dataColumnName("DATA"),
			_doImagePSF(false),
			_doSubtractModel(false),
			_addToModel(false),
			_smallInversion(false),
			_wLimit(0.0),
			_precalculatedWeightInfo(0),
			_polarization(Polarization::StokesI),
			_isComplex(false),
			_weighting(WeightMode::UniformWeighted),
			_verbose(false),
			_selection(),
			_antialiasingKernelSize(7),
			_overSamplingFactor(63)
		{
		}
		virtual ~InversionAlgorithm()
		{
		}
		
		size_t ImageWidth() const { return _imageWidth; }
		size_t ImageHeight() const { return _imageHeight; }
		double PixelSizeX() const { return _pixelSizeX; }
		double PixelSizeY() const { return _pixelSizeY; }
		bool HasWGridSize() const { return _wGridSize != 0; }
		size_t WGridSize() const { return _wGridSize; }
		void ClearMeasurementSetList() { _measurementSets.clear(); }
		class MSProvider &MeasurementSet(size_t index) const { return *_measurementSets[index]; }
		size_t MeasurementSetCount() const { return _measurementSets.size(); }
		const std::string &DataColumnName() const { return _dataColumnName; }
		bool DoImagePSF() const { return _doImagePSF; }
		bool DoSubtractModel() const { return _doSubtractModel; }
		bool AddToModel() const { return _addToModel; }
		bool SmallInversion() const { return _smallInversion; }
		PolarizationEnum Polarization() const { return _polarization; }
		WeightMode Weighting() const { return _weighting; }
		class ImageWeights* PrecalculatedWeightInfo() const { return _precalculatedWeightInfo; }
		const MSSelection& Selection() const { return _selection; }
		bool IsComplex() const { return _isComplex; }
		bool Verbose() const { return _verbose; }
		size_t AntialiasingKernelSize() const { return _antialiasingKernelSize; }
		size_t OverSamplingFactor() const { return _overSamplingFactor; }
		bool HasWLimit() const { return _wLimit != 0.0; }
		double WLimit() const { return _wLimit; }
		
		void SetImageWidth(size_t imageWidth)
		{
			_imageWidth = imageWidth;
		}
		void SetImageHeight(size_t imageHeight)
		{
			_imageHeight = imageHeight;
		}
		void SetPixelSizeX(double pixelSizeX)
		{
			_pixelSizeX = pixelSizeX;
		}
		void SetPixelSizeY(double pixelSizeY)
		{
			_pixelSizeY = pixelSizeY;
		}
		void SetWGridSize(size_t wGridSize)
		{
			_wGridSize = wGridSize;
		}
		void SetNoWGridSize()
		{
			_wGridSize = 0;
		}
		void AddMeasurementSet(class MSProvider* msProvider)
		{
			_measurementSets.push_back(msProvider);
		}
		void SetDataColumnName(const std::string &dataColumnName)
		{
			_dataColumnName = dataColumnName;
		}
		void SetDoImagePSF(bool doImagePSF)
		{
			_doImagePSF = doImagePSF;
		}
		void SetPolarization(PolarizationEnum polarization)
		{
			_polarization = polarization;
		}
		void SetIsComplex(bool isComplex)
		{
			_isComplex = isComplex;
		}
		void SetWeighting(WeightMode weighting)
		{
			_weighting = weighting;
		}
		void SetDoSubtractModel(bool doSubtractModel)
		{
			_doSubtractModel = doSubtractModel;
		}
		void SetAddToModel(bool addToModel)
		{
			_addToModel = addToModel;
		}
		void SetSmallInversion(bool smallInversion)
		{ 
			_smallInversion = smallInversion;
		}
		void SetPrecalculatedWeightInfo(class ImageWeights* precalculatedWeightInfo)
		{ 
			_precalculatedWeightInfo = precalculatedWeightInfo;
		}
		void SetSelection(const MSSelection& selection)
		{
			_selection = selection;
		}
		void SetVerbose(bool verbose)
		{
			_verbose = verbose;
		}
		void SetAntialiasingKernelSize(size_t kernelSize)
		{
			_antialiasingKernelSize = kernelSize;
		}
		void SetOverSamplingFactor(size_t factor)
		{
			_overSamplingFactor = factor;
		}
		void SetWLimit(double wLimit)
		{
			_wLimit = wLimit;
		}
		
		virtual void Invert() = 0;
		
		virtual void Predict(double* image) = 0;
		virtual void Predict(double* real, double* imaginary) = 0;
		
		virtual double *ImageRealResult() const = 0;
		virtual double *ImageImaginaryResult() const = 0;
		virtual double PhaseCentreRA() const = 0;
		virtual double PhaseCentreDec() const = 0;
		virtual bool HasDenormalPhaseCentre() const { return false; }
		virtual double PhaseCentreDL() const = 0;
		virtual double PhaseCentreDM() const = 0;
		virtual double HighestFrequencyChannel() const = 0;
		virtual double LowestFrequencyChannel() const = 0;
		virtual double BandStart() const = 0;
		virtual double BandEnd() const = 0;
		virtual double BeamSize() const = 0;
		virtual double StartTime() const = 0;
		virtual double ImageWeight() const = 0;
		
		/**
		 * Deallocate any data that is no longer necessary, but all methods
		 * will still return results from the imaging, with the exception of
		 * except ImageReal/ImagResult() and GetGriddingCorrectionImage().
		 */
		virtual void FreeImagingData() {	}
		
		virtual bool HasGriddingCorrectionImage() const = 0;
		virtual void GetGriddingCorrectionImage(double *image) const = 0;
		
		virtual size_t ActualInversionWidth() const { return _imageWidth; }
		virtual size_t ActualInversionHeight() const { return _imageHeight; }
	private:
		size_t _imageWidth, _imageHeight;
		double _pixelSizeX, _pixelSizeY;
		size_t _wGridSize;
		std::vector<MSProvider*> _measurementSets;
		std::string _dataColumnName;
		bool _doImagePSF, _doSubtractModel, _addToModel, _smallInversion;
		double _wLimit;
		class ImageWeights *_precalculatedWeightInfo;
		PolarizationEnum _polarization;
		bool _isComplex;
		WeightMode _weighting;
		bool _verbose;
		MSSelection _selection;
		size_t _antialiasingKernelSize, _overSamplingFactor;
};

#endif
