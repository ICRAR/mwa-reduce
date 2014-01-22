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
			_imageWidth(1024),
			_imageHeight(1024),
			_pixelSizeX((1.0 / 60.0) * M_PI / 180.0),
			_pixelSizeY((1.0 / 60.0) * M_PI / 180.0),
			_wGridSize(0),
			_measurementSetPaths(),
			_dataColumnName("DATA"),
			_doImagePSF(false),
			_doSubtractModel(false),
			_addToModel(false),
			_precalculatedWeightInfo(0),
			_polarization(Polarization::StokesI),
			_imaginaryPart(false),
			_weighting(WeightMode::UniformWeighted),
			_verbose(false),
			_selection(),
			_antialiasingKernelSize(7),
			_overSamplingFactor(15)
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
		const std::string &MeasurementSetPath(size_t index) const { return _measurementSetPaths[index]; }
		size_t MeasurementSetCount() const { return _measurementSetPaths.size(); }
		const std::string &DataColumnName() const { return _dataColumnName; }
		bool DoImagePSF() const { return _doImagePSF; }
		bool DoSubtractModel() const { return _doSubtractModel; }
		bool AddToModel() const { return _addToModel; }
		PolarizationEnum Polarization() const { return _polarization; }
		WeightMode Weighting() const { return _weighting; }
		class ImageWeights* PrecalculatedWeightInfo() const { return _precalculatedWeightInfo; }
		const MSSelection& Selection() const { return _selection; }
		bool ImaginaryPart() const { return _imaginaryPart; }
		bool Verbose() const { return _verbose; }
		size_t AntialiasingKernelSize() const { return _antialiasingKernelSize; }
		size_t OverSamplingFactor() const { return _overSamplingFactor; }
		
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
		void AddMeasurementSetPath(const std::string &measurementSetPath)
		{
			_measurementSetPaths.push_back(measurementSetPath);
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
		void SetImaginaryPart(bool imaginaryPart)
		{
			_imaginaryPart = imaginaryPart;
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
		
		virtual void Invert() = 0;
		
		virtual void InvertToVisibilities(const double *image) = 0;
		
		virtual const double *ImageResult() const = 0;
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
		
		virtual bool HasGriddingCorrectionImage() const = 0;
		virtual void GetGriddingCorrectionImage(double *image) const = 0;
	protected:
		size_t _imageWidth, _imageHeight;
		double _pixelSizeX, _pixelSizeY;
		size_t _wGridSize;
		std::vector<std::string> _measurementSetPaths;
		std::string _dataColumnName;
		bool _doImagePSF, _doSubtractModel, _addToModel;
		class ImageWeights *_precalculatedWeightInfo;
		PolarizationEnum _polarization;
		bool _imaginaryPart;
		WeightMode _weighting;
		bool _verbose;
		MSSelection _selection;
		size_t _antialiasingKernelSize, _overSamplingFactor;
};

#endif
