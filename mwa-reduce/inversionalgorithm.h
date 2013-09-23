#ifndef INVERSION_ALGORITHM_H
#define INVERSION_ALGORITHM_H

#include <cmath>
#include <string>
#include <vector>

class InversionAlgorithm
{
	public:
		enum PolarizationEnum { XX, XY, YX, YY, StokesI };
		enum WeightingEnum { NaturalWeighted, DistanceWeighted, UniformishWeighted };
		
		InversionAlgorithm() :
			_imageWidth(1024),
			_imageHeight(1024),
			_pixelSizeX((1.0 / 60.0) * M_PI / 180.0),
			_pixelSizeY((1.0 / 60.0) * M_PI / 180.0),
			_wGridSize(1),
			_intervalStart(0),
			_intervalStop(0),
			_measurementSetPaths(),
			_dataColumnName("DATA"),
			_doImagePSF(false),
			_doSubtractModel(false),
			_addToModel(false),
			_polarization(StokesI),
			_weighting(DistanceWeighted)
		{
		}
		virtual ~InversionAlgorithm()
		{
		}
		
		size_t ImageWidth() const { return _imageWidth; }
		size_t ImageHeight() const { return _imageHeight; }
		double PixelSizeX() const { return _pixelSizeX; }
		double PixelSizeY() const { return _pixelSizeY; }
		size_t WGridSize() const { return _wGridSize; }
		const std::string &MeasurementSetPath(size_t index) const { return _measurementSetPaths[index]; }
		size_t MeasurementSetCount() const { return _measurementSetPaths.size(); }
		const std::string &DataColumnName() const { return _dataColumnName; }
		bool DoImagePSF() const { return _doImagePSF; }
		bool DoSubtractModel() const { return _doSubtractModel; }
		bool AddToModel() const { return _addToModel; }
		PolarizationEnum Polarization() const { return _polarization; }
		WeightingEnum Weighting() const { return _weighting; }
		bool HasInterval() const { return _intervalStop != 0; }
		size_t IntervalStart() const { return _intervalStart; }
		size_t IntervalStop() const { return _intervalStop; }
		
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
		void SetWeighting(WeightingEnum weighting)
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
		void SetInterval(size_t intervalStart, size_t intervalStop)
		{
			_intervalStart = intervalStart;
			_intervalStop = intervalStop;
		}
		
		virtual void Invert() = 0;
		
		virtual void InvertToVisibilities(const double *image) = 0;
		
		virtual const double *ImageResult() const = 0;
		virtual double ImageResultRA() const = 0;
		virtual double ImageResultDec() const = 0;
		virtual double ImageHighestFrequencyChannel() const = 0;
		virtual double ImageLowestFrequencyChannel() const = 0;
		virtual double ImageBandStart() const = 0;
		virtual double ImageBandEnd() const = 0;
		virtual double ImageBeamSize() const = 0;
		virtual double ImageStartTime() const = 0;
		
		virtual bool HasGriddingCorrectionImage() const = 0;
		virtual void GetGriddingCorrectionImage(double *image) const = 0;
	protected:
		size_t _imageWidth, _imageHeight;
		double _pixelSizeX, _pixelSizeY;
		size_t _wGridSize, _intervalStart, _intervalStop;
		std::vector<std::string> _measurementSetPaths;
		std::string _dataColumnName;
		bool _doImagePSF, _doSubtractModel, _addToModel;
		enum PolarizationEnum _polarization;
		enum WeightingEnum _weighting;
};

#endif
