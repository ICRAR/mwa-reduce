#ifndef INVERSION_ALGORITHM_H
#define INVERSION_ALGORITHM_H

#include <string>
#include <cmath>

class InversionAlgorithm
{
	public:
		enum PolarizationEnum { XX, XY, YX, YY, StokesI };
		
		InversionAlgorithm() :
			_imageWidth(1024),
			_imageHeight(1024),
			_pixelSizeX((1.0 / 60.0) * M_PI / 180.0),
			_pixelSizeY((1.0 / 60.0) * M_PI / 180.0),
			_wGridSize(1),
			_measurementSetPath(),
			_dataColumnName("DATA"),
			_doImagePSF(false),
			_polarization(StokesI)
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
		const std::string &MeasurementSetPath() const { return _measurementSetPath; }
		const std::string &DataColumnName() const { return _dataColumnName; }
		bool DoImagePSF() const { return _doImagePSF; }
		PolarizationEnum Polarization() const { return _polarization; }
		
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
		void SetMeasurementSetPath(const std::string &measurementSetPath)
		{
			_measurementSetPath = measurementSetPath;
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
		
		virtual void Execute() = 0;
		
		virtual const double *ImageResult() const = 0;
		virtual double ImageResultRA() const = 0;
		virtual double ImageResultDec() const = 0;
		virtual double ImageFrequencyHigh() const = 0;
		virtual double ImageFrequencyLow() const = 0;
		virtual double ImageBeamSize() const = 0;
	protected:
		size_t _imageWidth, _imageHeight;
		double _pixelSizeX, _pixelSizeY;
		size_t _wGridSize;
		std::string _measurementSetPath, _dataColumnName;
		bool _doImagePSF;
		enum PolarizationEnum _polarization;
};

#endif
