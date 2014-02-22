
#include <iostream>

#include "modelrenderer.h"
#include "model.h"
#include "imagecoordinates.h"

template<typename T>
T ModelRenderer::gaus(T x, T sigma) const
{
	long double xi = x / sigma;
	return exp(T(-0.5) * xi * xi);// / (sigma * sqrt(T(2.0) * M_PIl));
}

template void ModelRenderer::Restore(double* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double beamSize, long double startFrequency, long double endFrequency, PolarizationEnum polarization);

template<typename NumType>
void ModelRenderer::Restore(NumType* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double beamSize, long double startFrequency, long double endFrequency, PolarizationEnum polarization)
{
	// Using the FWHM formula for a Gaussian:
	long double sigma = beamSize / (2.0L * sqrtl(2.0L * logl(2.0L)));
	
	int boundingBoxSize = ceil(sigma * 20.0 / std::min(_pixelScaleL, _pixelScaleM));
	for(Model::const_iterator src=model.begin(); src!=model.end(); ++src)
	{
		for(ModelSource::const_iterator comp=src->begin(); comp!=src->end(); ++comp)
		{
			long double
				posRA = comp->PosRA(),
				posDec = comp->PosDec(),
				sourceL, sourceM;
			ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec, sourceL, sourceM);
			const SpectralEnergyDistribution &sed = comp->SED();
			const long double intFlux = sed.IntegratedFlux(startFrequency, endFrequency, polarization);
			
			//std::cout << "Source: " << comp->PosRA() << "," << comp->PosDec() << " Phase centre: " << _phaseCentreRA << "," << _phaseCentreDec << " beamsize: " << beamSize << "\n";
				
			int sourceX, sourceY;
			ImageCoordinates::LMToXY<long double>(sourceL-_phaseCentreDL, sourceM-_phaseCentreDM, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);
			//std::cout << "Adding source " << comp->Name() << " at " << sourceX << "," << sourceY << " of "
			//	<< intFlux << " Jy ("
			//	<< startFrequency/1000000.0 << "-" << endFrequency/1000000.0 << " MHz).\n";
			int
				xLeft = sourceX - boundingBoxSize,
				xRight = sourceX + boundingBoxSize,
				yTop = sourceY - boundingBoxSize,
				yBottom = sourceY + boundingBoxSize;
			if(xLeft < 0) xLeft = 0;
			if(xLeft > (int) imageWidth) xLeft = (int) imageWidth;
			if(xRight < 0) xRight = 0;
			if(xRight > (int) imageWidth) xRight = (int) imageWidth;
			if(yTop < 0) yTop = 0;
			if(yTop > (int) imageHeight) yTop = (int) imageHeight;
			if(yBottom < 0) yBottom = 0;
			if(yBottom > (int) imageHeight) yBottom = (int) imageHeight;
			
			for(int y=yTop; y!=yBottom; ++y)
			{
				NumType *imageDataPtr = imageData + y*imageWidth+xLeft;
				for(int x=xLeft; x!=xRight; ++x)
				{
					long double l, m;
					ImageCoordinates::XYToLM<long double>(x, y, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, l, m);
					l += _phaseCentreDL; m += _phaseCentreDM;
					long double dist = sqrt((l-sourceL)*(l-sourceL) + (m-sourceM)*(m-sourceM));
					long double g = gaus(dist, sigma);
					(*imageDataPtr) += NumType(g * intFlux);
					++imageDataPtr;
				}
			}
		}
	}
}

template void ModelRenderer::RenderModel(double* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double startFrequency, long double endFrequency, PolarizationEnum polarization);

template<typename NumType>
void ModelRenderer::RenderModel(NumType* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double startFrequency, long double endFrequency, PolarizationEnum polarization)
{
	for(Model::const_iterator src=model.begin(); src!=model.end(); ++src)
	{
		std::cout << "Rendering " << src->Name() << '\n';
		for(ModelSource::const_iterator comp=src->begin(); comp!=src->end(); ++comp)
		{
			long double
				posRA = comp->PosRA(),
				posDec = comp->PosDec(),
				sourceL, sourceM;
			int sourceX, sourceY;
			ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec, sourceL, sourceM);
			ImageCoordinates::LMToXY<long double>(sourceL, sourceM, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);
			
			const long double intFlux = comp->SED().IntegratedFlux(startFrequency, endFrequency, polarization);
			
			if(sourceX >= 0 && sourceX < (int) imageWidth && sourceY >= 0 && sourceY < (int) imageHeight)
			{
				NumType *imageDataPtr = imageData + sourceY*imageWidth + sourceX;
				(*imageDataPtr) += NumType(intFlux);
			}
		}
	}
}

