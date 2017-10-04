
#include <iostream>

#include "fitsreader.h"
#include "modelrenderer.h"
#include "model/model.h"
#include "units/imagecoordinates.h"
#include "uvector.h"
#include "fftconvolver.h"

template<typename T>
T ModelRenderer::gaus(T x, T sigma)
{
	long double xi = x / sigma;
	return exp(T(-0.5) * xi * xi);// / (sigma * sqrt(T(2.0) * M_PIl));
}

/** Restore a circular beam*/
void ModelRenderer::Restore(double* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double beamSize, long double startFrequency, long double endFrequency, PolarizationEnum polarization)
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
			if(xRight < xLeft) xRight = xLeft;
			if(xRight > (int) imageWidth) xRight = (int) imageWidth;
			if(yTop < 0) yTop = 0;
			if(yTop > (int) imageHeight) yTop = (int) imageHeight;
			if(yBottom < yTop) yBottom = yTop;
			if(yBottom > (int) imageHeight) yBottom = (int) imageHeight;
			
			for(int y=yTop; y!=yBottom; ++y)
			{
				double *imageDataPtr = imageData + y*imageWidth+xLeft;
				for(int x=xLeft; x!=xRight; ++x)
				{
					long double l, m;
					ImageCoordinates::XYToLM<long double>(x, y, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, l, m);
					l += _phaseCentreDL; m += _phaseCentreDM;
					long double dist = sqrt((l-sourceL)*(l-sourceL) + (m-sourceM)*(m-sourceM));
					long double g = gaus(dist, sigma);
					(*imageDataPtr) += double(g * intFlux);
					++imageDataPtr;
				}
			}
		}
	}
}

void ModelRenderer::renderGaussianComponent(double* imageData, size_t imageWidth, size_t imageHeight, long double posRA, long double posDec, long double gausMaj, long double gausMin, long double gausPA, long double flux, bool normalizeIntegratedFlux)
{
	// Using the FWHM formula for a Gaussian:
	long double sigmaMaj = gausMaj / (2.0L * sqrtl(2.0L * logl(2.0L)));
	long double sigmaMin = gausMin / (2.0L * sqrtl(2.0L * logl(2.0L)));
	// TODO this won't work for non-equally spaced dimensions
	long double minPixelScale = std::min(_pixelScaleL, _pixelScaleM);
	// ...And this is not accurate, as the correlation between beam and source PA has
	// to be taken into account
	if(normalizeIntegratedFlux)
	{
		double factor = 2.0L * M_PI * sigmaMaj * sigmaMin / (minPixelScale * minPixelScale);
		if(factor > 1.0)
			flux /= factor;
	}
	
	// Make rotation matrix
	long double transf[4];
	// Position angle is angle from North: 
	long double angle = gausPA+0.5*M_PI;
	transf[2] = sin(angle);
	transf[0] = cos(angle);
	transf[1] = -transf[2];
	transf[3] = transf[0];
	double sigmaMax = std::max(std::fabs(sigmaMaj * transf[0]), std::fabs(sigmaMaj * transf[1]));
	// Multiply with scaling matrix to make variance 1.
	transf[0] = transf[0] / sigmaMaj;
	transf[1] = transf[1] / sigmaMaj;
	transf[2] = transf[2] / sigmaMin;
	transf[3] = transf[3] / sigmaMin;
	int boundingBoxSize = ceil(sigmaMax * 20.0 / minPixelScale);
	long double sourceL, sourceM;
	ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec, sourceL, sourceM);

	// Calculate the bounding box
	int sourceX, sourceY;
	ImageCoordinates::LMToXY<long double>(sourceL-_phaseCentreDL, sourceM-_phaseCentreDM, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);
	int
		xLeft = sourceX - boundingBoxSize,
		xRight = sourceX + boundingBoxSize,
		yTop = sourceY - boundingBoxSize,
		yBottom = sourceY + boundingBoxSize;
	if(xLeft < 0) xLeft = 0;
	if(xLeft > (int) imageWidth) xLeft = (int) imageWidth;
	if(xRight < xLeft) xRight = xLeft;
	if(xRight > (int) imageWidth) xRight = (int) imageWidth;
	if(yTop < 0) yTop = 0;
	if(yTop > (int) imageHeight) yTop = (int) imageHeight;
	if(yBottom < yTop) yBottom = yTop;
	if(yBottom > (int) imageHeight) yBottom = (int) imageHeight;
	
	for(int y=yTop; y!=yBottom; ++y)
	{
		double *imageDataPtr = imageData + y*imageWidth+xLeft;
		for(int x=xLeft; x!=xRight; ++x)
		{
			long double l, m;
			ImageCoordinates::XYToLM<long double>(x, y, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, l, m);
			l += _phaseCentreDL; m += _phaseCentreDM;
			long double
				lTransf = (l-sourceL)*transf[0] + (m-sourceM)*transf[1],
				mTransf = (l-sourceL)*transf[2] + (m-sourceM)*transf[3];
			long double dist = sqrt(lTransf*lTransf + mTransf*mTransf);
			long double g = gaus(dist, 1.0L);
			(*imageDataPtr) += double(g * flux);
			++imageDataPtr;
		}
	}
}

void ModelRenderer::renderPointComponent(double* imageData, size_t imageWidth, size_t imageHeight, long double posRA, long double posDec, long double flux)
{
	long double sourceL, sourceM;
	ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec, sourceL, sourceM);
	sourceL -= _phaseCentreDL; sourceM -= _phaseCentreDM;
	
	int sourceX, sourceY;
	ImageCoordinates::LMToXY<long double>(sourceL, sourceM, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);
	
	if(sourceX >= 0 && sourceX < (int) imageWidth && sourceY >= 0 && sourceY < (int) imageHeight)
	{
		double *imageDataPtr = imageData + sourceY*imageWidth + sourceX;
		(*imageDataPtr) += double(flux);
	}
}

/** Restore a model with an elliptical beam */
void ModelRenderer::Restore(double* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double beamMaj, long double beamMin, long double beamPA, long double startFrequency, long double endFrequency, PolarizationEnum polarization)
{
	ao::uvector<double> renderedWithoutBeam(imageWidth * imageHeight, 0.0);
	renderModel(renderedWithoutBeam.data(), imageWidth, imageHeight, model, startFrequency, endFrequency, polarization, true);
	Restore(imageData, renderedWithoutBeam.data(), imageWidth, imageHeight, beamMaj, beamMin, beamPA, _pixelScaleL, _pixelScaleM);
}

/**
 * Restore a diffuse image (e.g. produced with multi-scale clean)
 */
void ModelRenderer::Restore(double* imageData, const double* modelData, size_t imageWidth, size_t imageHeight, long double beamMaj, long double beamMin, long double beamPA, long double pixelScaleL, long double pixelScaleM)
{
	if(beamMaj == 0.0 && beamMin == 0.0)
	{
		for(size_t j=0; j!=imageWidth*imageHeight; ++j)
			imageData[j] += modelData[j];
	}
	else {
		// Using the FWHM formula for a Gaussian:
		long double sigmaMaj = beamMaj / (2.0L * sqrtl(2.0L * logl(2.0L)));
		long double sigmaMin = beamMin / (2.0L * sqrtl(2.0L * logl(2.0L)));
		
		// Make rotation matrix
		long double transf[4];
		// Position angle is angle from North: 
		long double angle = beamPA+0.5*M_PI;
		transf[2] = sin(angle);
		transf[0] = cos(angle);
		transf[1] = -transf[2];
		transf[3] = transf[0];
		double sigmaMax = std::max(std::fabs(sigmaMaj * transf[0]), std::fabs(sigmaMaj * transf[1]));
		// Multiple with scaling matrix to make variance 1.
		transf[0] = transf[0] / sigmaMaj;
		transf[1] = transf[1] / sigmaMaj;
		transf[2] = transf[2] / sigmaMin;
		transf[3] = transf[3] / sigmaMin;
		
		size_t minDimension = std::min(imageWidth, imageHeight);
		size_t boundingBoxSize = std::min<size_t>(ceil(sigmaMax * 40.0 / std::min(pixelScaleL, pixelScaleM)), minDimension);
		if(boundingBoxSize%2!=0) ++boundingBoxSize;
		ao::uvector<double> kernel(boundingBoxSize*boundingBoxSize);
		typename ao::uvector<double>::iterator i=kernel.begin();
		for(size_t y=0; y!=boundingBoxSize; ++y)
		{
			for(size_t x=0; x!=boundingBoxSize; ++x)
			{
				long double l, m;
				ImageCoordinates::XYToLM<long double>(x, y, pixelScaleL, pixelScaleM, boundingBoxSize, boundingBoxSize, l, m);
				long double
					lTransf = l*transf[0] + m*transf[1],
					mTransf = l*transf[2] + m*transf[3];
				long double dist = sqrt(lTransf*lTransf + mTransf*mTransf);
				*i = gaus(dist, (long double) 1.0);
				++i;
			}
		}
		
		ao::uvector<double> convolvedModel(imageWidth*imageHeight);
		memcpy(convolvedModel.data(), modelData, sizeof(double)*imageWidth*imageHeight);
		
		FFTConvolver::Convolve(convolvedModel.data(), imageWidth, imageHeight, kernel.data(), boundingBoxSize);
		for(size_t j=0; j!=imageWidth*imageHeight; ++j)
			imageData[j] += convolvedModel[j];
	}
}

/**
* Render without beam convolution, such that each point-source is one pixel.
*/
void ModelRenderer::renderModel(double* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double startFrequency, long double endFrequency, PolarizationEnum polarization, bool normalizeIntegratedFlux)
{
	for(Model::const_iterator src=model.begin(); src!=model.end(); ++src)
	{
		for(ModelSource::const_iterator comp=src->begin(); comp!=src->end(); ++comp)
		{
			const long double
				posRA = comp->PosRA(),
				posDec = comp->PosDec(),
				intFlux = comp->SED().IntegratedFlux(startFrequency, endFrequency, polarization);
			
			if(comp->Type() == ModelComponent::GaussianSource)
			{
				long double
					gausMaj = comp->MajorAxis(),
					gausMin = comp->MinorAxis(),
					gausPA = comp->PositionAngle();
				renderGaussianComponent(imageData, imageWidth, imageHeight, posRA, posDec, gausMaj, gausMin, gausPA, intFlux, normalizeIntegratedFlux);
			}
			else
				renderPointComponent(imageData, imageWidth, imageHeight, posRA, posDec, intFlux);
		}
	}
}

void ModelRenderer::RenderInterpolatedSource(double* image, size_t width, size_t height, double flux, double x, double y)
{
	ao::uvector<double>
		hSinc(std::min<size_t>(width,128)+1),
		vSinc(std::min<size_t>(height,128)+1);

	int midH = hSinc.size()/2;
	int midV = vSinc.size()/2;
	
	double xr = x - floor(x);
	double yr = y - floor(y);
	
	for(size_t i=0; i!=hSinc.size(); ++i)
	{
		double xi = (int(i) - midH - xr) * M_PI;
		if(xi == 0.0)
			hSinc[i] = 1.0;
		else
			hSinc[i] = sin(xi) / xi;
	}
	for(size_t i=0; i!=vSinc.size(); ++i)
	{
		double yi = (int(i) - midV - yr) * M_PI;
		if(yi == 0.0)
			vSinc[i] = 1.0;
		else
			vSinc[i] = sin(yi) / yi;
	}
	
	size_t
		xOffset = floor(x) - midH,
		yOffset = floor(y) - midV,
		startX = std::max<int>(xOffset, 0),
		startY = std::max<int>(yOffset, 0),
		endX = std::min<size_t>(xOffset + hSinc.size(), width),
		endY = std::min<size_t>(yOffset + vSinc.size(), height);
	for(size_t yi=startY; yi!=endY; ++yi)
	{
		double* ptr = &image[yi * width];
		double vFlux = flux * vSinc[yi - yOffset];
		for(size_t xi=startX; xi!=endX; ++xi)
		{
			ptr[xi] += vFlux * hSinc[xi - xOffset];
		}
	}
}
