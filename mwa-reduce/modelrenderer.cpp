
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

template void ModelRenderer::Render(double* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double beamSize, long double startFrequency, long double endFrequency);

template<typename NumType>
void ModelRenderer::Render(NumType* imageData, size_t imageWidth, size_t imageHeight, const Model& model, long double beamSize, long double startFrequency, long double endFrequency)
{
	int boundingBoxSize = ceil(beamSize * 5.0 / (0.5 * (_pixelScaleL + _pixelScaleM)));
	std::cout << "Render bounding box: " << boundingBoxSize << '\n';
	for(Model::const_iterator src=model.begin(); src!=model.end(); ++src)
	{
		long double
			posRA = src->PosRA(),
			posDec = src->PosDec(),
			sourceL, sourceM;
		ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec, sourceL, sourceM);
		const SourceSDF<long double> &brightness = src->Brightness();
		const long double intFlux = brightness.IntegratedFlux(startFrequency, endFrequency);
		
		//std::cout << "Source: " << src->PosRA() << "," << src->PosDec() << " Phase centre: " << _phaseCentreRA << "," << _phaseCentreDec << " beamsize: " << beamSize << "\n";
			
		int sourceX, sourceY;
		ImageCoordinates::LMToXY<long double>(sourceL, sourceM, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);
		std::cout << "Adding source " << src->Name() << " at " << sourceX << "," << sourceY << " of "
			<< intFlux << " Jy ("
			<< startFrequency/1000000.0 << "-" << endFrequency/1000000.0 << " MHz).\n";
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
				long double dist = sqrt((l-sourceL)*(l-sourceL) + (m-sourceM)*(m-sourceM));
				long double g = gaus(dist, beamSize);
				(*imageDataPtr) += NumType(g * g * intFlux);
				++imageDataPtr;
			}
		}
	}
}

