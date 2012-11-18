
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

template void ModelRenderer::Render(double* imageData, size_t imageSize, const Model& model, long double beamSize, long double startFrequency, long double endFrequency);

template<typename NumType>
void ModelRenderer::Render(NumType* imageData, size_t imageSize, const Model& model, long double beamSize, long double startFrequency, long double endFrequency)
{
	long double midX = (long double) imageSize / 2.0, midY = (long double) imageSize / 2.0;
	for(Model::const_iterator src=model.begin(); src!=model.end(); ++src)
	{
		long double
			posRA = src->PosRA(),
			posDec = src->PosDec(),
			sourceL, sourceM;
		ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec, sourceL, sourceM);
		const SourceStrength<long double> &brightness = src->Brightness();
		NumType *imageDataPtr = imageData;
		const long double intFlux = brightness.IntegratedFlux(startFrequency, endFrequency);
		
		//std::cout << "Source: " << src->PosRA() << "," << src->PosDec() << " Phase centre: " << _phaseCentreRA << "," << _phaseCentreDec << " beamsize: " << beamSize << "\n";
		std::cout << "Adding source " << src->Name() << " at " << sourceL << "," << sourceM << " of "
			<< intFlux << " Jy ("
			<< startFrequency/1000000.0 << "-" << endFrequency/1000000.0 << " MHz).\n";
		for(size_t y=0; y!=imageSize; ++y)
		{
			for(size_t x=0; x!=imageSize; ++x)
			{
				long double l = ((NumType) x - midX) * -_pixelScaleL;
				long double m = (midY - (NumType) y) * -_pixelScaleM;
				long double dist = sqrt((l-sourceL)*(l-sourceL) + (m-sourceM)*(m-sourceM));
				long double g = gaus(dist, beamSize);
				(*imageDataPtr) += NumType(g * g * intFlux);
				++imageDataPtr;
			}
		}
	}
}

