#ifndef MODELRENDERER_H
#define MODELRENDERER_H

#include <cstring>

class ModelRenderer
{
	public:
		ModelRenderer(long double phaseCentreRA, long double phaseCentreDec, long double pixelScaleL, long double pixelScaleM) :
			_phaseCentreRA(phaseCentreRA), _phaseCentreDec(phaseCentreDec), _pixelScaleL(pixelScaleL), _pixelScaleM(pixelScaleM)
		{
		}
		
		template<typename NumType>
		void Render(NumType* imageData, size_t imageWidth, size_t imageHeight, const class Model& model, long double beamSize, long double startFrequency, long double endFrequency, size_t polarizationIndex);
	private:
		long double _phaseCentreRA;
		long double _phaseCentreDec;
		long double _pixelScaleL;
		long double _pixelScaleM;
		template<typename T>
		T gaus(T x, T sigma) const;
		
		ModelRenderer(const ModelRenderer &) { }
		void operator=(const ModelRenderer &) { };
};

#endif
