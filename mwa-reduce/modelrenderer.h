#ifndef MODELRENDERER_H
#define MODELRENDERER_H

#include <cstring>
#include "polarizationenum.h"

class ModelRenderer
{
	public:
		ModelRenderer(long double phaseCentreRA, long double phaseCentreDec, long double pixelScaleL, long double pixelScaleM, long double phaseCentreDL=0.0, long double phaseCentreDM=0.0) :
			_phaseCentreRA(phaseCentreRA), _phaseCentreDec(phaseCentreDec), _pixelScaleL(pixelScaleL), _pixelScaleM(pixelScaleM), _phaseCentreDL(phaseCentreDL), _phaseCentreDM(phaseCentreDM)
		{
		}
		
		template<typename NumType>
		void Restore(NumType* imageData, size_t imageWidth, size_t imageHeight, const class Model& model, long double beamSize, long double startFrequency, long double endFrequency, PolarizationEnum polarization);
		
		template<typename NumType>
		void RenderModel(NumType* imageData, size_t imageWidth, size_t imageHeight, const class Model& model, long double startFrequency, long double endFrequency, PolarizationEnum polarization);
	private:
		long double _phaseCentreRA;
		long double _phaseCentreDec;
		long double _pixelScaleL, _pixelScaleM;
		long double _phaseCentreDL, _phaseCentreDM;
		template<typename T>
		T gaus(T x, T sigma) const;
		
		ModelRenderer(const ModelRenderer &) { }
		void operator=(const ModelRenderer &) { };
};

#endif
