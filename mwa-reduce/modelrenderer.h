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
		
		/**
		 * Restore with circular beam
		 */
		void Restore(double* imageData, size_t imageWidth, size_t imageHeight, const class Model& model, long double beamSize, long double startFrequency, long double endFrequency, PolarizationEnum polarization);
		
		/**
		 * Restore with elliptical beam
		 */
		void Restore(double* imageData, size_t imageWidth, size_t imageHeight, const class Model& model, long double beamMaj, long double beamMin, long double beamPA, long double startFrequency, long double endFrequency, PolarizationEnum polarization);
		
		void Restore(double* imageData, double* modelData, size_t imageWidth, size_t imageHeight, long double beamMaj, long double beamMin, long double beamPA);

		/**
		 * Render each point-source as one pixel
		 */
		void RenderModel(double* imageData, size_t imageWidth, size_t imageHeight, const class Model& model, long double startFrequency, long double endFrequency, PolarizationEnum polarization);
		
		/**
		 * This will render a source and sinc-interpolate it so it
		 * can be on non-integer positions.
		 */
		static void RenderInterpolatedSource(double* image, size_t width, size_t height, double flux, double x, double y);
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
