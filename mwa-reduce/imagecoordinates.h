#ifndef IMAGE_COORDINATES_H
#define IMAGE_COORDINATES_H

#include <cmath>

/**
 * This class collects all the LM coordinate transform as defined in
 * Perley (1999)'s "imaging with non-coplaner arrays".
 */
class ImageCoordinates
{
	public:
		template<typename T>
		static void RaDecToLM(T ra, T dec, T phaseCentreRa, T phaseCentreDec, T &destL, T &destM)
		{
			const T
				deltaAlpha = ra - phaseCentreRa,
				sinDeltaAlpha = sin(deltaAlpha),
				cosDeltaAlpha = cos(deltaAlpha),
				sinDec = sin(dec),
				cosDec = cos(dec),
				sinDec0 = sin(phaseCentreDec),
				cosDec0 = cos(phaseCentreDec);
			
			destL = cosDec * sinDeltaAlpha;
			destM = sinDec*cosDec0 - cosDec*sinDec0*cosDeltaAlpha;
		}
		
		template<typename T>
		static T RaDecToN(T ra, T dec, T phaseCentreRa, T phaseCentreDec)
		{
			const T
				cosDeltaAlpha = cos(ra - phaseCentreRa),
				sinDec = sin(dec),
				cosDec = cos(dec),
				sinDec0 = sin(phaseCentreDec),
				cosDec0 = cos(phaseCentreDec);
			
			return sinDec*sinDec0 + cosDec*cosDec0*cosDeltaAlpha;
		}
		
		template<typename T>
		static void LMToRaDec(T l, T m, T phaseCentreRa, T phaseCentreDec, T &destRa, T &destDec)
		{
			const T
				cosDec0 = cos(phaseCentreDec),
				sinDec0 = sin(phaseCentreDec),
				lmTerm = sqrt((T) 1.0 - l*l - m*m),
				deltaAlpha = atan2(l, lmTerm*cosDec0 - m*sinDec0);
				
			destRa = deltaAlpha + phaseCentreRa;
			destDec = asin(m*cosDec0 + lmTerm*sinDec0);
		}
		
	private:
		ImageCoordinates();
};

#endif
