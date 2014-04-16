#include <iostream>
#include <vector>
#include "tilebeam.h"
#include "../fitswriter.h"
#include "../matrix2x2.h"

/**
	* Make azimuth and zenith angle arrays for a square image of side npix
	* Projection is sine, all-sky. Returns (az,za). Angles are in radian.
*/
void makeAZZA(size_t npix, double* zas, double* azs)
{
	// build az and za arrays
	for(size_t yi=0; yi!=npix; ++yi)
	{
		for(size_t xi=0; xi!=npix; ++xi)
		{
			size_t index = xi + yi*npix;
			double x = yi - npix*0.5;
			double y = xi - npix*0.5;
			double d = sqrt(x*x + y*y)/(0.5*npix);
			// only select pixels above horizon
			if(d <= 1.0)
				zas[index] = asin(d);
			else
				zas[index] = M_PI/2.0;
			azs[index] = atan2(y, x);
		}
	}
}

// plot the output of tile Jones matrices
void plotArrayJones(size_t npix, const std::vector<std::complex<double>>& js, double za, double freq, bool invert)
{
	const size_t size = npix*npix;
	
	std::vector<std::complex<double>> copy(js);
	std::string extension;
	if(invert)
	{
		for(size_t i=0; i!=size; ++i)
		{
			if(!Matrix2x2::Invert(&copy[i*4]))
			{
				for(size_t p=0; p!=4; ++p) copy[i*4+p] = std::numeric_limits<double>::quiet_NaN();
			}
		}
		extension = "-inv.fits";
	}
	else extension = ".fits";
	
	std::vector<double> image(size);
	FitsWriter writer;
	writer.SetImageDimensions(npix, npix, 0, 0, 0, 0);
	
	for(size_t i=0; i!=size; ++i)
		image[i] = std::abs(copy[i*4 + 0]);
	std::ostringstream s1;
	s1 << "/tmp/MWA_J00_voltage_mag_" << (freq/1e6) << "MHz_ZA" << za << extension;
	writer.Write(s1.str(), image.data());

	for(size_t i=0; i!=size; ++i)
		image[i] = std::abs(copy[i*4 + 1]);
	std::ostringstream s2;
	s2 << "/tmp/MWA_J01_voltage_mag_" << (freq/1e6) << "MHz_ZA" << za << extension;
	writer.Write(s2.str(), image.data());

	for(size_t i=0; i!=size; ++i)
		image[i] = std::abs(copy[i*4 + 2]);
	std::ostringstream s3;
	s3 << "/tmp/MWA_J10_voltage_mag_" << (freq/1e6) << "MHz_ZA" << za << extension;
	writer.Write(s3.str(), image.data());

	for(size_t i=0; i!=size; ++i)
		image[i] = std::abs(copy[i*4 + 3]);
	std::ostringstream s4;
	s4 << "/tmp/MWA_J11_voltage_mag_" << (freq/1e6) << "MHz_ZA" << za << extension;
	writer.Write(s4.str(), image.data());
}

void plotVisResponse(size_t npix, const std::vector<std::complex<double>>& js, double za, double freq)
{
	size_t size = npix*npix;
	std::vector<double> image(size);
	FitsWriter writer;
	writer.SetImageDimensions(npix, npix, 0, 0, 0, 0);
	
	std::vector<std::complex<double>> response(size*4);
	for(size_t i=0; i!=size; ++i)
	{
		size_t xx = i*4, xy = xx+1, yx=xx+2, yy=xx+3;
    response[xx] = js[xx]*std::conj(js[xx]) + js[xy]*std::conj(js[xy]);
    response[yy] = js[yx]*std::conj(js[yx]) + js[yy]*std::conj(js[yy]);
    response[xy] = js[xx]*std::conj(js[yx]) + js[xy]*std::conj(js[yy]);
    response[yx] = js[yx]*std::conj(js[xx]) + js[yy]*std::conj(js[xy]);
	}
	
	for(size_t i=0; i!=size; ++i)
		image[i] = std::abs(response[i*4 + 0]);
	std::ostringstream s1;
	s1 << "/tmp/MWA_XX_mag_" << (freq/1e6) << "MHz_ZA" << za << ".fits";
	writer.Write(s1.str(), image.data());

	for(size_t i=0; i!=size; ++i)
		image[i] = std::abs(response[i*4 + 1]);
	std::ostringstream s2;
	s2 << "/tmp/MWA_XY_mag_" << (freq/1e6) << "MHz_ZA" << za << ".fits";
	writer.Write(s2.str(), image.data());

	for(size_t i=0; i!=size; ++i)
		image[i] = std::abs(response[i*4 + 2]);
	std::ostringstream s3;
	s3 << "/tmp/MWA_YX_mag_" << (freq/1e6) << "MHz_ZA" << za << ".fits";
	writer.Write(s3.str(), image.data());

	for(size_t i=0; i!=size; ++i)
		image[i] = std::abs(response[i*4 + 3]);
	std::ostringstream s4;
	s4 << "/tmp/MWA_YY_mag_" << (freq/1e6) << "MHz_ZA" << za << ".fits";
	writer.Write(s4.str(), image.data());
}

int main(int argc, char* argv[])
{
	double frequencies[2] = { 216e6, 88e6 };
	for(size_t freqIndex = 0; freqIndex!=2; ++freqIndex)
	{
		double lat = -26.7;
		double freq = frequencies[freqIndex];
		const size_t npix = 512;
		
		std::vector<double> zas(npix*npix), azs(npix*npix);
		makeAZZA(npix, zas.data(), azs.data());

		// parallactic angle
		//(ha,dec) = h2e(az,za,lat)
		//pa = calcParallacticAngle(ha,dec,lat)

		//plotDipoleJones(freq);

		// tests with full tile
		// az=0, za=14 degs
		std::vector<double> delays[3];
		delays[0] = std::vector<double>{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		delays[1] = std::vector<double>{6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0,6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0};
		delays[2] = std::vector<double>{12,12,12,12,8,8,8,8,4,4,4,4,0,0,0,0,12,12,12,12,8,8,8,8,4,4,4,4,0,0,0,0};
		double zaOfDelays[3] = {0, 14, 28};
		
		for(size_t i=0; i!=3; ++i)
		{
			std::cout << "ZA is: " << zaOfDelays[i] << ". First delay is: " << delays[i][0] << '\n';
			
			TileBeam2014 tileBeam(delays[i].data());
			//tileBeam.ArrayResponse(zas[i], az, freq, gains);
			
			//std::cout << "plotting Array factor voltage for ZA " << zaOfDelays[i] << '\n';
			//plotArrayFactors(ax0,ay0,za_delay);
			
			std::cout << "plotting tile Jones response for ZA " << zaOfDelays[i] << '\n';
			std::vector<std::complex<double>> js(4 * npix*npix);
			for(size_t j=0; j!=npix*npix; ++j)
				tileBeam.ArrayResponse(zas[j], azs[j], freq, &js[j*4]);
			plotArrayJones(npix, js, zaOfDelays[i], freq, false);
			std::cout << "plotting inverse tile Jones response for ZA " << zaOfDelays[i] << '\n';
			plotArrayJones(npix, js, zaOfDelays[i], freq, true);
			std::cout << "Plotting visibility response for two identical tiles ZA " << zaOfDelays[i] << '\n';
			plotVisResponse(npix, js, zaOfDelays[i], freq);
		}
	}
}
