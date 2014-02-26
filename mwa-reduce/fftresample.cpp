#include "fftresampler.h"
#include "fitsreader.h"

#include "uvector.h"
#include "fitswriter.h"

#include <iostream>

#include <unistd.h>

int main(int argc, char* argv[])
{
	if(argc < 5)
	{
		std::cout << "Syntax: fftresample <input> <output> <width> <height>\n";
	}
	else {
		FitsReader reader(argv[1]);
		const size_t width = reader.ImageWidth(), height = reader.ImageHeight();
		const size_t outWidth = atoi(argv[3]), outHeight = atoi(argv[4]);
		double* inputData = reinterpret_cast<double*>(fftw_malloc(width*height*sizeof(double)));
		double* outputData = reinterpret_cast<double*>(fftw_malloc(outWidth*outHeight*sizeof(double)));
		FFTResampler resampler(width, height, outWidth, outHeight, sysconf(_SC_NPROCESSORS_ONLN));
		reader.Read(inputData);
		
		resampler.RunSingle(inputData, outputData);
		fftw_free(inputData);
		
		FitsWriter writer(reader);
		writer.SetImageDimensions(outWidth, outHeight);
		writer.Write(argv[2], outputData);
		
		fftw_free(outputData);
	}
}
