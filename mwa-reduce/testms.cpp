
#include "cleanalgorithms/multiscaleclean.h"
#include "fitswriter.h"
#include "fftresampler.h"

int main(int argc, char* argv[])
{
	ao::uvector<double> shape;
	size_t n;
	MultiScaleClean<>::MakeShapeFunction(200.0, shape, n);
	FitsWriter fitsWriter;
	fitsWriter.SetImageDimensions(n, n);
	fitsWriter.Write("testms-shape.fits", shape.data());
	
	size_t width = 1024, height = 1024;
	fitsWriter.SetImageDimensions(width, height);
	
	ao::uvector<double> kernel(width * height, 0.0);
	MultiScaleClean<>::PrepareKernel(kernel.data(), width, height, shape.data(), n);
	fitsWriter.Write("testms-kernel.fits", kernel.data());
	
	ao::uvector<double> image(width * height, 0.0);
	image[width/2 + width*(height/2)] = 1;
	image[100 + 100*width] = 1;
	MultiScaleClean<>::ConvolveSameSize(image.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-convolved.fits", image.data());
	
	MultiScaleClean<>::ConvolveSameSize(image.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-convolved2x.fits", image.data());
	
	const std::string& imageName(argv[1]), psfName(argv[2]);
	FitsReader imgReader(imageName), psfReader(psfName);
	width = imgReader.ImageWidth();
	height = imgReader.ImageHeight();
	ao::uvector<double> psf(width*height);
	image.resize(width*height);
	imgReader.Read(image.data());
	psfReader.Read(psf.data());
	
	size_t rescaledWidth = width/4, rescaledHeight = height/4;
	ao::uvector<double> rescaledImage(rescaledWidth * rescaledHeight), rescaledPsf(rescaledWidth * rescaledHeight);
	FFTResampler resampler(width, height, rescaledWidth, rescaledHeight, 2);
	resampler.AddTask(image.data(), rescaledImage.data());
	resampler.AddTask(psf.data(), rescaledPsf.data());
	resampler.Start();
	resampler.Finish();
	
	width = rescaledWidth; height = rescaledHeight;
	kernel.resize(width*height);
	MultiScaleClean<>::PrepareKernel(kernel.data(), width, height, shape.data(), n);
	
	fitsWriter.SetImageDimensions(width, height);
	fitsWriter.Write("testms-input.fits", rescaledImage.data());
	MultiScaleClean<>::ConvolveSameSize(rescaledImage.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-input-convolved.fits", rescaledImage.data());
	fitsWriter.Write("testms-psf.fits", rescaledPsf.data());
	MultiScaleClean<>::ConvolveSameSize(rescaledPsf.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-psf-convolved.fits", rescaledPsf.data());
	MultiScaleClean<>::ConvolveSameSize(rescaledPsf.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-psf-convolved2x.fits", rescaledPsf.data());
}
