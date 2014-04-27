
#include "cleanalgorithms/multiscaleclean.h"
#include "fitswriter.h"
#include "fftconvolver.h"
#include "fftresampler.h"

void testSomeSpecifics(const std::string& imageName, const std::string& psfName)
{
	double scaleA = 100.0;
	
	ao::uvector<double> shape;
	size_t n;
	MultiScaleClean<>::MakeShapeFunction(200.0, shape, n);
	FitsWriter fitsWriter;
	fitsWriter.SetImageDimensions(n, n);
	fitsWriter.Write("testms-shape.fits", shape.data());
	
	size_t width = 1024, height = 1024;
	fitsWriter.SetImageDimensions(width, height);
	
	ao::uvector<double> kernel(width * height, 0.0);
	FFTConvolver::PrepareKernel(kernel.data(), width, height, shape.data(), n);
	fitsWriter.Write("testms-kernel.fits", kernel.data());
	
	ao::uvector<double> image(width * height, 0.0);
	image[width/2 + width*(height/2)] = 1;
	image[100 + 100*width] = 1;
	FFTConvolver::ConvolveSameSize(image.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-convolved.fits", image.data());
	
	FFTConvolver::ConvolveSameSize(image.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-convolved2x.fits", image.data());
	
	FitsReader imgReader(imageName), psfReader(psfName);
	width = imgReader.ImageWidth();
	height = imgReader.ImageHeight();
	ao::uvector<double> psf(width*height);
	image.resize(width*height);
	imgReader.Read(image.data());
	psfReader.Read(psf.data());
	
	size_t rescaledWidth = width/4, rescaledHeight = height/4;
	ao::uvector<double> rescaledImageA(rescaledWidth * rescaledHeight), rescaledPsf(rescaledWidth * rescaledHeight);
	FFTResampler resampler(width, height, rescaledWidth, rescaledHeight, 2);
	resampler.AddTask(image.data(), rescaledImageA.data());
	resampler.AddTask(psf.data(), rescaledPsf.data());
	resampler.Start();
	resampler.Finish();
	ao::uvector<double> rescaledImageB(rescaledImageA);
	
	width = rescaledWidth; height = rescaledHeight;
	kernel.resize(width*height);
	
	MultiScaleClean<>::MakeShapeFunction(scaleA, shape, n);
	FFTConvolver::PrepareKernel(kernel.data(), width, height, shape.data(), n);
	
	fitsWriter.SetImageDimensions(width, height);
	fitsWriter.Write("testms-input.fits", rescaledImageA.data());
	fitsWriter.Write("testms-psf.fits", rescaledPsf.data());
	
	FFTConvolver::ConvolveSameSize(rescaledImageA.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-input-convolvedA.fits", rescaledImageA.data());
	FFTConvolver::ConvolveSameSize(rescaledPsf.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-psf-convolvedA.fits", rescaledPsf.data());
	FFTConvolver::ConvolveSameSize(rescaledPsf.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-psf-convolvedA2x.fits", rescaledPsf.data());
	
	MultiScaleClean<>::MakeShapeFunction(scaleA*0.5, shape, n);
	FFTConvolver::PrepareKernel(kernel.data(), width, height, shape.data(), n);
	
	FFTConvolver::ConvolveSameSize(rescaledImageB.data(), kernel.data(), width, height);
	fitsWriter.Write("testms-input-convolvedB.fits", rescaledImageB.data());
}

void testClean(const std::string& imageName, const std::string& psfName)
{
	FitsReader imgReader(imageName), psfReader(psfName);
	
	/**
	 * The real cleaning
	 */
	size_t width = imgReader.ImageWidth();
	size_t height = imgReader.ImageHeight();
	ao::uvector<double> psf(width*height);
	psfReader.Read(psf.data());
	
	ImageBufferAllocator<double> allocator;
	MultiScaleClean<clean_algorithms::SingleImageSet> msClean(imgReader.BeamMinorAxisRad(), imgReader.PixelSizeX(), imgReader.PixelSizeY());
	clean_algorithms::SingleImageSet
		imageSet(width*height, allocator),
		modelSet(width*height, allocator);
	imgReader.Read(imageSet.GetImage(0));
	memset(modelSet.GetImage(0), 0, sizeof(double)*width*height);
	std::vector<double*> psfs(1, psf.data());
	
	bool reachedStopGain = false;
	msClean.SetSubtractionGain(0.1);
	msClean.SetThreshold(0.8);
	msClean.SetMaxNIter(10000);
	msClean.ExecuteMajorIteration(imageSet, modelSet, psfs, width, height, reachedStopGain);
	
	FitsWriter writer;
	writer.SetImageDimensions(width, height);
	writer.Write("multiscale-model.fits", modelSet.GetImage(0));
	writer.Write("multiscale-residual.fits", imageSet.GetImage(0));
}

void testMultiClean(const std::string& imageName, const std::string& psfName)
{
	FitsReader imgReader(imageName), psfReader(psfName);
	
	/**
	 * The real cleaning
	 */
	size_t width = imgReader.ImageWidth();
	size_t height = imgReader.ImageHeight();
	ao::uvector<double> psf(width*height);
	psfReader.Read(psf.data());
	
	ImageBufferAllocator<double> allocator;
	MultiScaleClean<clean_algorithms::MultiImageSet> msClean(imgReader.BeamMinorAxisRad(), imgReader.PixelSizeX(), imgReader.PixelSizeY());
	clean_algorithms::MultiImageSet
		imageSet(width*height, 2, allocator),
		modelSet(width*height, 2, allocator);
	for(size_t i=0; i!=8; ++i)
	{
		imgReader.Read(imageSet.GetImage(i));
		memset(modelSet.GetImage(i), 0, sizeof(double)*width*height);
	}
	std::vector<double*> psfs(2, psf.data());
	
	bool reachedStopGain = false;
	msClean.SetSubtractionGain(0.1);
	msClean.SetThreshold(0.8);
	msClean.SetMaxNIter(10000);
	msClean.ExecuteMajorIteration(imageSet, modelSet, psfs, width, height, reachedStopGain);
	
	FitsWriter writer;
	writer.SetImageDimensions(width, height);
	writer.Write("multiscale-model.fits", modelSet.GetImage(0));
	writer.Write("multiscale-residual.fits", imageSet.GetImage(0));
}

int main(int argc, char* argv[])
{
	const std::string& imageName(argv[1]), psfName(argv[2]);
	
	//testSomeSpecifics(imageName, psfName);
	//testClean(imageName, psfName);
	testMultiClean(imageName, psfName);
	
}
