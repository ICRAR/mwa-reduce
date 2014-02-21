#include "fitsreader.h"
#include "fitswriter.h"
#include "matrix2x2.h"

#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include <complex>
#include <limits>

void read(const FitsReader& templateReader, const char* filename, std::vector<double>& data)
{
	FitsReader inpReader(filename);
	size_t
		width = inpReader.ImageWidth(),
		height = inpReader.ImageHeight();
		
	if(width != templateReader.ImageWidth() || height != templateReader.ImageHeight())
		throw std::runtime_error("Not all images had same size");
	data.resize(width * height);
	inpReader.Read<double>(&data[0]);
}

int main(int argc, char *argv[])
{
	if(argc != 17)
		std::cout << "Syntax:\n"
			"pbcorrect <xx.fits> <xy.fits> <xyi.fits> <yx.fits> <beamxx.fits> <beamxxi.fits> <beamxy.fits> <beamxyi.fits> <beamyx.fits> <beamyxi.fits> <beamyy.fits> <beamyyi.fits> <outi.fits> <outq.fits> <outu.fits> <outv.fits>\n";
		
	const char
		*inpFilename[4] = { argv[1], argv[2], argv[3], argv[4] },
		*beamFilename[8] = { argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], argv[11], argv[12] },
		*outFilename[4] = { argv[13], argv[14], argv[15], argv[16] };
	
	FitsReader firstImage(argv[1]);
	std::vector<double> inputData[4], beamData[8];
	for(size_t i=0; i!=4; ++i)
		read(firstImage, inpFilename[i], inputData[i]);
	for(size_t i=0; i!=8; ++i)
		read(firstImage, beamFilename[i], beamData[i]);
	
	const size_t 
		width = firstImage.ImageWidth(),
		height = firstImage.ImageHeight(),
		imgSize = width * height;
	for(size_t i=0; i!=imgSize; ++i)
	{
		std::complex<double> imgValues[4], beamValues[4];
		imgValues[0] = inputData[0][i];
		imgValues[1] = std::complex<double>(inputData[1][i], inputData[2][i]);
		imgValues[2] = std::conj(imgValues[1]);
		imgValues[3] = inputData[3][i];
		
		beamValues[0] = std::complex<double>(beamData[0][i], beamData[1][i]);
		beamValues[1] = std::complex<double>(beamData[2][i], beamData[3][i]);
		beamValues[2] = std::complex<double>(beamData[4][i], beamData[5][i]);
		beamValues[3] = std::complex<double>(beamData[6][i], beamData[7][i]);
		
		if(Matrix2x2::Invert(beamValues))
		{
			std::complex<double> tempValues[4];
			Matrix2x2::ATimesB(tempValues, beamValues, imgValues);
			Matrix2x2::ATimesHermB(imgValues, tempValues, beamValues);
		}
		else {
			for(size_t p=0; p!=4; ++p)
				inputData[p][i] = std::numeric_limits<double>::quiet_NaN();
		}
			
		inputData[0][i] = 0.5 * (imgValues[0].real() + imgValues[3].real());
		inputData[1][i] = 0.5 * (imgValues[0].real() - imgValues[3].real());
		inputData[2][i] = 0.5 * (imgValues[1].real() + imgValues[2].real());
		inputData[3][i] = 0.5 * (-imgValues[1].imag() + imgValues[2].imag());
	}
	
	FitsWriter writer(firstImage);
	writer.SetPolarization(Polarization::StokesI);
	writer.Write<double>(outFilename[0], &inputData[0][0]);
	writer.SetPolarization(Polarization::StokesQ);
	writer.Write<double>(outFilename[1], &inputData[1][0]);
	writer.SetPolarization(Polarization::StokesU);
	writer.Write<double>(outFilename[2], &inputData[2][0]);
	writer.SetPolarization(Polarization::StokesV);
	writer.Write<double>(outFilename[3], &inputData[3][0]);
}
