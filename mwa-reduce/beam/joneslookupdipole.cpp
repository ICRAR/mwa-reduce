#include "joneslookupdipole.h"
#include "system.h"

#ifdef HAVE_ALGLIB

#include <fitsio.h>

#include <iostream>
#include <set>
#include <string>
#include <stdexcept>
#include <limits>

alglib::real_1d_array JonesLookupDipole::_zaValues;
alglib::real_1d_array JonesLookupDipole::_phValues;
std::map<double, JonesLookupDipole::FrequencyTable> JonesLookupDipole::_tables;

void JonesLookupDipole::loadLookupTable()
{
	std::string filename = System::FindPythonFilePath("mwapy/pb/Jmatrix.fits");
	// Load a dipole Jones response lookup table (FITS file)
	// Data are a direct conversion of the output of the simulation so has redundant info.
	// We do all the conversion and stuff here
	// Columns are:
	// theta phi  real(Jxt(t,p)) imag(Jxt(t,p)) real(Jxp(t,p)) imag(Jxp(t,p)) real(Jyt(t,p)) imag(Jyt(t,p)) real(Jyp(t,p)) imag(Jyp(t,p)))
	
	int status = 0;
	fitsfile *fitsPtr;
	fits_open_file(&fitsPtr, filename.c_str(), READONLY, &status);
	checkStatus(status, filename, "Open Jones matrix file");
	
	int nFreqs = 0;
	fits_get_num_hdus(fitsPtr, &nFreqs, &status);
	checkStatus(status, filename);
	
	int naxis;
	fits_get_img_dim(fitsPtr, &naxis, &status);
	checkStatus(status, filename);
	if(naxis != 2)
		throw std::runtime_error("Fits file had unexpected dimensions");
	
	long naxes[2];
	fits_get_img_size(fitsPtr, naxis, naxes, &status);
	checkStatus(status, filename);
	size_t rowCount = naxes[1], columnCount = naxes[0];
	size_t valueCount = rowCount * columnCount;

	std::cout << "Loading ZA/PH values...\n";
	ao::uvector<double> data(valueCount);
	long fpixel[3] = {1, 1};
	fits_read_pix(fitsPtr, TDOUBLE, fpixel, valueCount, 0, data.data(), 0, &status);
	checkStatus(status, filename);
	ao::uvector<double> zaVector(rowCount), phVector(rowCount);
	for(size_t row=0; row!=rowCount; ++row)
	{
		zaVector[row] = data[row * columnCount];
		phVector[row] = data[row * columnCount + 1];
	}
	
	std::set<double> zaSet(zaVector.begin(), zaVector.end()), phSet(phVector.begin(), phVector.end());
	zaVector.assign(zaSet.begin(), zaSet.end());
	phVector.assign(phSet.begin(), phSet.end());
	size_t uniqueZAValues = zaSet.size(), uniquePHValues = phSet.size();
	if(uniquePHValues * uniqueZAValues != rowCount)
		throw std::runtime_error("Something is wrong with the jones lookup table: nrrows != nzavalues x nphvalues");
	_zaValues.setcontent(uniqueZAValues, zaVector.data());
	_phValues.setcontent(uniquePHValues, phVector.data());
	
	for(int f=0; f!=nFreqs; ++f)
	{
		fits_movabs_hdu(fitsPtr, f+1, NULL, &status);
		checkStatus(status, filename);
	
		char keyValue[FLEN_VALUE];
		fits_read_keyword(fitsPtr, "FREQ", keyValue, NULL, &status);
		double frequency = atof(keyValue);
		std::cout << "Loading J lookup matrix for frequency " << round(frequency*100.0*1e-6)/100.0 << " MHz\n";
		if(f != 0)
			fits_read_pix(fitsPtr, TDOUBLE, fpixel, valueCount, 0, data.data(), 0, &status);
		
		FrequencyTable& table = _tables.insert(std::make_pair(frequency, FrequencyTable())).first->second;
		
		ao::uvector<double> values(uniqueZAValues * uniquePHValues * 8);
		for(size_t row=0; row!=uniqueZAValues * uniquePHValues; ++row)
		{
			for(size_t i=0; i!=8; ++i)
				values[row*8+i] = data[row * columnCount + 2 + i];
		}
		
		table.values.setcontent(uniqueZAValues*uniquePHValues*8, values.data());
		
		alglib::spline2dbuildbicubicv(
			_zaValues, uniqueZAValues,
			_phValues, uniquePHValues,
			table.values, 8,
			table.spline);

		// Determine the indices of 0 and 90 degrees of phi in the tabulated values
		size_t phZero=size_t(-1), ph90=size_t(-1), zaZero=size_t(-1);
		for(size_t phIndex=0; phIndex!=uniquePHValues; ++phIndex)
		{
			if(phVector[phIndex] == 0.0)
				phZero = phIndex;
			else if(phVector[phIndex] == 90.0)
				ph90 = phIndex;
		}
		for(size_t zaIndex=0; zaIndex!=uniqueZAValues; ++zaIndex)
		{
			if(zaVector[zaIndex] == 0.0)
				zaZero = zaIndex;
		}

		table.norm[0] = d2c(values[4*(phZero*uniqueZAValues + zaZero) + 0]);  // (zaZero, phZero)
		table.norm[1] = -d2c(values[4*(ph90*uniqueZAValues + zaZero) + 1]);   // (zaZero, ph90)
		table.norm[2] = d2c(values[4*(ph90*uniqueZAValues + zaZero) + 2]);    // (zaZero, ph90)
		table.norm[3] = d2c(values[4*(phZero*uniqueZAValues + zaZero) + 3]);  // (zaZero, phZero)
	}
}

// Return the Jones matrix for a given az, za and freq value
void JonesLookupDipole::Interpolate(std::complex<double>* jonesMatrix, double az, double za, double freq, bool zenithNorm)
{
	/* this method interpolates from the tablulated numerical results loaded
	by the constructor.*/
	
	// find the nearest freq lookup table
	std::map<double, FrequencyTable>::iterator nearest = _tables.lower_bound(freq);
	if(nearest == _tables.end())
		--nearest;
	else if(nearest != _tables.begin()) {
		std::map<double, FrequencyTable>::iterator previous = nearest;
		--previous;
		if(freq - previous->first < nearest->first - freq)
			nearest = previous;
	}
	//if(std::fabs(nearest->first - freq) > 2e6)
	//	std::cout << "Warning: Nearest tabulated impedance matrix freq is more than 2 MHz away from desired freq.\n";
	//std::cout << "Calculating beam at " << round(freq*100*1e-6)/100.0 << " from matrix for frequency " << round(nearest->first*100*1e-6)/100.0 << " MHz.\n";
	FrequencyTable& table = nearest->second;
	
	double
		phDeg = 90.0 - az*180.0/M_PI,
		zaDeg = za*180.0/M_PI;
	if(phDeg < 0.0) phDeg += 360.0;
	if(!std::isfinite(phDeg) || !std::isfinite(zaDeg))
	{
		for(size_t i=0; i!=4; ++i)
			jonesMatrix[i] = std::numeric_limits<double>::quiet_NaN();
	}
	else {
		alglib::real_1d_array jonesMatrixArray;
		alglib::spline2dcalcv(table.spline, zaDeg, phDeg, jonesMatrixArray);
		for(size_t i=0; i!=4; ++i)
		{
			jonesMatrix[i] = d2c(jonesMatrixArray.getcontent()[i*2]) / table.norm[i];
		}
	}
}

#endif // HAVE_ALGLIB
