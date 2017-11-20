#include <iostream>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>

#include "units/radeccoord.h"

#include "nlplfitter.h"
#include "matrix2x2.h"

#include "beam/tilebeam.h"

#include "units/imagecoordinates.h"

std::string raStr(const std::string &str)
{
	std::ostringstream out;
	out << RaDecCoord::ParseRA(str) << " = " << RaDecCoord::RAToString(RaDecCoord::ParseRA(str));
	return out.str();
}

void testRaDecCoord()
{
	std::cout.precision(12);
	std::cout << "ParseRA(\"00h00m00.0s\") = " << raStr("00h00m00.0s") << '\n';
	std::cout << "ParseRA(\"12h00m00.0s\") = " << raStr("12h00m00.0s") << '\n';
	std::cout << "ParseRA(\"00:00:00.0\") = " << raStr("00:00:00.0") << '\n';
	std::cout << "ParseRA(\"-00:00:00.0\") = " << raStr("-00:00:00.0") << '\n';
	std::cout << "ParseRA(\"-00:01:00.0\") = " << raStr("-00:01:00.0") << '\n';
	std::cout << "ParseRA(\"-23:59:59.9\") = " << raStr("-23:59:59.9") << '\n';
	std::cout << "ParseRA(\"11:59:59.9\") = " << raStr("11:59:59.9") << '\n';
	std::cout << "ParseRA(\"-11:59:59.9\") = " << raStr("-11:59:59.9") << '\n';
	std::cout << "ParseRA(\"23:59:59.9\") = " << raStr("23:59:59.9") << '\n';
	std::cout << "ParseDec(\"00d00m00.0s\")=" << RaDecCoord::ParseDec("00d00m00.0s") << '\n';
	std::cout << "ParseDec(\"90d00m00.0s\")=" << RaDecCoord::ParseDec("90d00m00.0s") << '\n';
	std::cout << "ParseDec(\"00.00.00.0\")=" << RaDecCoord::ParseDec("00.00.00.0") << '\n';
	std::cout << "ParseDec(\"89.59.59.9\")=" << RaDecCoord::ParseDec("89.59.59.9") << '\n';
	std::cout << "ParseDec(\"-89.59.59.9\")=" << RaDecCoord::ParseDec("-89.59.59.9") << '\n';
}

int baselineindex(size_t a1, size_t a2, size_t n)
{
	return (a1*(2*n - a1 - 3) + 2*a2 - 2)/2;
}

void testBaselineindex()
{
	size_t n = 5;
	for(size_t a1=0; a1!=n; ++a1)
	{
		for(size_t a2=a1+1; a2!=n; ++a2)
		{
			std::cout << a1 << '\t' << a2 << '\t' << baselineindex(a1, a2, n) << '/' << (n * (n-1) / 2) << '\n';
		}
	}
}

void testNLPLFitter()
{
	double e, f;
	NonLinearPowerLawFitter fitter1;
	fitter1.AddDataPoint(1.0, 1.0);
	fitter1.AddDataPoint(4.0, 1.0);
	fitter1.Fit(e, f);
	std::cout << "1 = " << f << " * x^" << e << '\n';
	
	NonLinearPowerLawFitter fitter2;
	fitter2.AddDataPoint(0.5, -1.0);
	fitter2.AddDataPoint(0.5, 3.0);
	fitter2.AddDataPoint(2.0, 0.0);
	fitter2.AddDataPoint(2.0, 8.0);
	
	fitter2.Fit(e, f);
	std::cout << "2x = " << f << " * x^" << e << '\n';
	
	NonLinearPowerLawFitter fitter3;
	fitter3.AddDataPoint(1.0, -1.0);
	fitter3.AddDataPoint(1.0, 3.0);
	fitter3.AddDataPoint(4.0, 0.0);
	fitter3.AddDataPoint(4.0, 1.0);
	
	fitter3.Fit(e, f);
	std::cout << "x^-0.5 = " << f << " * x^" << e << '\n';

	NonLinearPowerLawFitter fitter4;
	fitter4.AddDataPoint(1000000.0, -1.0);
	fitter4.AddDataPoint(1000000.0, 3.0);
	fitter4.AddDataPoint(4000000.0, 0.0);
	fitter4.AddDataPoint(4000000.0, 1.0);
	
	fitter4.Fit(e, f);
	std::cout << "1000 x^-0.5 = " << f << " * x^" << e << '\n';
	
	NonLinearPowerLawFitter fitter5;
	fitter5.AddDataPoint(1000000.0, -3.0);
	fitter5.AddDataPoint(1000000.0, 1.0);
	fitter5.AddDataPoint(4000000.0, 0.0);
	fitter5.AddDataPoint(4000000.0, -1.0);
	
	fitter5.Fit(e, f);
	std::cout << "-1000 x^-0.5 = " << f << " * x^" << e << '\n';
}

std::string matrixToStr(std::complex<double> *m)
{
	std::ostringstream str;
	str
		<< '[' << m[0].real() << ", " << m[1].real()
		<< "; " << m[2].real() << ", " << m[3].real() << ']';
	return str.str();
}

void testRotationAngle()
{
	std::complex<double> r[4];
	r[0] = 1.0; r[1] = 0.0;
	r[2] = 0.0; r[3] = 1.0;
	std::cout << "angle(" << matrixToStr(r) << ") = " << Matrix2x2::RotationAngle(r) << '\n';
	r[0] = 1.0; r[1] = 0.0;
	r[2] = 0.0; r[3] = -1.0;
	std::cout << "angle(" << matrixToStr(r) << ") = " << Matrix2x2::RotationAngle(r) << '\n';
	r[0] = -1.0; r[1] = 0.0;
	r[2] = 0.0; r[3] = -1.0;
	std::cout << "angle(" << matrixToStr(r) << ") = " << Matrix2x2::RotationAngle(r) << '\n';
	r[0] = 0.5; r[1] = -0.5;
	r[2] = 0.5; r[3] = 0.5;
	std::cout << "angle(" << matrixToStr(r) << ") = " << Matrix2x2::RotationAngle(r) << '\n';
	r[0] = 1.0; r[1] = -1.0;
	r[2] = 1.0; r[3] = 1.0;
	std::cout << "angle(" << matrixToStr(r) << ") = " << Matrix2x2::RotationAngle(r) << '\n';
	r[0] = 0.0; r[1] = 1.0;
	r[2] = -1.0; r[3] = 0.0;
	std::cout << "angle(" << matrixToStr(r) << ") = " << Matrix2x2::RotationAngle(r) << '\n';
	r[0] = 0.5; r[1] = -0.5;
	r[2] = 0.0; r[3] = 1.0;
	std::cout << "angle(" << matrixToStr(r) << ") = " << Matrix2x2::RotationAngle(r) << '\n';
}

void showBeam(const double *delays, double ra, double dec)
{
	TileBeam beam(delays);
	std::complex<double> gains[4];
	
	dec = dec *(M_PI/180.0);
	ra = ra *(M_PI/180.0);
	double frequency = 150000000.0;
	casacore::MEpoch time = casacore::MEpoch(casacore::MVEpoch(casacore::Quantity(4.88193e+09, "s")));
	std::cout << "time=" << time << '\n';
	casacore::MPosition arrayPos = casacore::MPosition(casacore::MVPosition(-2.55952e+06, 5.09585e+06, -2.84899e+06)); // pos of tile 011
	std::cout << "Pos=" << arrayPos << '\n';
	beam.ArrayResponse(time, arrayPos, ra, dec, frequency, gains);
	std::cout << "Gains: " << matrixToStr(gains) << '\n';
}

void testBeam()
{
	double delays1[16] = {
		0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0 };
	showBeam(delays1, 0.0, -30.0);

	double delays2[16] = {
		0.0, 1.0, 2.0, 3.0,
		0.0, 1.0, 2.0, 3.0,
		0.0, 1.0, 2.0, 3.0,
		0.0, 1.0, 2.0, 3.0 };
	showBeam(delays2, 0.0, -30.0);
	
	double delays3[16] = {
		0.0, 0.0, 0.0, 0.0,
		1.0, 1.0, 1.0, 1.0,
		2.0, 2.0, 2.0, 2.0,
		3.0, 3.0, 3.0, 3.0 };
	showBeam(delays3, 0.0, -30.0);
	
	showBeam(delays1, 0.0, -20.0);

	showBeam(delays1, 180.0, -20.0);

}

void testCentreRA()
{
	std::vector<double> ras(5);
	for(size_t i=0; i!=5; ++i)
		ras[i] = i;
	std::cout << "Mean RA 0,1,2,3,4: " << ImageCoordinates::MeanRA(ras) << '\n';
	
	for(size_t i=0; i!=5; ++i)
		ras[i] = double(i) + 3.0;
	std::cout << "Mean RA 3,4,5,6,7: " << ImageCoordinates::MeanRA(ras) << '\n';
}

int main(int argc, char *argv[])
{
	testCentreRA();
	testBaselineindex();
	testRaDecCoord();
	testNLPLFitter();
	testRotationAngle();
	testBeam();
}
