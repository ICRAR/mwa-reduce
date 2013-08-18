#include <iostream>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MPosition.h>

#include "sourcesdfwithsamples.h"
#include "radeccoord.h"
#include "nlplfitter.h"
#include "matrix2x2.h"
#include "tilebeam.h"

void testSourceSDFWithSamples()
{
	SourceSDFWithSamples<long double> sdf1;
	sdf1.AddSample(1.0, 2.0);
	sdf1.AddSample(1.0, 3.0);
	std::cout << "SDF1@1 : expecting 1, is " << sdf1.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF1@2 : expecting 1, is " << sdf1.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF1@3 : expecting 1, is " << sdf1.FluxAtFrequency(3.0) << '\n';
	std::cout << "SDF1@4 : expecting 1, is " << sdf1.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF1@2.5 : expecting 1, is " << sdf1.FluxAtFrequency(2.5) << '\n';
	
	std::cout << "SDF1 int[2,3] : expecting 1, is " << sdf1.IntegratedFlux(2.0, 3.0) << '\n';
	std::cout << "SDF1 int[1,3] : expecting 1, is " << sdf1.IntegratedFlux(1.0, 3.0) << '\n';
	std::cout << "SDF1 int[1,2] : expecting 1, is " << sdf1.IntegratedFlux(1.0, 2.0) << '\n';
	std::cout << "SDF1 int[3,4] : expecting 1, is " << sdf1.IntegratedFlux(3.0, 4.0) << '\n';
	std::cout << "SDF1 int[2,4] : expecting 1, is " << sdf1.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF1 int[1,4] : expecting 1, is " << sdf1.IntegratedFlux(1.0, 4.0) << '\n';
	std::cout << "SDF1 int[1,1] : expecting 1, is " << sdf1.IntegratedFlux(1.0, 1.0) << '\n';
	
	SourceSDFWithSamples<long double> sdf2;
	sdf2.AddSample(1.0, 2.0);
	sdf2.AddSample(2.0, 4.0);
	std::cout << "SDF2@1 : expecting 0.5, is " << sdf2.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF2@2 : expecting 1, is " << sdf2.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF2@4 : expecting 2, is " << sdf2.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF2@8 : expecting 4, is " << sdf2.FluxAtFrequency(8.0) << '\n';
	std::cout << "SDF2@3 : expecting 1.5, is " << sdf2.FluxAtFrequency(3.0) << '\n';
	
	std::cout << "SDF2 int[2,4] : expecting 1.5, is " << sdf2.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF2 int[1,2] : expecting 0.75, is " << sdf2.IntegratedFlux(1.0, 2.0) << '\n';
	std::cout << "SDF2 int[3,4] : expecting 1.75, is " << sdf2.IntegratedFlux(3.0, 4.0) << '\n';
	std::cout << "SDF2 int[1,4] : expecting 1.25, is " << sdf2.IntegratedFlux(1.0, 4.0) << '\n';
	std::cout << "SDF2 int[1,1] : expecting 0.5, is " << sdf2.IntegratedFlux(1.0, 1.0) << '\n';
	std::cout << "SDF2 int[3,5] : expecting 2, is " << sdf2.IntegratedFlux(3.0, 5.0) << '\n';
	std::cout << "SDF2 int[1,8] : expecting 2.25, is " << sdf2.IntegratedFlux(1.0, 8.0) << '\n';
	
	SourceSDFWithSamples<long double> sdf3(sdf2);
	sdf3.AddSample(2.0, 6.0);
	std::cout << "SDF3@1 : expecting 0.5, is " << sdf3.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF3@2 : expecting 1, is " << sdf3.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF3@4 : expecting 2, is " << sdf3.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF3@6 : expecting 2, is " << sdf3.FluxAtFrequency(6.0) << '\n';
	std::cout << "SDF3@8 : expecting 2, is " << sdf3.FluxAtFrequency(8.0) << '\n';
	std::cout << "SDF3@3 : expecting 1.5, is " << sdf3.FluxAtFrequency(3.0) << '\n';

	std::cout << "SDF3 int[2,4] : expecting 1.5, is " << sdf3.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF3 int[1,2] : expecting 0.75, is " << sdf3.IntegratedFlux(1.0, 2.0) << '\n';
	std::cout << "SDF3 int[3,4] : expecting 1.75, is " << sdf3.IntegratedFlux(3.0, 4.0) << '\n';
	std::cout << "SDF3 int[1,4] : expecting 1.25, is " << sdf3.IntegratedFlux(1.0, 4.0) << '\n';
	std::cout << "SDF3 int[1,1] : expecting 0.5, is " << sdf3.IntegratedFlux(1.0, 1.0) << '\n';
	std::cout << "SDF3 int[3,5] : expecting 1.875, is " << sdf3.IntegratedFlux(3.0, 5.0) << '\n';
	std::cout << "SDF3 int[1,8] : expecting 1.679, is " << sdf3.IntegratedFlux(1.0, 8.0) << '\n';
	
	SourceSDFWithSamples<long double> sdf4;
	sdf4.AddSample(-1.0, 2.0);
	sdf4.AddSample(-2.0, 4.0);
	sdf4.AddSample(-2.0, 6.0);
	std::cout << "SDF4@1 : expecting -0.5, is " << sdf4.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF4@2 : expecting -1, is " << sdf4.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF4@4 : expecting -2, is " << sdf4.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF4@6 : expecting -2, is " << sdf4.FluxAtFrequency(6.0) << '\n';
	std::cout << "SDF4@8 : expecting -2, is " << sdf4.FluxAtFrequency(8.0) << '\n';
	std::cout << "SDF4@3 : expecting -1.5, is " << sdf4.FluxAtFrequency(3.0) << '\n';
	
	std::cout << "SDF4 int[2,4] : expecting -1.5, is " << sdf4.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF4 int[1,2] : expecting -0.75, is " << sdf4.IntegratedFlux(1.0, 2.0) << '\n';
	std::cout << "SDF4 int[3,4] : expecting -1.75, is " << sdf4.IntegratedFlux(3.0, 4.0) << '\n';
	std::cout << "SDF4 int[1,4] : expecting -1.25, is " << sdf4.IntegratedFlux(1.0, 4.0) << '\n';
	std::cout << "SDF4 int[1,1] : expecting -0.5, is " << sdf4.IntegratedFlux(1.0, 1.0) << '\n';
	std::cout << "SDF4 int[3,5] : expecting -1.875, is " << sdf4.IntegratedFlux(3.0, 5.0) << '\n';
	std::cout << "SDF4 int[1,8] : expecting -1.679, is " << sdf4.IntegratedFlux(1.0, 8.0) << '\n';
	
	SourceSDFWithSamples<long double> sdf5;
	sdf5.AddSample(2.0, 2.0);
	sdf5.AddSample(-2.0, 4.0);
	sdf5.AddSample(-2.0, 6.0);
	sdf5.AddSample(2.0, 8.0);
	std::cout << "SDF5@1 : expecting 4, is " << sdf5.FluxAtFrequency(1.0) << '\n';
	std::cout << "SDF5@2 : expecting 2, is " << sdf5.FluxAtFrequency(2.0) << '\n';
	std::cout << "SDF5@3 : expecting 0, is " << sdf5.FluxAtFrequency(3.0) << '\n';
	std::cout << "SDF5@4 : expecting -2, is " << sdf5.FluxAtFrequency(4.0) << '\n';
	std::cout << "SDF5@6 : expecting -2, is " << sdf5.FluxAtFrequency(6.0) << '\n';
	std::cout << "SDF5@7 : expecting 0, is " << sdf5.FluxAtFrequency(7.0) << '\n';
	std::cout << "SDF5@8 : expecting 2, is " << sdf5.FluxAtFrequency(8.0) << '\n';
	std::cout << "SDF5@9 : expecting 4, is " << sdf5.FluxAtFrequency(9.0) << '\n';
	
	std::cout << "SDF5 int[2,4] : expecting 0, is " << sdf5.IntegratedFlux(2.0, 4.0) << '\n';
	std::cout << "SDF5 int[0,2] : expecting 4, is " << sdf5.IntegratedFlux(0.0, 2.0) << '\n';
	std::cout << "SDF5 int[0,4] : expecting 2, is " << sdf5.IntegratedFlux(0.0, 4.0) << '\n';
}

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
	casa::MEpoch time = casa::MEpoch(casa::MVEpoch(casa::Quantity(4.88193e+09, "s")));
	std::cout << "time=" << time << '\n';
	casa::MPosition arrayPos = casa::MPosition(casa::MVPosition(-2.55952e+06, 5.09585e+06, -2.84899e+06)); // pos of tile 011
	std::cout << "Pos=" << arrayPos << '\n';
	beam.AnalyticJones(time, arrayPos, ra, dec, frequency, gains);
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

int main(int argc, char *argv[])
{
	testBaselineindex();
	testRaDecCoord();
	testSourceSDFWithSamples();
	testNLPLFitter();
	testRotationAngle();
	testBeam();
}
