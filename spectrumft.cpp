#include "spectrumft.h"
#include "windowfunction.h"

#include "model/measuredsed.h"

#include "universe.h"
#include "gnuplot.h"

#include <cmath>

SpectrumFT::SpectrumFT(size_t channelCount) : _channelCount(channelCount), _outputFit(false)
{
	_inputData = reinterpret_cast<std::complex<double>*>(fftw_malloc(channelCount * sizeof(std::complex<double>)));
	_outputData = reinterpret_cast<std::complex<double>*>(fftw_malloc(channelCount * sizeof(std::complex<double>)));
	_fftPlan = fftw_plan_dft_1d(channelCount, reinterpret_cast<fftw_complex*>(_inputData), reinterpret_cast<fftw_complex*>(_outputData), FFTW_FORWARD, FFTW_ESTIMATE);
}

SpectrumFT::~SpectrumFT()
{
	fftw_free(_outputData);
	fftw_free(_inputData);
}

void SpectrumFT::GetFTPower(ao::uvector<double>& destination, const MeasuredSED& sed, size_t nTermsFitted, bool toTemperature, bool useLomb) const
{
	if(sed.MeasurementCount() != _channelCount)
		throw std::runtime_error("Number of measurements in SED != specified channel count");
	
	ao::uvector<double> realData(_channelCount);
	ao::uvector<double>::iterator dataPtr = realData.begin();
	
	ao::uvector<double> terms;
	if(nTermsFitted!=0)
		sed.FitLogPolynomial(terms, nTermsFitted, Polarization::StokesI);
	
	size_t sampleCount = 0;
	double missing = useLomb ? std::numeric_limits<double>::quiet_NaN() : 0.0;
	for(MeasuredSED::const_iterator i=sed.begin(); i!=sed.end(); ++i)
	{
		const Measurement& m = i->second;
		double flux = m.FluxDensityFromIndex(0);
		double fitVal = nTermsFitted!=0 ?  (NonLinearPowerLawFitter::Evaluate(m.FrequencyHz(), terms)) : 0.0;
		if(_outputFit && nTermsFitted>0)
		{
			if(std::isfinite(flux))
			{
				*dataPtr = fitVal;
				++sampleCount;
			}
			else
				*dataPtr = missing;
		}
		else {
			if(std::isfinite(flux))
			{
				*dataPtr = flux-fitVal;
				++sampleCount;
			}
			else
				*dataPtr = missing;
		}
		++dataPtr;
	}
	
	if(useLomb)
	{
		interpolate(realData);
		twoPoint(destination, realData);
		//lombPeriodogram(destination, realData);
	}
	else {
		ao::uvector<std::complex<double>> data;
		applyWindow(realData);
		data = ao::uvector<std::complex<double>>(realData.begin(), realData.end());
		transformSpectrum(data, sampleCount);
		destination.resize(_channelCount);
		for(size_t i=0; i!=_channelCount; ++i)
		{
			destination[i] = data[i].real()*data[i].real() + data[i].imag()*data[i].imag();
		}
	}
	if(toTemperature)
	{
		ConvertFluxToTemperature(destination, _frequencies.data(), _pBeamOmegas.data());
		//double kParResolution = GetMaxKParallell(sed, false) / sed.MeasurementCount();
		// TODO
	}
}

double SpectrumFT::GetAveragePowerInKRange(double kLow, double kHigh, const MeasuredSED& sed, size_t nTermsFitted, bool toTemperature, bool useLomb)
{
	double kMax = GetMaxKParallell(sed, false);
	ao::uvector<double> values;
	GetFTPower(values, sed, nTermsFitted, toTemperature, useLomb);
	double powerSum = 0.0;
	size_t n = 0;
	for(size_t i=0; i!=values.size(); ++i)
	{
		double k = kMax*double(i)/values.size();
		if(k >= kLow && k < kHigh)
		{
			powerSum += values[i];
			n++;
		}
	}
	return powerSum / n;
}

void SpectrumFT::applyWindow(ao::uvector<double>& data) const
{
	double windowSum = 0.0;
	for(size_t i=0; i!=data.size(); ++i)
	{
		double windowValue =
			WindowFunction::EvaluateBlackmanNutall(data.size(), i);
		windowSum += windowValue;
		data[i] *= windowValue;
	}
	double factor = data.size() / windowSum;
	for(size_t i=0; i!=data.size(); ++i)
	{
		data[i] *= factor;
	}
}

void SpectrumFT::transformSpectrum(ao::uvector<std::complex<double>>& data, size_t sampleCount) const
{
	const double factor = 1.0 / sqrt(sampleCount);
	
	for(size_t i=0; i!=_channelCount; ++i)
		_inputData[i] = data[i];
	fftw_execute_dft(_fftPlan, reinterpret_cast<fftw_complex*>(_inputData), reinterpret_cast<fftw_complex*>(_outputData));
	for(size_t i=0; i!=_channelCount; ++i)
		data[i] = _outputData[i] * factor;
}

void SpectrumFT::interpolate(ao::uvector<double>& data) const
{
	size_t i=0;
	while(i!=data.size() && !std::isfinite(data[i]))
		++i;
	
  while(i!=data.size())
  {
    if(!std::isfinite(data[i]))
      interpolate(data, i);
		++i;
  }
}

void SpectrumFT::interpolate(ao::uvector<double>& data, size_t i) const
{
	size_t n = data.size();
	size_t p = i;
	while(p<n && !std::isfinite(data[p]))
	{
		++p;
	}
	if(p >= n)
		return;
	double startVal = data[i-1], endVal = data[p];
	double l = double(p) - double(i);
	for(size_t x=i; x!=p; ++x)
	{
		data[x] = startVal + (endVal-startVal) * (double(x)-double(i))/l;
		//std::cout << x << ": " << data[x] << '\n';
	}	
}

void SpectrumFT::twoPoint(ao::uvector<double>& output, const ao::uvector<double>& input) const
{
	size_t n = 0;
	double sum = 0.0;
	for(size_t i=0; i!=input.size(); ++i)
	{
		if(std::isfinite(input[i]))
		{
			sum += input[i];
			n++;
		}
	}
	
	ao::uvector<double> weightSums(input.size(), 0.0);
	output.assign(input.size(), 0.0);
	const double mu = sum/n;
	
	double sumSq = 0.0;
	for(size_t i=0; i!=input.size(); ++i)
	{
		if(std::isfinite(input[i]))
			sumSq += (input[i]-mu) * (input[i]-mu);
	}
	const double sigmaSq = sumSq/n;
	
	for(size_t i=0; i!=input.size(); ++i)
	{
		for(size_t j=0; j!=input.size(); ++j)
		{
			double v1 = input[i]-mu, v2 = input[j]-mu;
			if(std::isfinite(v1) && std::isfinite(v2))
			{
				size_t index = std::abs(int(j)-int(i));
				output[index] += v1 * v2;
				weightSums[index] += 1;
			}
		}
	}
	
	size_t count = 0;
	for(size_t i=0; i!=input.size(); ++i)
	{
		if(weightSums[i] == 0.0)
			output[i] = 0.0;
		else {
			output[i] = output[i]/(weightSums[i]*sigmaSq);
			++count;
		}
	}
	applyWindow(output);
	ao::uvector<std::complex<double>> data(input.size());
	for(size_t i=0; i!=input.size(); ++i)
		data[i] = output[i];
	transformSpectrum(data, count);
	for(size_t i=0; i!=input.size(); ++i)
		output[i] = std::abs(data[i]);
}

void SpectrumFT::lombPeriodogram(ao::uvector<double>& output, const ao::uvector<double>& input) const
{
	ao::uvector<double> weights(input.size(), 1.0);
	applyWindow(weights);
	
	double n = 0;
	double mu = 0.0;
	for(size_t i=0; i!=input.size(); ++i)
	{
		if(std::isfinite(input[i]))
		{
			mu += input[i] * weights[i];
			n += weights[i];
		}
		else weights[i] = 0.0;
	}
	
	if(n == 0 || n==1)
	{
		output.assign(input.size(), 0.0);
		return;
	}
	mu /= n;
	
	
	ao::uvector<double> scaledInput(input.size());
	for(size_t i=0; i!=input.size(); ++i)
		scaledInput[i] = input[i] - mu;
	
	double sigmaSq = 0.0;
	for(size_t i=0; i!=scaledInput.size(); ++i)
	{
		if(std::isfinite(scaledInput[i]))
		{
			sigmaSq += scaledInput[i]*scaledInput[i]*weights[i];
		}
	}
	sigmaSq = sigmaSq / (double(n) - 1.0);
	double sigma = sqrt(sigmaSq);
	
	//for(size_t i=0; i!=input.size(); ++i)
	//	scaledInput[i] = (input[i] - mu) / sigma;
	
	std::cout << "Average signal: " << mu << " Jy, sigma=" << sigma << "\n";
	
	output.resize(scaledInput.size());
	for(size_t i=0; i!=output.size(); ++i)
	{
		double omega = 2.0 * M_PI * double(i) / double(output.size());
		
		// Find tau
		double tau;
		if(i*2 == output.size())
		{
			tau = 0.0;
		}
		else {
			double tauNum = 0.0, tauDen = 0.0;
			for(size_t j=0; j!=scaledInput.size(); ++j)
			{
				if(std::isfinite(scaledInput[j]))
				{
					double s, c;
					sincos(2.0 * omega * j, &s, &c);
					tauNum += s * weights[j];
					tauDen += c * weights[j];
				}
			}
			tau = atan2(tauNum, tauDen) / (2.0 * omega);
		}
		
		// evaluate sums of P
		double
			sumCos = 0.0, sumSin = 0.0,
			sumSqCos = 0.0, sumSqSin = 0.0;
		for(size_t j=0; j!=scaledInput.size(); ++j)
		{
			double h = scaledInput[j];
			double w = weights[j];
			if(std::isfinite(h))
			{
				double s, c;
				sincos(omega * (j - tau), &s, &c);
				sumCos += h * c * w;
				sumSin += h * s * w;
				sumSqCos += c * c * w;
				sumSqSin += s * s * w;
			}
		}
		if(i*2 == output.size())
			output[i] = (sumCos*sumCos / sumSqCos);// / (2.0 /* * sigmaSq*/);
		else
			output[i] = ((sumCos*sumCos / sumSqCos) + (sumSin*sumSin / sumSqSin));// / (2.0 /* * sigmaSq*/);
	}
}

double SpectrumFT::getKForFreqStep(double freqStepHz, double centralFreq)
{
	double redshiftA = Universe::HIRedshift(centralFreq-freqStepHz/2.0);
	double redshiftB = Universe::HIRedshift(centralFreq+freqStepHz/2.0);
	double distA = Universe::ComovingDistanceInMPCoverH(redshiftA);
	double distB = Universe::ComovingDistanceInMPCoverH(redshiftB);
	return 2.0*M_PI/(distA-distB);
}

double SpectrumFT::FrequencyToDistance(double frequency)
{
	double redshift = Universe::HIRedshift(frequency);
	return Universe::ComovingDistanceInMPCoverH(redshift);
}

double SpectrumFT::GetMaxKParallell(const MeasuredSED& sed, bool verbose)
{
	double freqStep = (sed.HighestFrequency()-sed.LowestFrequency())/sed.MeasurementCount();
	double redshiftLo = Universe::HIRedshift(sed.HighestFrequency());
	double redshiftLoNext = Universe::HIRedshift(sed.HighestFrequency() - freqStep);
	double redshiftHi = Universe::HIRedshift(sed.LowestFrequency());
	double redshiftHiNext = Universe::HIRedshift(sed.LowestFrequency() + freqStep);
	double distLo = Universe::ComovingDistanceInMPCoverH(redshiftLo);
	double distLoNext = Universe::ComovingDistanceInMPCoverH(redshiftLoNext);
	double distHi = Universe::ComovingDistanceInMPCoverH(redshiftHi);
	double distHiNext = Universe::ComovingDistanceInMPCoverH(redshiftHiNext);
	//double smallestDistance = 1.0/(distHi-distLo);
	double largestDistance = 2.0*M_PI*sed.MeasurementCount()/(distHi-distLo);
	double timeRange = sed.MeasurementCount()/(sed.HighestFrequency()-sed.LowestFrequency());
	if(verbose)
	{
		std::cout << "- Frequency range: " << round(sed.LowestFrequency()*1e-6) << " - " << round(sed.HighestFrequency()*1e-6) << " MHz\n";
		std::cout << "- Time range: " << round(timeRange*1.0e11)*1e-2 << " ns\n";
		std::cout << "- Redshift range: " << redshiftLo << "," << redshiftLoNext << ",...," << redshiftHiNext << ',' << redshiftHi << '\n';
		std::cout << "- Distance: " << distLo << ',' << distLoNext << ",...," << distHiNext << ',' << distHi << " Mpc/h, or " << Universe::ComovingDistanceInMPC(redshiftLo) << " - " << Universe::ComovingDistanceInMPC(redshiftHi) << " Mpc\n";
		std::cout << "- Largest probed distance: " << largestDistance << " h/Mpc (with " << sed.MeasurementCount() << " steps)\n";
	}
	return largestDistance;
}

double SpectrumFT::GetMaxDelayInMuSec(const MeasuredSED& sed)
{
	return sed.MeasurementCount()/(sed.HighestFrequency()-sed.LowestFrequency()) * 1e6;
}

void SpectrumFT::PreparePSPlot(GNUPlot& plot, double centralFrequency, bool useDelay)
{
	double lengths[] = { 0.74, // in micro s
		1.2,
		1.9,
		2.6,
		3.3,
		4.3
	};
	std::cout << "- Cable reflections:\n";
	if(useDelay)
	{
		for(size_t i=1; i!=6; ++i)
		{
			plot.VerticalLine(lengths[i], 2);
		}
		for(size_t i=1; i!=12; ++i)
		{
			plot.VerticalLine(lengths[0]*i, 1);
			plot.VerticalLine(1e6/1280000.0*i, 3);
		}
	}
	else {
		double k0 = getKForFreqStep(1000000.0/lengths[0], centralFrequency);
		double kSB = getKForFreqStep(1000000.0/(1e6/1280000.0), centralFrequency);
		std::cout << lengths[0] << " -> k = " << k0 << '\n';
		for(size_t i=1; i!=12; ++i)
		{
			plot.VerticalLine(k0*i, 1);
			plot.VerticalLine(kSB*i, 1);
		}
		
		for(size_t i=1; i!=6; ++i)
		{
			double k = getKForFreqStep(1000000.0/lengths[i], centralFrequency);
			plot.VerticalLine(k, 2);
			std::cout << lengths[i] << " -> k = " << k << '\n';
		}
	}
}
