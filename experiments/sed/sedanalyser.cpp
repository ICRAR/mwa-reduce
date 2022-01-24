#include "sedanalyser.h"

#include "../../gnuplot.h"
#include "../../spectrumft.h"
#include "../../gnustatplot.h"
#include "../../logbinnedplot.h"

void SEDAnalyser::Process()
{
	const Model& model = *_model;
	long double centreFrequency = CentralFrequency();
	
	std::vector<long double> ras;
	long double meanDec = 0.0;
	for(size_t i=0; i!=model.SourceCount(); ++i)
	{
		const ModelSource& s = model.Source(i);
		ras.push_back(s.MeanRA());
		meanDec += s.MeanDec();
	}
	long double meanRA = aocommon::ImageCoordinates::MeanRA(ras);
	meanDec /= model.SourceCount();
	std::cout << "Catalogue centre: " << RaDecCoord::RAToString(meanRA) << " " << RaDecCoord::DecToString(meanDec) << '\n';
	
	if(_edgeChannels != 0)
		removeEdgeChannels();
	if(_flaggedEndChannels != 0)
		removeEndChannels();
	removeBadChannels();
	
	ProgressBar progress("Analysing SEDs");
	for(size_t i=0; i!=model.SourceCount(); ++i)
	{
		const ModelSource& source = model.Source(i);
		SourceInfo newSource;
		const MeasuredSED sed = source.GetIntegratedMSED();
		sed.FitPowerlaw(newSource.plFactor, newSource.plExponent, aocommon::Polarization::StokesI);
		newSource.index = i;
		newSource.modalFlux = newSource.plFactor * powl(centreFrequency, newSource.plExponent);
		newSource.componentCount = source.ComponentCount();
		newSource.rms = sourceRMS(sed, newSource.plFactor, newSource.plExponent);
		newSource.qrms = sourceRMS(sed, aocommon::Polarization::StokesQ);
		newSource.urms = sourceRMS(sed, aocommon::Polarization::StokesU);
		newSource.vrms = sourceRMS(sed, aocommon::Polarization::StokesV);
		newSource.diffRMS = diffRMS(sed);
		if(!std::isfinite(newSource.diffRMS))
			newSource.diffRMS = newSource.rms;
		newSource.maxStep = maxStep(sed);
		newSource.stokesVFrac = std::abs(sed.AverageFlux(aocommon::Polarization::StokesV) / sed.AverageFlux(aocommon::Polarization::StokesI));
		newSource.distance = aocommon::ImageCoordinates::AngularDistance<long double>(meanRA, meanDec, source.MeanRA(), source.MeanDec());
		
		sed.FitLogPolynomial(newSource.terms, _nTermsInPS, aocommon::Polarization::StokesI, _referenceFrequency);
		newSource.pl2ndOrder = newSource.terms[2];
		newSource.rms2ndOrder = sourceRMS(sed, newSource.terms);
		
		newSource.maxError = 0.0;
		newSource.maxAbsError = 0.0;
		newSource.minAbsError = std::numeric_limits<long double>::max();
		size_t channel = 0;
		for(MeasuredSED::const_iterator m=sed.begin(); m!=sed.end(); ++m)
		{
			// There seems to be some consistent high values in the last two channels... Systematics?
			// Ignore for now
			if(m->first > 197.61*1e6)
				break;
			double flux = m->second.FluxDensity(aocommon::Polarization::StokesI);
			double freq = m->first;
			double model = NonLinearPowerLawFitter::Evaluate(freq, newSource.terms, _referenceFrequency);
			double absErr = (flux - model);
			double err = absErr / newSource.rms;
			//bool edgeChannel =
			//	channel%32<=2 || channel%32>=29;
			if(std::abs(err) > std::abs(newSource.maxError))// && !edgeChannel)
			{
				newSource.maxError = err;
				newSource.maxErrorFrequency = m->first;
			}
			if(absErr > newSource.maxAbsError)
			{
				newSource.maxAbsError = absErr;
			}
			if(absErr < newSource.minAbsError)
			{
				newSource.minAbsError = absErr;
			}
			++channel;
		}
		//long double aTemp, bTemp;
		//sed.FitPowerlaw2ndOrder(aTemp, bTemp, newSource.pl2ndOrder, Polarization::StokesI);
		//newSource.rms2ndOrder = sourceRMS(sed, aTemp, bTemp, newSource.pl2ndOrder);
		//std::cout << "First fit=" << newSource.plExponent << ", 2nd=" << aTemp << ", c=" << newSource.pl2ndOrder << '\n';
		//std::cout << newSource.terms[0] << " (" << log(newSource.plFactor) - newSource.plExponent*log(1e-8)<< ")";
		//for(size_t i=1; i!=newSource.terms.size(); ++i)
		//	std::cout << ',' << newSource.terms[i];
		//std::cout << '\n';
		//double multiTermRMS = sourceRMS(sed, newSource.terms);
		//newSource.terms.assign(newSource.terms.size(), 0.0);
		//newSource.terms[0] = log(newSource.plFactor) - newSource.plExponent*log(1e-8);
		//newSource.terms[1] = newSource.plExponent;
		//std::cout << "rms=" << newSource.rms << ", " << newSource.rms2ndOrder << ", " << multiTermRMS << " ," << sourceRMS(sed, newSource.terms) << '\n';
		
		if(_analysePassband)
			newSource.subbandCorrelation = getSubbandShapeCorrelation(sed);
		else
			newSource.subbandCorrelation = 0.0;
		if(source.ComponentCount() == 1) // skip complex sources
			_sources.push_back(newSource);
		
		progress.SetProgress(i+1, model.SourceCount());
	}
}

void SEDAnalyser::ProcessRM()
{
	const Model& model = *_model;
	ProgressBar progress("Doing RM synthesis");
	for(size_t i=0; i!=_sources.size(); ++i)
	{
		SourceInfo& newSource = _sources[i];
		const ModelSource& source = model.Source(newSource.index);
		const MeasuredSED sed = source.GetIntegratedMSED();
		RMSynthesis rmSynthesis(sed);
		rmSynthesis.Synthesize();
		rmSynthesis.PeakWithNonZeroRM(5, newSource.rmPeakPos, newSource.rmPeakValue);
		double tmp;
		rmSynthesis.PeakWithSmallRM(5, tmp, newSource.rmZeroValue);
		progress.SetProgress(i+1, model.SourceCount());
	}
}

void SEDAnalyser::SaveResults()
{
	std::cout << "Sorting... " << std::flush;
	std::sort(_sources.rbegin(), _sources.rend());
	std::cout << "DONE\n";
	
	std::cout << "Name\tcomplex?\tS" << CentralFrequency()*1e-6 << "\tSI\tRMS\n";
	
	double spectralIndexSum = 0.0;
	size_t siAvgCount = 0;
	
	const Model& model = *_model;
	_brightestSourcesModel = Model();
	for(std::vector<SourceInfo>::const_iterator i=_sources.begin(); i!=_sources.end(); ++i)
	{
		const ModelSource& ms = model.Source(i->index);
		std::cout << ms.Name() << '\t' << i->modalFlux << '\t' << i->plExponent << '\t' << i->rms << '\n';
		spectralIndexSum += i->plExponent;
		++siAvgCount;
		if(siAvgCount < 20) _brightestSourcesModel.AddSource(ms);
	}

	std::sort(_sources.rbegin(), _sources.rend());
	_brightestSourcesModel.Save("brightest-sources-model.txt");
	
	std::cout << "\nSorting on RMS... " << std::flush;
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerRMS);
	std::cout << "DONE\n";
	std::cout << "High RMS sources:\n";
	outputList(_sources, *_model, "highrms", false, 20);
	std::cout << "\nLow RMS sources:\n";
	outputList(_sources, *_model, "lowrms", true, 20);
	
	std::cout << "\nSorting on spectral index... " << std::flush;
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerSI);
	std::cout << "DONE\n";
	std::cout << "High SI sources:\n";
	outputList(_sources, *_model, "highsi", false, 40);
	std::cout << "\nLow SI sources:\n";
	outputList(_sources, *_model, "lowsi", true, 10);
	
	std::cout << "\nSorting on Stokes V frac... " << std::flush;
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerVFrac);
	std::cout << "DONE\n";
	std::cout << "High V-fractional sources:\n";
	outputList(_sources, *_model, "highvfrac", false, 20);

	std::cout << "\nSorting on max abs error... " << std::flush;
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerMaxAbsError);
	std::cout << "DONE\n";
	std::cout << "Max abs errors:\n";
	outputList(_sources, *_model, "maxabserror", false, 3);
	
	std::cout << "\nSorting on min abs error... " << std::flush;
	std::sort(_sources.begin(), _sources.end(), &SourceInfo::hasLowerMinAbsError);
	std::cout << "DONE\n";
	std::cout << "Min abs errors:\n";
	outputList(_sources, *_model, "minabserror", false, 3);
	
	std::cout << "\nSorting on max error... " << std::flush;
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerMaxError);
	std::cout << "DONE\n";
	std::cout << "Max errors:\n";
	outputList(_sources, *_model, "maxerror", false, 40);
	
	std::cout << "\nSorting on RM value... " << std::flush;
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerRMValue);
	std::cout << "DONE\n";
	std::cout << "Relatively high RM values:\n";
	outputList(_sources, *_model, "highrm", false, 20);
	
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerSubbandCorrelation);
	std::cout << "High sub-band correlation:\n";
	outputList(_sources, *_model, "highsubbandcorrelation", false, 20);
	std::cout << "Low sub-band correlation:\n";
	outputList(_sources, *_model, "lowsubbandcorrelation", true, 20);
	
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLower2ndOrder);
	std::cout << "High 2nd order value:\n";
	outputList(_sources, *_model, "high2ndorder", false, 20);
	std::cout << "Low 2nd order value:\n";
	outputList(_sources, *_model, "low2ndorder", true, 20);
	
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLower2ndRMSDiff);
	std::cout << "High 2nd order RMS diff:\n";
	outputList(_sources, *_model, "high2ndrmsdiff", false, 20);
	std::cout << "Low 2nd order RMS diff:\n";
	outputList(_sources, *_model, "low2ndrmsdiff", true, 20);
	
	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerDiffRMS);
	std::cout << "High diff RMS:\n";
	outputList(_sources, *_model, "highdiffrms", false, 20);
	std::cout << "Low diff RMS:\n";
	outputList(_sources, *_model, "lowdiffrms", true, 20);

	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerQRMS);
	std::cout << "High Q RMS:\n";
	outputList(_sources, *_model, "highqrms", false, 10);

	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerURMS);
	std::cout << "High U RMS:\n";
	outputList(_sources, *_model, "highurms", false, 10);

	std::sort(_sources.rbegin(), _sources.rend(), &SourceInfo::hasLowerVRMS);
	std::cout << "High V RMS:\n";
	outputList(_sources, *_model, "highvrms", false, 10);
	
	make_RA_vs_SI_plot();
	make_Dec_vs_SI_plot();
	make_SI_vs_brightness_plot();
	OutputSIStatistics();
	OutputSIBrightnessStatistics();
	write_csv_file();
	write_cath_csv_file();
	
	if(_sourceCountAnalysis)
		measureSourceCountNoise();
	
	if(_powerSpectrumAnalysis)
	{
		MakeAveragedPSPlot();
		MakePSPlot();
		MakeNoisePSPlot();
	}
}

void SEDAnalyser::outputList(const std::vector< SEDAnalyser::SourceInfo >& sortedList, const Model& model, const std::string& name, bool reverse, size_t count)
{
	std::ofstream latexTable("latex-"+name+"-table.tex");
	latexTable <<
		"\\begin{tabular}{rllrlr}%\n"
		"\\textbf{\\#}&\\textbf{RA}&\\textbf{Dec}&\\textbf{$S_{168}$}&"
		"\\textbf{$\\sigma_\\textrm{40~kHz}$}&\\textbf{SI}\\\\\n"
		"\\hline%\n";
	Model outputModel;
	for(size_t i=0; i<std::min(size_t(count), sortedList.size()); ++i)
	{
		size_t index = reverse ? (sortedList.size()-i-1) : i;
		const SourceInfo& source = sortedList[index];
		outputModel.AddSource(model.Source(source.index));
		std::cout << source.Description(model) << '\n';
		latexTable << source.LatexEntry(model) << "\\\\\n";
	}
	latexTable << "\\end{tabular}%\n";
	outputModel.Save(name + "-sources-model.txt");
}

void SEDAnalyser::outputGeneralSIStats()
{
	size_t steepCount = 0, flatCount = 0;
	double curvature = 0.0;
	for(const SEDAnalyser::SourceInfo& s : _sources)
	{
		if(s.plExponent < -1.3)
			++steepCount;
		if(s.plExponent > 0.0)
			++flatCount;
		curvature += s.pl2ndOrder;
	}
	std::cout << "#sources with SI < -1.3 : " << steepCount << '\n';
	std::cout << "#sources with SI > 0 : " << flatCount << '\n';

	curvature = curvature / _sources.size();
	std::cout << "Average curvature: " << curvature;
	double cSq = 0.0;
	for(const SEDAnalyser::SourceInfo& s : _sources)
	{
		cSq += (s.pl2ndOrder - curvature) * (s.pl2ndOrder - curvature);
	}
	std::cout << " ± " << sqrt(cSq/_sources.size()) << '\n';
}


void SEDAnalyser::OutputSIBrightnessStatistics()
{
	std::cout << "\nbinStart binEnd N SI_avg SI_median SI_stddev SI_min SI_max lowestFlux\n";
	std::sort(_sources.rbegin(), _sources.rend());
	double brightest = _sources.front().modalFlux;
	double faintest = _sources.back().modalFlux;
	double step = pow(brightest/faintest, 1.0/double(_siBinCount));
	double bin = faintest;
	for(size_t i=0; i!=_siBinCount; ++i)
	{
		double binEnd = bin * step;
		outputSIBrightnessStatistics(bin, binEnd);
		bin = binEnd;
	}
}

void SEDAnalyser::outputSIBrightnessStatistics(double binStart, double binEnd)
{
	double step = binEnd / binStart;
	std::vector<SEDAnalyser::SourceInfo> bin;
	for(const SEDAnalyser::SourceInfo& s : _sources)
	{
		if(s.modalFlux >= binStart && s.modalFlux < binEnd)
		{
			bin.push_back(s);
		}
	}
	if(bin.size()>=2) {
		std::cout << binStart*sqrt(step) << '\t' << binStart << '\t' << binEnd << '\t';
		outputSIStats(bin, bin.size());
	}
}

void SEDAnalyser::outputSIStats(const std::vector<SEDAnalyser::SourceInfo>& sortedList, size_t n)
{
	size_t actualN = n > sortedList.size() ? sortedList.size() : n;
	std::vector<long double> sis(actualN);
	long double sum = 0.0;
	for(size_t i=0; i!=actualN; ++i)
	{
		const SourceInfo& source = sortedList[i];
		sis[i] = source.plExponent;
		sum += source.plExponent;
	}
	long double lowestFlux = sortedList[actualN-1].modalFlux;
	long double average = sum / actualN;
	long double squaredDistSum = 0.0;
	for(size_t i=0; i!=actualN; ++i)
	{
		long double dist = sis[i]-average;
		squaredDistSum += dist*dist;
	}
	long double stddev = sqrtl(squaredDistSum/actualN);
	
	std::sort(sis.begin(), sis.end());
	long double median;
	if(actualN%2 == 0)
		median = (sis[actualN/2-1] + sis[actualN/2]) * 0.5;
	else
		median = sis[actualN/2];
	long double highestSI = sis.back(), lowestSI = sis.front();
	
	std::cout << actualN << ' ' << average << ' ' << median << ' ' << stddev << ' ' << lowestSI << ' ' << highestSI << ' ' << lowestFlux << ' ' << '\n';
	
	long double firstBinValue = roundl(lowestSI*10.0)/10.0;
	std::vector<size_t> histData(10.0*(roundl(highestSI*10.0)/10.0-firstBinValue)+1, 0);
	for(size_t i=0; i!=actualN; ++i)
	{
		long double val = sis[i];
		long double bindex = (roundl(val*10.0)/10.0 - firstBinValue) * 10.0;
		histData[size_t(roundl(bindex))]++;
	}
	
	std::ostringstream histFilename;
	histFilename << "histogram-n" << actualN << ".txt";
	std::ofstream histFile(histFilename.str());
	for(size_t i=0; i!=histData.size(); ++i)
	{
		long double binCentre = firstBinValue + i/10.0;
		histFile << i << '\t' << binCentre << '\t' << histData[i] << '\n';
	}
}

double SEDAnalyser::BestRMS()
{
	std::sort(_sources.begin(), _sources.end(), &SourceInfo::hasLowerDiffRMS);
	const SourceInfo& source = _sources.front();
	return source.diffRMS;
}

double SEDAnalyser::sourceRMS(const MeasuredSED& sed, long double factor, long double exponent)
{
	long double RMSsum = 0.0L;
	size_t count = 0;
	for(MeasuredSED::const_iterator i=sed.begin(); i!=sed.end(); ++i)
	{
		long double flux = i->second.FluxDensity(aocommon::Polarization::StokesI);
		long double powerLawValue = factor * powl(i->second.FrequencyHz(), exponent);
		if(std::isfinite(flux))
		{
			long double diff = flux - powerLawValue;
			++count;
			RMSsum += diff * diff;
		}
	}
	return sqrtl(RMSsum / (long double)(count));
}

double SEDAnalyser::sourceRMS(const MeasuredSED& sed, aocommon::PolarizationEnum pol)
{
	long double RMSsum = 0.0L;
	size_t count = 0;
	for(MeasuredSED::const_iterator i=sed.begin(); i!=sed.end(); ++i)
	{
		long double flux = i->second.FluxDensity(pol);
		if(std::isfinite(flux))
		{
			++count;
			RMSsum += flux * flux;
		}
	}
	return sqrtl(RMSsum / (long double)(count));
}

double SEDAnalyser::sourceRMS(const MeasuredSED& sed, long double a, long double b, long double c)
{
	long double RMSsum = 0.0L;
	size_t count = 0;
	for(MeasuredSED::const_iterator i=sed.begin(); i!=sed.end(); ++i)
	{
		long double flux = i->second.FluxDensity(aocommon::Polarization::StokesI);
		long double x = i->second.FrequencyHz();
		long double powerLawValue = powl(b*x + c*x*x, a);
		if(std::isfinite(flux))
		{
			long double diff = flux - powerLawValue;
			++count;
			RMSsum += diff * diff;
		}
	}
	return sqrtl(RMSsum / (long double)(count));
}

double SEDAnalyser::sourceRMS(const MeasuredSED& sed, const aocommon::UVector<float>& terms)
{
	long double RMSsum = 0.0L;
	size_t count = 0;
	for(MeasuredSED::const_iterator i=sed.begin(); i!=sed.end(); ++i)
	{
		long double flux = i->second.FluxDensity(aocommon::Polarization::StokesI);
		long double x = i->second.FrequencyHz();
		long double powerLawValue = NonLinearPowerLawFitter::Evaluate(x, terms, _referenceFrequency);
		if(std::isfinite(flux))
		{
			long double diff = flux - powerLawValue;
			++count;
			RMSsum += diff * diff;
		}
	}
	return sqrtl(RMSsum / (long double)(count));
}

double SEDAnalyser::diffRMS(const MeasuredSED& sed)
{
	long double RMSsum = 0.0L;
	size_t count = 0;
	long double prevFlux=std::numeric_limits<long double>::quiet_NaN();
	for(MeasuredSED::const_iterator i=sed.begin(); i!=sed.end(); ++i)
	{
		long double flux = i->second.FluxDensity(aocommon::Polarization::StokesI);
		if(std::isfinite(flux) && std::isfinite(prevFlux))
		{
			long double diff = flux - prevFlux;
			++count;
			RMSsum += diff * diff;
		}
		prevFlux = flux;
	}
	// Since we have calculated the RMS of the difference between
	// two Gaussian distributions, the estimate of the RMS^2 of
	// the single distribution is half that of the difference dist.
	return sqrtl(0.5 * RMSsum / (long double)(count));
}

void SEDAnalyser::make_RA_vs_SI_plot()
{
	GNUPlot plot("si-vs-ra", "Right Ascension", "Spectral index");
	GNUPlot::Line* line = plot.AddPointSet("si-vs-ra-data.txt", "");
	double firstRA = _model->Source(0).MeanRA()*12.0/M_PI;
	bool aroundZero = firstRA<3.0 || firstRA>21.0;
	for(std::vector<SourceInfo>::const_iterator i=_sources.begin();
			i!=_sources.end(); ++i)
	{
		const SourceInfo& sInfo = *i;
		const ModelSource s = _model->Source(sInfo.index);
		double ra = s.MeanRA()*12.0/M_PI;
		if(aroundZero && ra > 12.0) ra -= 24.0;
		line->AddPoint(ra, sInfo.plExponent);
	}
}

void SEDAnalyser::make_Dec_vs_SI_plot()
{
	GNUPlot plot("si-vs-dec", "Declination", "Spectral index");
	GNUPlot::Line* line = plot.AddPointSet("si-vs-dec-data.txt", "");
	for(std::vector<SourceInfo>::const_iterator i=_sources.begin();
			i!=_sources.end(); ++i)
	{
		const SourceInfo& sInfo = *i;
		const ModelSource s = _model->Source(sInfo.index);
		line->AddPoint(s.MeanDec()*180.0/M_PI, sInfo.plExponent);
	}
}

void SEDAnalyser::make_SI_vs_brightness_plot()
{
	GNUPlot plot("si-vs-brightness", "Brightness", "Spectral index", true, false);
	GNUPlot::Line* line = plot.AddPointSet("si-vs-brightness-data.txt", "");
	for(std::vector<SourceInfo>::const_iterator i=_sources.begin();
			i!=_sources.end(); ++i)
	{
		const SourceInfo& sInfo = *i;
		line->AddPoint(sInfo.modalFlux, sInfo.plExponent);
	}
}

void SEDAnalyser::write_csv_file()
{
	std::ofstream csvFile("si.csv");
	
	csvFile << "RA,Dec,RA2,Brightness,SI,RMS,DiffRMS,FracRMS,RMS2,Term2,MaxErr,FracMaxErr,MaxStep,FracMaxStep,SubbandCorrelation,Distance,Q_RMS,U_RMS,V_RMS\n";
	
	for(std::vector<SourceInfo>::const_iterator i=_sources.begin();
			i!=_sources.end(); ++i)
	{
		const SourceInfo& sInfo = *i;
		const ModelSource& s = _model->Source(sInfo.index);

		double ra = s.MeanRA()*180.0/M_PI;
		
		csvFile
			<< ra << ','
			<< s.MeanDec()*180.0/M_PI << ','
			<< (ra<180.0 ? ra : (ra-360.0)) << ','
			<< sInfo.modalFlux << ','
			<< sInfo.plExponent << ','
			<< sInfo.rms << ','
			<< sInfo.diffRMS << ','
			<< sInfo.rms/sInfo.modalFlux << ','
			<< sInfo.rms2ndOrder << ','
			<< sInfo.pl2ndOrder << ','
			<< std::abs(sInfo.maxError) << ','
			<< std::abs(sInfo.maxError)/sInfo.modalFlux << ','
			<< sInfo.maxStep << ','
			<< sInfo.maxStep/sInfo.modalFlux << ','
			<< sInfo.subbandCorrelation << ','
			<< sInfo.distance*180.0/M_PI << ','
			<< sInfo.qrms << ','
			<< sInfo.urms << ','
			<< sInfo.vrms
			<< '\n';
	}
}

double SEDAnalyser::maxStep(const MeasuredSED& sed)
{
	double maxStepValue = 0.0;
	if(sed.MeasurementCount() == 0)
		return 0.0;
	for(MeasuredSED::const_iterator m=sed.begin();
			std::next(m)!=sed.end(); ++m)
	{
		double l = m->second.FluxDensityFromIndex(0);
		MeasuredSED::const_iterator mr = std::next(m);
		double r = mr->second.FluxDensityFromIndex(0);
		double step = std::abs(l-r);
		if(step > maxStepValue)
			maxStepValue = step;
	}
	return maxStepValue;
}

void SEDAnalyser::MakePSPlot()
{
	const Model& model = *_model;
	const MeasuredSED firstSED = model.Source(0).GetIntegratedMSED();
	aocommon::UVector<double>
		xValues(firstSED.MeasurementCount());
	double maxX = SpectrumFT::GetMaxKParallell(firstSED);
	MeasuredSED::const_iterator sedIter = firstSED.begin();
	for(size_t i=0; i!=firstSED.MeasurementCount(); ++i)
	{
		xValues[i] = maxX*double(i)/firstSED.MeasurementCount();
		++sedIter;
	}
	GNUStatPlot
		ftStatNoNormPlot, ftStatNormPlot,
		ftStatNoNormNoFitPlot, ftStatNormNoFitPlot,
		ftStatBinnedPlot, ftSmoothPlot;
	ftStatNoNormPlot.SetXValues(xValues);
	ftStatNormPlot.SetXValues(xValues);
	//ftStatNormPlot.SetDrawExtremes(true, true);
	ftStatNoNormNoFitPlot.SetXValues(xValues);
	ftStatNormNoFitPlot.SetXValues(xValues);
	ftStatBinnedPlot.SetDrawExtremes(true, true);
	ftStatBinnedPlot.SetDrawStdDev(true, false);
	ftSmoothPlot.SetXValues(xValues);
	SpectrumFT sft(xValues.size()), sftSmooth(xValues.size());
	sft.SetMetaData(firstSED, 2900.0, model.SourceCount());
	sftSmooth.SetMetaData(firstSED, 2900.0, model.SourceCount());
	sftSmooth.SetOutputFit(true);
	//sft.SetOutputFit(true); //testing
	for(size_t i=0; i!=model.SourceCount(); ++i)
	{
		LogBinnedPlot binned;
		
		const ModelSource& source = model.Source(i);
		const MeasuredSED sed = source.GetIntegratedMSED();
		aocommon::UVector<double> ftValues;
		sft.GetFTPower(ftValues, sed, 3, _psInTemperature, _useLombPeriodogram);
		ftStatNoNormPlot.AddYSet(ftValues);
		binned.Add(xValues, ftValues);
		
		aocommon::UVector<double> binnedXVal, binnedYVal;
		binned.GetData(binnedXVal, binnedYVal);
		if(!ftStatBinnedPlot.HasXValues())
			ftStatBinnedPlot.SetXValues(binnedXVal);
		ftStatBinnedPlot.AddYSet(binnedYVal);
		
		SpectrumFT::Normalize(ftValues);
		ftStatNormPlot.AddYSet(ftValues);
		
		sft.GetFTPower(ftValues, sed, 0, _psInTemperature, _useLombPeriodogram);
		ftStatNoNormNoFitPlot.AddYSet(ftValues);
		SpectrumFT::Normalize(ftValues);
		ftStatNormNoFitPlot.AddYSet(ftValues);
		
		sftSmooth.GetFTPower(ftValues, sed, 3, _psInTemperature, _useLombPeriodogram);
		ftSmoothPlot.AddYSet(ftValues);
	}
	
	const bool logX = true;
	
	GNUPlot binnedPlot("ps-binned", "k\u2225 (h/Mpc)", psYAxisDesc(), logX, true);
	binnedPlot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
	ftStatBinnedPlot.Plot(binnedPlot);
	
	GNUPlot noNormPlot("ps-fitted-unnormalized", "k\u2225 (h/Mpc)", psYAxisDesc(), logX, true);
	noNormPlot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
	ftStatNoNormPlot.Plot(noNormPlot);
	
	GNUPlot normPlot("ps-fitted-normalized", "k\u2225 (h/Mpc)", "Flux (Jy)", logX, true);
	normPlot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
	normPlot.SetYRange(1e-1, 10);
	normPlot.DrawLine(SpectrumFT::LowerKParallellEoRWindow(),1e-1,SpectrumFT::LowerKParallellEoRWindow(),10.0);
	normPlot.DrawLine(SpectrumFT::UpperKParallellEoRWindow(),1e-1,SpectrumFT::UpperKParallellEoRWindow(),10.0);
	ftStatNormPlot.Plot(normPlot);
	
	GNUPlot noNormNoFitPlot("ps-nofit-unnormalized", "k\u2225 (h/Mpc)", psYAxisDesc(), logX, true);
	noNormNoFitPlot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
	ftStatNoNormNoFitPlot.Plot(noNormNoFitPlot);
	
	GNUPlot normNoFitPlot("ps-nofit-normalized", "k\u2225 (h/Mpc)", "Flux (Jy)", logX, true);
	normNoFitPlot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
	normNoFitPlot.DrawLine(SpectrumFT::LowerKParallellEoRWindow(),1e-1,SpectrumFT::LowerKParallellEoRWindow(),10.0);
	normNoFitPlot.DrawLine(SpectrumFT::UpperKParallellEoRWindow(),1e-1,SpectrumFT::UpperKParallellEoRWindow(),10.0);
	ftStatNormNoFitPlot.Plot(normNoFitPlot);
	
	GNUPlot smoothPlot("ps-smooth", "k\u2225 (h/Mpc)", psYAxisDesc(), logX, true);
	sftSmooth.PreparePSPlot(smoothPlot, firstSED.CentreFrequency(), false);
	smoothPlot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
	ftSmoothPlot.Plot(smoothPlot);
}

double SEDAnalyser::getSourceSetWeightSum(const std::set<size_t>& sourceIndices)
{
	double weightSum = 0.0;
	for(std::set<size_t>::const_iterator iter=sourceIndices.begin(); iter!=sourceIndices.end(); ++iter)
	{
		const SourceInfo& info = _sources[*iter];
		double weight = 1.0/(info.diffRMS*info.diffRMS);
		weightSum += weight;
	}
	return weightSum;
}

void SEDAnalyser::getAveragedSED(aocommon::UVector<double>& averagedValues, const std::set<size_t>& sourceIndices, bool subtractModel, bool useSmoothSpectra)
{
	const Model& model = *_model;
	aocommon::UVector<double>
		sum(_model->Source(0).GetIntegratedMSED().MeasurementCount(), 0.0),
		weightSum(sum.size(), 0.0);
	for(std::set<size_t>::const_iterator iter=sourceIndices.begin(); iter!=sourceIndices.end(); ++iter)
	{
		const SourceInfo& info = _sources[*iter];
		const ModelSource& source = model.Source(info.index);
		const MeasuredSED sed = source.GetIntegratedMSED();
		aocommon::UVector<float> terms;
		if(subtractModel)
			terms = info.terms;
		if(terms.size() > _nTermsInPS)
			terms.resize(_nTermsInPS);
		else if(terms.size() < _nTermsInPS)
		{
			sed.FitLogPolynomial(terms, _nTermsInPS, aocommon::Polarization::StokesI, _referenceFrequency);
		}
		size_t channel = 0;
		double weight;
		weight = 1.0/(info.diffRMS*info.diffRMS);
		for(MeasuredSED::const_iterator m=sed.begin(); m!=sed.end(); ++m)
		{
			double fitVal = (subtractModel || useSmoothSpectra) ? (NonLinearPowerLawFitter::Evaluate(m->first, terms, _referenceFrequency)) : 0.0;
			const Measurement& meas = m->second;
			double flux = meas.FluxDensityFromIndex(0);
			if(std::isfinite(flux) || (useSmoothSpectra && _noMissing))
			{
				double v = useSmoothSpectra ? fitVal : (flux-fitVal);
				sum[channel] += v * weight;
				weightSum[channel] += weight;
			}
			++channel;
		}
	}
	averagedValues.resize(sum.size());
	for(size_t ch=0; ch!=averagedValues.size(); ++ch)
	{
		averagedValues[ch] = sum[ch] / weightSum[ch];
	}
}

void SEDAnalyser::getAveragedSED(MeasuredSED& averagedSED, const std::set<size_t>& sourceIndices, bool subtractModel, bool useSmoothSpectra)
{
	// Initialize it so it has appropriate channels
	averagedSED = _model->Source(0).GetIntegratedMSED();
	
	aocommon::UVector<double> averagedValues;
	getAveragedSED(averagedValues, sourceIndices, subtractModel, useSmoothSpectra);
	size_t ch = 0;
	for(MeasuredSED::iterator m=averagedSED.begin(); m!=averagedSED.end(); ++m)
	{
		m->second.SetFluxDensityFromIndex(0, averagedValues[ch]);
		++ch;
	}
}

void SEDAnalyser::getAveragedPSData(aocommon::UVector<double>& power, aocommon::UVector<double>& powerSigma, size_t nTermsFitted, const SpectrumFT& sft, bool useSmoothSpectra)
{
	std::set<size_t> sourceIndices;
	for(size_t i=0; i!=_sources.size(); ++i)
		sourceIndices.insert(i);
	
	MeasuredSED averagedSED;
	getAveragedSED(averagedSED, sourceIndices, nTermsFitted!=0, useSmoothSpectra);
	
	double measuredSigma = diffRMS(averagedSED);
	std::cout << "Measured noise level=" << measuredSigma << " Jy (RMS=" << sourceRMS(averagedSED) << " Jy).\n";
	sft.GetFTPower(power, averagedSED, 0, _psInTemperature, _useLombPeriodogram);
}

void SEDAnalyser::getAveragedSEDData(aocommon::UVector<double>& averageResidualFlux, aocommon::UVector<double>& fluxSigma, size_t nTermsFitted, bool useSmoothSpectra)
{
	std::set<size_t> sourceIndices;
	for(size_t i=0; i!=_sources.size(); ++i)
		sourceIndices.insert(i);
	
	getAveragedSED(averageResidualFlux, sourceIndices, nTermsFitted!=0, useSmoothSpectra);
	
	fluxSigma.assign(averageResidualFlux.size(), 0.0);
	const Model& model = *_model;
	double weightSum = 0.0;
	for(size_t i=0; i!=_sources.size(); ++i)
	{
		const SourceInfo& info = _sources[i];
		const ModelSource& source = model.Source(info.index);
		const MeasuredSED sed = source.GetIntegratedMSED();
		aocommon::UVector<float> terms;
		if(nTermsFitted != 0)
			terms = info.terms;
		size_t channel = 0;
		double weight = 1.0/(info.diffRMS*info.diffRMS);
		weightSum += weight;
		for(MeasuredSED::const_iterator m=sed.begin(); m!=sed.end(); ++m)
		{
			double fitVal = nTermsFitted!=0 ? (NonLinearPowerLawFitter::Evaluate(m->first, terms, _referenceFrequency)) : 0.0;
			const Measurement& meas = m->second;
			double val = meas.FluxDensityFromIndex(0)-fitVal;
			double mean = averageResidualFlux[channel];
			fluxSigma[channel] += (val-mean) * (val-mean) * weight;
			++channel;
		}
	}
	for(size_t i=0; i!=fluxSigma.size(); ++i)
		fluxSigma[i] = sqrt(fluxSigma[i] / weightSum);
}

void SEDAnalyser::MakeAveragedPSPlot()
{
	for(size_t imageSmoothed=0; imageSmoothed!=2; ++imageSmoothed)
	{
		bool doImgSmooth = imageSmoothed==1;
		if(doImgSmooth)
			std::cout << "Making smooth PS...\n";
		else
			std::cout << "Making averaged PS...\n";
		const std::string prefixName = doImgSmooth ? "smoothps-averaged" : "ps-averaged";
		const std::string sedPrefix = doImgSmooth ? "smoothsed-averaged" : "sed-averaged";
		
		const Model& model = *_model;
		MeasuredSED firstSED = model.Source(0).GetIntegratedMSED();
		aocommon::UVector<double>
			xValues(firstSED.MeasurementCount()),
			frequencies(firstSED.MeasurementCount());
		double maxX = SpectrumFT::GetMaxKParallell(firstSED);
		MeasuredSED::const_iterator sedIter = firstSED.begin();
		for(size_t i=0; i!=firstSED.MeasurementCount(); ++i)
		{
			xValues[i] = maxX*double(i)/firstSED.MeasurementCount();
			frequencies[i] = sedIter->first;
			++sedIter;
		}
		SpectrumFT sft(xValues.size());
		sft.SetMetaData(frequencies, 2900.0, model.SourceCount());
		//sft.SetOutputFit(true); //testing
		aocommon::UVector<double> power, powerSigma;
		getAveragedPSData(power, powerSigma, _nTermsInPS, sft, doImgSmooth);
		GNUPlot plot(prefixName, "k\u2225 (h/Mpc)", psYAxisDesc(), true, true);
		LogBinnedPlot bins;
		plot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
		GNUPlot::Line* line = plot.AddLine(prefixName + "-power.txt", "power");
		for(size_t i=0; i!=xValues.size(); ++i)
		{
			line->AddPoint(xValues[i], power[i]);
			bins.Add(xValues[i], power[i], 1.0);
		}
		
		GNUPlot binnedPlot(prefixName + "-binned", "k\u2225 (h/Mpc)", psYAxisDesc(), true, true);
		binnedPlot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
		GNUPlot::Line* binnedLine = plot.AddLine(prefixName + "-binned-power.txt", "power");
		bins.Plot(*binnedLine);
		
		aocommon::UVector<double> avgSED;
		getAveragedSEDData(avgSED, powerSigma, _nTermsInPS, doImgSmooth);
		GNUPlot plotAvgSED(sedPrefix, "Frequency (MHz)", "Averaged residual (mJy)", false, false);
		plotAvgSED.SetXRange(firstSED.begin()->first*1e-6, firstSED.rbegin()->first*1e-6);
		line = plotAvgSED.AddErrorSet(sedPrefix + "-data.txt", "");
		MeasuredSED::const_iterator m=firstSED.begin();
		for(size_t i=0; i!=avgSED.size(); ++i, ++m)
		{
			line->AddPoint(m->first*1e-6, avgSED[i]*1000.0, powerSigma[i]/sqrt(model.SourceCount())*1000.0);
		}
	}
}

void SEDAnalyser::runNoiseSimulation(aocommon::UVector<double>& simPower, size_t nTermsFitted, const MeasuredSED& templateSED, const SpectrumFT& sft, bool makePsf, bool simulateWithMissingData) const
{
	const Model& model = *_model;
	std::normal_distribution<double> normalDist(0.0, 1.0);
	
	aocommon::UVector<double> simSum(templateSED.MeasurementCount(), 0.0);
	double weightSum = 0.0;
	for(size_t i=0; i!=_sources.size(); ++i)
	{
		const SourceInfo& info = _sources[i];
		const ModelSource& source = model.Source(info.index);
		MeasuredSED simulatedSED = source.GetIntegratedMSED();
		
		// Simulate a new SED
		for(MeasuredSED::iterator m=simulatedSED.begin(); m!=simulatedSED.end(); ++m)
		{
			Measurement& meas = m->second;
			double plVal = (nTermsFitted==0) ? 0.0 : NonLinearPowerLawFitter::Evaluate(m->first, info.terms, _referenceFrequency);
			double noiseVal = normalDist(_rng) * info.diffRMS;
			double newVal = plVal + noiseVal;
			if(!simulateWithMissingData || std::isfinite(meas.FluxDensityFromIndex(0)))
				meas.SetFluxDensityFromIndex(0, newVal);
		}
		
		// continue generating normal PS
		aocommon::UVector<float> terms;
		if(nTermsFitted != 0)
			simulatedSED.FitLogPolynomial(terms, nTermsFitted, aocommon::Polarization::StokesI, _referenceFrequency);
		size_t channel = 0;
		double weight = 1.0/(info.diffRMS * info.diffRMS);
		for(MeasuredSED::const_iterator m=simulatedSED.begin(); m!=simulatedSED.end(); ++m)
		{
			double fitVal = nTermsFitted!=0 ? (NonLinearPowerLawFitter::Evaluate(m->first, terms, _referenceFrequency)) : 0.0;
			const Measurement& meas = m->second;
			double newVal = meas.FluxDensityFromIndex(0)-fitVal;
			simSum[channel] += newVal * weight;
			++channel;
		}
		weightSum += weight;
	}
	size_t channel = 0;
	MeasuredSED averageSED = templateSED;
	for(MeasuredSED::iterator m=averageSED.begin(); m!=averageSED.end(); ++m)
	{
		m->second.SetFluxDensityFromIndex(0, simSum[channel]/weightSum);
		++channel;
	}
	double simulatedSigma = diffRMS(averageSED);
	std::cout << "Simulated noise level=" << simulatedSigma << " Jy.\n";
	if(makePsf)
	{
		if(!std::isfinite(simulatedSigma))
			simulatedSigma = 1.0;
		for(MeasuredSED::iterator m=averageSED.begin(); m!=averageSED.end(); ++m)
		{
			Measurement& meas = m->second;
			if(std::isfinite(meas.FluxDensityFromIndex(0)))
				meas.SetFluxDensityFromIndex(0, simulatedSigma);
		}
	}
	sft.GetFTPower(simPower, averageSED, 0, _psInTemperature, _useLombPeriodogram);
}

void SEDAnalyser::MakeNoisePSPlot()
{
	const Model& model = *_model;
	MeasuredSED firstSED = model.Source(0).GetIntegratedMSED();
	aocommon::UVector<double>
		xValues(firstSED.MeasurementCount()),
		frequencies(firstSED.MeasurementCount());
	double maxX = SpectrumFT::GetMaxKParallell(firstSED);
	MeasuredSED::const_iterator sedIter = firstSED.begin();
	for(size_t i=0; i!=firstSED.MeasurementCount(); ++i)
	{
		xValues[i] = maxX*double(i)/firstSED.MeasurementCount();
		frequencies[i] = sedIter->first;
		++sedIter;
	}
	
	double weightSum = 0.0;
	double noiseSum = 0.0, rmsSum = 0.0;
	for(size_t i=0; i!=_sources.size(); ++i)
	{
		const SourceInfo& info = _sources[i];
		const ModelSource& source = model.Source(info.index);
		MeasuredSED noiseSED = source.GetIntegratedMSED();
		double weight;
		weight = 1.0/(info.diffRMS*info.diffRMS);
		noiseSum += info.diffRMS*info.diffRMS * weight* weight;
		rmsSum += info.rms*info.rms * weight*weight;
		weightSum += weight;
	}
	double noiseSigma = sqrt(noiseSum)/weightSum;
	double rmsSigmaNoDiff = sqrt(rmsSum)/weightSum;
	std::cout << "Theoretical noise level=" << noiseSigma << " Jy.\n";
	std::cout << "( " << rmsSigmaNoDiff << " Jy without using differential RMS).\n";
	MeasuredSED noiseSED = firstSED;
	size_t channel = 0;
	aocommon::UVector<double> noisePower(firstSED.MeasurementCount());
	for(MeasuredSED::iterator m=noiseSED.begin(); m!=noiseSED.end(); ++m)
	{
		m->second.SetFluxDensityFromIndex(0, noiseSigma);
		noisePower[channel] = noiseSigma;
		++channel;
	}
	
	SpectrumFT sft(xValues.size());
	sft.SetMetaData(frequencies, 2900.0, model.SourceCount());
	//sft.SetOutputFit(true); //testing
	for(size_t i=0; i!=noisePower.size(); ++i)
		noisePower[i] *= noisePower[i] * 2.0;
	if(_psInTemperature)
		sft.ConvertPSToTemperature(noisePower);
	
	aocommon::UVector<double> measuredData, sigmas;
	getAveragedPSData(measuredData, sigmas, _nTermsInPS, sft, false);
		
	aocommon::UVector<double> summedPower(firstSED.MeasurementCount(), 0.0);
	aocommon::UVector<double> summedPowerSq(firstSED.MeasurementCount(), 0.0);
	aocommon::UVector<double> summedPowerNoiseOnly(firstSED.MeasurementCount(), 0.0);
	aocommon::UVector<double> summedPowerNoiseOnlySq(firstSED.MeasurementCount(), 0.0);
	aocommon::UVector<double> summedPowerPSF(firstSED.MeasurementCount(), 0.0);
	
	const size_t NumberOfSimulations = 50;
	for(size_t simIter=0; simIter!=NumberOfSimulations; ++simIter)
	{
		std::cout << "==SIMULATION " << (simIter+1) << "==\n";
		aocommon::UVector<double> simPower, simPowerNoiseOnly, simPowerPSF;
		runNoiseSimulation(simPower, _nTermsInPS, firstSED, sft, false, true);
		runNoiseSimulation(simPowerNoiseOnly, 0, firstSED, sft, false, true);
		runNoiseSimulation(simPowerPSF, 0, firstSED, sft, true, true);
		for(size_t i=0; i!=simPower.size(); ++i)
		{
			summedPower[i] += simPower[i];
			summedPowerSq[i] += simPower[i]*simPower[i];
			summedPowerNoiseOnly[i] += simPowerNoiseOnly[i];
			summedPowerNoiseOnlySq[i] += simPowerNoiseOnly[i]*simPowerNoiseOnly[i];
			summedPowerPSF[i] += simPowerPSF[i];
		}
		
		GNUPlot plot("ps-averaged-simulated", "k\u2225 (h/Mpc)", psYAxisDesc(), true, true);
		sft.PreparePSPlot(plot, firstSED.CentreFrequency(), false);
		plot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
		GNUPlot::Line* simulatedLine = plot.AddLine("ps-averaged-simulated-power.txt", "Simulated");
		GNUPlot::Line* simNoiseOnlyLine = plot.AddLine("ps-averaged-simulated-nofit.txt", "Noise-only simulation");
		GNUPlot::Line* noiseLine = plot.AddLine("ps-averaged-simulated-noise.txt", "Noise estimate");
		size_t binsPer10 = 25;
		LogBinnedPlot simulatedBin(binsPer10), simNoiseOnlyBin(binsPer10), noiseBin(binsPer10), measuredBin(binsPer10), psfBin(binsPer10);
		for(size_t i=0; i!=xValues.size(); ++i)
		{
			double n = (simIter+1);
			double sum = summedPower[i], sumSq = summedPowerSq[i];
			double val = sum / n;
			double variance =
				(sumSq - (sum * sum / n)) / (n-1.0);
				
			double sumNO = summedPowerNoiseOnly[i], sumSqNO = summedPowerNoiseOnlySq[i];
			double valNO = sumNO / n;
			double varianceNO =
				(sumSqNO - (sumNO * sumNO / n)) / (n-1.0);
				
			double valPSF = summedPowerPSF[i] / n;
				
			simulatedLine->AddPoint(xValues[i], val, sqrt(variance));
			simulatedBin.Add(xValues[i], val, variance);
			simNoiseOnlyBin.Add(xValues[i], valNO, sqrt(varianceNO));
			simNoiseOnlyLine->AddPoint(xValues[i], valNO, sqrt(varianceNO));
			noiseLine->AddPoint(xValues[i], noisePower[i]);
			noiseBin.Add(xValues[i], noisePower[i], 1.0);
			measuredBin.Add(xValues[i], measuredData[i], 1.0);
			psfBin.Add(xValues[i], valPSF, 1.0);
		}
		plot.AddLineFromExistingFile("ps-averaged-power.txt", "Measured power");
		
		GNUPlot binnedPlot("ps-averaged-simulated-binned", "k\u2225 (h/Mpc)", psYAxisDesc(), true, true);
		sft.PreparePSPlot(binnedPlot, firstSED.CentreFrequency(), false);
		binnedPlot.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
		GNUPlot::Line* measBinLine = binnedPlot.AddErrorSet("ps-averaged-simulated-binned-measured.txt", "Measured power", 5.0);
		measuredBin.PlotWithOtherErrorBars(*measBinLine, simulatedBin);
		GNUPlot::Line* simBinLine = binnedPlot.AddLine("ps-averaged-simulated-binned-power.txt", "Simulated power");
		simulatedBin.Plot(*simBinLine);
		//GNUPlot::Line* simOnlyBinLine = binnedPlot.AddLine("ps-averaged-simulated-no-binned-power.txt", "Simulated power (noise only)");
		//simNoiseOnlyBin.Plot(*simOnlyBinLine);
		GNUPlot::Line* noiseBinLine = binnedPlot.AddLine("ps-averaged-simulated-binned-noise.txt", "Noise power");
		noiseBin.Plot(*noiseBinLine);
		//GNUPlot::Line* psfBinLine = binnedPlot.AddLine("ps-averaged-simulated-binned-psf.txt", "PSF");
		//psfBin.Plot(*psfBinLine);
	}
}

void SEDAnalyser::removeEdgeChannels()
{
	for(size_t i=0; i!=_model->SourceCount(); ++i)
	{
		ModelSource& source = _model->Source(i);
		for(ModelSource::iterator c=source.begin(); c!=source.end(); ++c)
		{
			MeasuredSED& sed = c->MSED();
			size_t channel = 0;
			for(MeasuredSED::iterator m=sed.begin(); m!=sed.end(); ++m)
			{
				bool isEdgeChannel =
					channel%32<_edgeChannels || channel%32>=32-_edgeChannels;
				if(isEdgeChannel)
				{
					Measurement& meas = m->second;
					meas.SetFluxDensityFromIndex(0, std::numeric_limits<long double>::quiet_NaN());
				}
				++channel;
			}
		}
	}
}

void SEDAnalyser::removeEndChannels()
{
	for(size_t i=0; i!=_model->SourceCount(); ++i)
	{
		ModelSource& source = _model->Source(i);
		for(ModelSource::iterator c=source.begin(); c!=source.end(); ++c)
		{
			MeasuredSED& sed = c->MSED();
			size_t channel = 0;
			for(MeasuredSED::iterator m=sed.begin(); m!=sed.end(); ++m)
			{
				bool isEdgeChannel = channel >= sed.MeasurementCount()-_flaggedEndChannels;
				if(isEdgeChannel)
				{
					Measurement& meas = m->second;
					meas.SetFluxDensityFromIndex(0, std::numeric_limits<long double>::quiet_NaN());
				}
				++channel;
			}
		}
	}
}

void SEDAnalyser::removeBadChannels()
{
	for(size_t i=0; i!=_model->SourceCount(); ++i)
	{
		ModelSource& source = _model->Source(i);
		for(ModelSource::iterator c=source.begin(); c!=source.end(); ++c)
		{
			MeasuredSED& sed = c->MSED();
			for(MeasuredSED::iterator m=sed.begin(); m!=sed.end(); ++m)
			{
				Measurement& meas = m->second;
				bool isBad = std::abs(meas.FluxDensity(aocommon::Polarization::StokesI)) > 100.0;
				if(isBad)
				{
					meas.SetFluxDensityFromIndex(0, std::numeric_limits<long double>::quiet_NaN());
				}
			}
		}
	}
}

void SEDAnalyser::write_cath_csv_file()
{
	std::ofstream csvFile("seds-cath.txt");
	std::ofstream csvFileNO("seds-cath-noise-only.txt");
	std::ofstream csvFileModel("seds-cath-smoothmodel.txt");
	std::ofstream csvFileSIModel("seds-cath-smoothmodel-si-only.txt");
	csvFile.precision(16);
	csvFileNO.precision(16);
	csvFileModel.precision(16);
	csvFileSIModel.precision(16);
	std::normal_distribution<double> normalDist(0.0, 1.0);
	
	// Columns: RA Dec L M RMS Freq1(Hz) Ch1(mJy) Freq2(Hz) Ch2 ... FreqN ChN
	
	for(std::vector<SourceInfo>::const_iterator i=_sources.begin();
			i!=_sources.end(); ++i)
	{
		const SourceInfo& sInfo = *i;
		const ModelSource& s = _model->Source(sInfo.index);

		double raRad = s.MeanRA(), decRad = s.MeanDec();
		double phaseCentreRA = 0.0, phaseCentreDec = -27.0*M_PI/180.0;
		double l,m;
		size_t nTermsFitted = sInfo.terms.size();
		aocommon::ImageCoordinates::RaDecToLM<double>(raRad, decRad, phaseCentreRA, phaseCentreDec, l, m);
		
		csvFile
			<< raRad*180.0/M_PI << ' '
			<< decRad*180.0/M_PI << ' '
			<< l << ' ' << m << ' '
			<< sInfo.diffRMS*1000.0;
		csvFileNO
			<< raRad*180.0/M_PI << ' '
			<< decRad*180.0/M_PI << ' '
			<< l << ' ' << m << ' '
			<< sInfo.diffRMS*1000.0;
		csvFileModel
			<< raRad*180.0/M_PI << ' '
			<< decRad*180.0/M_PI << ' '
			<< l << ' ' << m << ' '
			<< sInfo.diffRMS*1000.0;
		csvFileSIModel
			<< raRad*180.0/M_PI << ' '
			<< decRad*180.0/M_PI << ' '
			<< l << ' ' << m << ' '
			<< sInfo.diffRMS*1000.0;
		MeasuredSED sed = s.GetIntegratedMSED();
		for(MeasuredSED::const_iterator mi=sed.begin(); mi!=sed.end(); ++mi)
		{
			double fitVal = nTermsFitted!=0 ? (NonLinearPowerLawFitter::Evaluate(mi->first, sInfo.terms, _referenceFrequency)) : 0.0;
			aocommon::UVector<float> termsWithoutSI(sInfo.terms);
			if(termsWithoutSI.size()>2)
				termsWithoutSI.resize(2);
			double fitValSI = nTermsFitted!=0 ? (NonLinearPowerLawFitter::Evaluate(mi->first, termsWithoutSI, _referenceFrequency)) : 0.0;
			const Measurement& m = mi->second;
			double v = m.FluxDensityFromIndex(0)-fitVal;
			csvFile << ' ' << m.FrequencyHz() << ' ' << (v*1000.0);
			double noiseVal = normalDist(_rng) * sInfo.diffRMS;
			if(std::isfinite(v))
			{
				csvFileNO << ' ' << m.FrequencyHz() << ' ' << (noiseVal*1000.0);
				csvFileModel << ' ' << m.FrequencyHz() << ' ' << (fitVal*1000.0);
				csvFileSIModel << ' ' << m.FrequencyHz() << ' ' << (fitValSI*1000.0);
			}
			else
			{
				// v is nan, so just write nans.
				csvFileNO << ' ' << m.FrequencyHz() << ' ' << (v*1000.0);
				csvFileModel << ' ' << m.FrequencyHz() << ' ' << (v*1000.0);
				csvFileSIModel << ' ' << m.FrequencyHz() << ' ' << (v*1000.0);
			}
		}
		csvFile << '\n';
		csvFileNO << '\n';
		csvFileModel << '\n';
		csvFileSIModel << '\n';
	}
}

void SEDAnalyser::measureSourceCountNoise()
{
	size_t sourceCount = 1;
	double sourceCountIncrease = 1.1;
	
	GNUPlot plot("source-count-behaviour", "Effective nr. of sources", "RMS", true, true);
	GNUPlot::Line *diffRMSLine = plot.AddPointSet("source-count-behaviour-diff-rms.txt", "Differential RMS");
	GNUPlot::Line *rmsLine = plot.AddPointSet("source-count-behaviour-residual-rms.txt", "Residual RMS");
	GNUPlot::Line *powerLine = plot.AddPointSet("source-count-behaviour-residual-power.txt", "Power");
	
	std::set<size_t> fullSet;
	for(size_t i=0; i!=_sources.size(); ++i)
		fullSet.insert(i);
	double fullWeight = getSourceSetWeightSum(fullSet);
	
	ProgressBar progress("Calculating noise for different source counts");
	
	double fittedRMS = 0.0, fittedRMSSourceCount = 0.0, fittedPower = 0.0, fittedPowerSourceCount = 0.0, fittedResidualRMS = 0.0, fittedResidualRMSSourceCount = 0.0, lastRMS = 0.0, lastResidualRMS = 0.0, lastPower = 0.0;

	MeasuredSED firstSED = _model->Source(0).GetIntegratedMSED();
	SpectrumFT sft(firstSED.MeasurementCount());
	sft.SetMetaData(firstSED, 2900.0, sourceCount);
			
	while(sourceCount <= _sources.size())
	{
		double avgRMS = 0.0, avgPower = 0.0, avgSourceCount = 0.0, avgResidualRMS = 0.0;
		size_t repeatCount = 10;
		for(size_t repeat=0; repeat!=repeatCount; ++repeat)
		{
			std::uniform_int_distribution<size_t> dist(0, _sources.size()-1);
			std::set<size_t> indexPermutation;
			while(indexPermutation.size() < sourceCount)
				indexPermutation.insert(dist(_rng));
			
			MeasuredSED averagedSED;
			getAveragedSED(averagedSED, indexPermutation, true, false);
			
			double permutationWeight = getSourceSetWeightSum(indexPermutation);
			double effectiveSourceCount = double(_sources.size()) * permutationWeight / fullWeight;
			
			double residualRMS = sourceRMS(averagedSED);
			double differentialRMS = diffRMS(averagedSED);
			
			aocommon::UVector<double> values;
			double stPower = sqrt(sft.GetAveragePowerInKRange(0.1, 0.2, averagedSED, 0, _psInTemperature, _useLombPeriodogram));
			
			diffRMSLine->AddPoint(effectiveSourceCount, differentialRMS);
			rmsLine->AddPoint(effectiveSourceCount, residualRMS);
			powerLine->AddPoint(effectiveSourceCount, stPower);
			
			avgRMS += differentialRMS;
			avgPower += stPower;
			avgSourceCount += effectiveSourceCount;
			avgResidualRMS += residualRMS;
		}
		avgRMS /= repeatCount;
		avgPower /= repeatCount;
		avgSourceCount /= repeatCount;
		avgResidualRMS /= repeatCount;
		if(sourceCount >= 8 && sourceCount <= 15)
		{
			fittedRMS = avgRMS;
			fittedRMSSourceCount = avgSourceCount;
			fittedPower = avgPower;
			fittedPowerSourceCount = avgSourceCount;
			fittedResidualRMSSourceCount = avgSourceCount;
			fittedResidualRMS = avgResidualRMS;
		}
		lastRMS = avgRMS;
		lastResidualRMS = avgResidualRMS;
		lastPower = avgPower;
		
		if(sourceCount != _sources.size())
		{
			size_t newSourceCount = floor(double(sourceCount) * sourceCountIncrease);
			if(newSourceCount > sourceCount)
				sourceCount = newSourceCount;
			else
				sourceCount++;
			if(sourceCount > _sources.size())
				sourceCount = _sources.size();
		}
		else ++sourceCount;
		
		progress.SetProgress(sourceCount, _sources.size());
	}
	
	double fittedLastRMS = fittedRMS / sqrt(double(_sources.size())/fittedRMSSourceCount);
	double fittedLastResidualRMS = fittedResidualRMS / sqrt(double(_sources.size())/fittedResidualRMSSourceCount);
	double fittedLastPower = fittedPower / sqrt(double(_sources.size())/fittedPowerSourceCount);
	std::cout <<
		"Ratio between best RMS and fitted RMS: " << lastRMS << '/' << fittedLastRMS << " = " << (lastRMS/fittedLastRMS) << 
		"\nRatio between best residual RMS and fitted residual RMS: " << lastResidualRMS << '/' << fittedLastResidualRMS << " = " << (lastResidualRMS/fittedLastResidualRMS) <<
		"\nRatio between best power and fitted power: " << lastPower << '/' << fittedLastPower << " = " << (lastPower/fittedLastPower) <<
		'\n';
	
	std::ostringstream s1;
	s1 << fittedRMS << "/sqrt(x/" << fittedRMSSourceCount << ")";
	plot.AddFunction(s1.str(), "{/Symbol \265} 1/\\{/Symbol \326}t");
	
	std::ostringstream s2;
	s2 << fittedPower << "/sqrt(x/" << fittedPowerSourceCount << ")";
	plot.AddFunction(s2.str(), "");
	
	std::ostringstream s3;
	s3 << fittedResidualRMS << "/sqrt(x/" << fittedResidualRMSSourceCount << ")";
	plot.AddFunction(s3.str(), "");
}
