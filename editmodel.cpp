#include "banddata.h"
#include "delaunay.h"
#include "units/imagecoordinates.h"
#include "loghistogram.h"
#include "rmsynthesis.h"
#include "uvector.h"

#include "model/model.h"
#include "model/powerlawsed.h"
#include "fitsreader.h"
#include "matrix2x2.h"
#include "gnuplot.h"
#include "spectrumft.h"
#include "ndppp.h"

#include "deconvolution/spectralfitter.h"

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <set>

std::string getSagecalSourceName(const ModelSource& s, const ModelComponent& c, size_t compIndex)
{
	std::stringstream str;
	switch(c.Type())
	{
		case ModelComponent::PointSource:
			str << std::string("P_") << s.Name() << "_C" << compIndex;
			return str.str();
		case ModelComponent::GaussianSource:
			str << std::string("G_") << s.Name() << "_C" << compIndex;
			return str.str();
		default:
			throw std::runtime_error("Invalid component type");
	}
}

void readImage(ao::uvector<double>& image, FitsReader& templateReader, const std::string& filename)
{
	FitsReader reader(filename);
	size_t
		width = reader.ImageWidth(),
		height = reader.ImageHeight();
	if(width != templateReader.ImageWidth() || height != templateReader.ImageHeight())
		throw std::runtime_error("Beam images have different sizes");
	image.resize(width*height);
	reader.Read(image.data());
}

bool compareWithBeam(const ModelSource& lhs, const ModelSource& rhs)
{
	double lhsVal = *static_cast<double*>(lhs.UserData());
	double rhsVal = *static_cast<double*>(rhs.UserData());
	return lhsVal < rhsVal;
}

struct SplitLine
{
	double ra1, dec1;
	double ra2, dec2;
	
	bool IsPointAboveLine(double ra, double dec) const
	{
		if(ra1 == ra2)
			return ( dec2 >= dec1 && ra > ra1 )  ||  ( dec1 > dec2 && ra <= ra1 );
		else {
			double alpha = (dec2 - dec1) / (ra2 - ra1);
			if(ra2 > ra1)
				return dec >= alpha * (ra - ra1) + dec1;
			else
				return dec < alpha * (ra - ra1) + dec1;
		}
	}
};

int main(int argc, char *argv[])
{
	if(argc == 1)
	{
		std::cout << "editmodel -- Interpolation, extrapolation, plotting and scaling of the \n"
		"spectral energy distribution. Usage:\n"
		"\teditmodel [-p [-ft]] [-rmp] [-m <output model>] [-o] [-s <scale>] [-s-to <freq> <flux>] [-sp <peakflux A> <freq A> <peakflux B> <freq B>] [-sc <intflux A> <freq A> <intflux B> <freq B>] [-set0/1/2/3 <flux>] [-unpolarized] [-pl] [-t <threshold>] [-tbeam <beamprefix> <threshold>] [-tc <compthreshold>] [-tcl <cluster threshold] [-r <new-nr-channels>] [-ravg <new-nr-channels>] [-near/outside <ra> <dec> <dist>] [-combine-diff-meas] [-collect <name>] [-uncollect] [-rnd <n> <ra> <dec> <dist>] [-sort] [-sortbeam <beamprefix>] [-lognlogs <frequency> <bincount>] [-stats] [-setfrequency <val>] [-delnans] [-sagecal <prefix> <chunksize>] [-dppp-model <filename>] [-skymodel <filename>] [-toapp <beamprefix>] [-save-clusters <clusters.ann>] [-list] [-evaluate <freq>] [-scale-to <model> <freq-start> <freq-end> <terms>] [-select <name>] [-search <count> <ra> <dec>] [-simuniform N RA Dec dist flux] [-simpopulation RA Dec dist lowflux highflux] [-min-separation <angle>] [-replace-si <n> <si> <terms...>] [-set-cluster <name>] "
		"[-rts <filename>] [-to-powerlaw <n>] [-shift <ra-angle> <dec-angle>] [-kvis <kvis .ann file>] [-split/-split2 <ra1> <dec1> <ra2> <dec2>] "
		"<model> [<more models...>]\n";
		return 0;
	}
	int argi = 1;
	bool outputPlot = false, doFT = false, outputRMPlot = false, outputCsv = false, outputSICsv = false, powerlaw = false, optimize = false, applyThreshold = false, applyIndexThreshold = false, applyComponentThreshold = false, applyClusterThreshold = false, resample = false, resampleByAveraging = false, unpolarized = false, delNaNs = false, outputList = false;
	double evaluateFrequency = 0.0;
	bool setPolarization[4] = {false, false, false, false};
	long double setPolFlux[4] = {0.0, 0.0, 0.0, 0.0};
	long double scale = 1.0, scaleToFreq = 0.0, scaleToFlux = 0.0, threshold = 0.0, logNlogSFrequency = 0.0;
	long double scalePeakA = 1.0, scaleFreqA = 0.0, scalePeakB = 1.0, scaleFreqB = 0.0;
	size_t newChannelCount = 0, logNlogSBinCount = 0, sagecalFirstClusterChunkSize = 1;
	std::string outputModel, collectName, csvFilename, plotTitle, sagecalPrefix, toAppBeamPrefix, dpppModelFilename, kvisAnnFile, sortBeamPrefix;
	bool nearFilter = false, outsideFilter = false, scalePeak = false, scaleSource = false, scaleWithSI = false, doCollect = false, doSort = false, setFrequency = false, uncollect = false;
	long double nearFilterRA = 0.0, nearFilterDec = 0.0, nearFilterDist = 0.0, setFrequencyValue = 0.0, powerSumLimit = 0.0, shiftRA = 0.0, shiftDec = 0.0;
	enum { AddFluxes, AverageFluxes, DifferentFrequencies } combineStrategy = AddFluxes;
	bool logNlogS = false, sourceStats = false;
	size_t rndN = 0, minMeasurements = 0, thresholdIndex = 0;
	double rndRA = 0.0, rndDec = 0.0, rndDist = 0.0;
	bool doSearch = false;
	size_t searchCount = 3;
	double searchRA = 0.0, searchDec = 0.0;
	size_t simUniformN = 0;
	double simUniformRa=0.0, simUniformDec=0.0, simUniformDist=0.0, simUniformFlux=0.0;
	double simPopulationRa=0.0, simPopulationDec=0.0, simPopulationDist=0.0, simPopulationFluxLo=0.0, simPopulationFluxHi=0.0;
	double simSIAverage = -0.7, simSIStdDev = 0.2;
	std::string saveClusters;
	std::string scaleToModelFilename;
	double scaleToFreqStart = 0.0, scaleToFreqEnd = 0.0;
	size_t scaleToNTerms = 0;
	double minSeparation = 0.0;
	bool replaceSI = false;
	ao::uvector<double> replaceSIvalues;
	std::string setClusterName;
	std::string rtsFilename;
	size_t toPowerlaw = 0;
	std::vector<SplitLine> splits;
	
	std::string selectSourceName;
	while(argi!=argc && argv[argi][0]=='-')
	{
		const std::string option(&argv[argi][1]);
		if(option == "p")
		{
			outputPlot = true;
		} else if(option == "ptitle") {
			++argi;
			outputPlot = true;
			plotTitle = argv[argi];
		} else if(option == "rmp")
		{
			outputRMPlot = true;
		} else if(option == "csv") {
			++argi;
			outputCsv = true;
			csvFilename = argv[argi];
		} else if(option == "si-csv") {
			++argi;
			outputSICsv = true;
			csvFilename = argv[argi];
		} else if(option == "list") {
			outputList = true;
		} else if(option == "evaluate") {
			++argi;
			evaluateFrequency = atof(argv[argi]);
		} else if(option == "scale-to") {
			scaleToModelFilename = argv[argi+1];
			scaleToFreqStart = atof(argv[argi+2]);
			scaleToFreqEnd = atof(argv[argi+3]);
			scaleToNTerms = atoi(argv[argi+4]);
			argi += 4;
		} else if(option == "save-clusters") {
			++argi;
			saveClusters = argv[argi];
		} else if(option == "search") {
			++argi;
			searchCount = atoi(argv[argi]);
			++argi;
			searchRA = RaDecCoord::ParseRA(argv[argi]);
			++argi;
			searchDec = RaDecCoord::ParseDec(argv[argi]);
			doSearch = true;
		} else if(option == "near")
		{
			++argi;
			nearFilterRA = RaDecCoord::ParseRA(argv[argi]);
			++argi;
			nearFilterDec = RaDecCoord::ParseDec(argv[argi]);
			++argi;
			nearFilterDist = Angle::Parse(argv[argi], "near filter distance", Angle::Degrees);
			nearFilter = true;
		} else if(option == "outside")
		{
			++argi;
			nearFilterRA = RaDecCoord::ParseRA(argv[argi]);
			++argi;
			nearFilterDec = RaDecCoord::ParseDec(argv[argi]);
			++argi;
			nearFilterDist = Angle::Parse(argv[argi], "outside filter distance", Angle::Degrees);
			outsideFilter = true;
		} else if(option == "rnd")
		{
			++argi;
			rndN = atof(argv[argi]);
			++argi;
			rndRA = RaDecCoord::ParseRA(argv[argi]);
			++argi;
			rndDec = RaDecCoord::ParseDec(argv[argi]);
			++argi;
			rndDist = Angle::Parse(argv[argi], "rnd", Angle::Degrees);
		} else if(option == "uncollect")
		{
			uncollect = true;
		} else if(option == "collect")
		{
			doCollect = true;
			++argi;
			collectName = argv[argi];
		} else if(option == "m")
		{
			++argi;
			outputModel = argv[argi];
		} else if(option == "s")
		{
			++argi;
			scale = atof(argv[argi]);
		} else if(option == "s-to")
		{
			++argi;
			scaleToFreq = atof(argv[argi]);
			++argi;
			scaleToFlux = atof(argv[argi]);
		} else if(option == "sp") // scale w.r.t. peak flux
		{
			scalePeak = true;
			++argi; scalePeakA = atof(argv[argi]);
			++argi; scaleFreqA = atof(argv[argi]) * 1000000.0;
			++argi; scalePeakB = atof(argv[argi]);
			++argi; scaleFreqB = atof(argv[argi]) * 1000000.0;
		} else if(option == "sc") // scale w.r.t. whole source flux
		{
			scaleSource = true;
			++argi; scalePeakA = atof(argv[argi]);
			++argi; scaleFreqA = atof(argv[argi]) * 1000000.0;
			++argi; scalePeakB = atof(argv[argi]);
			++argi; scaleFreqB = atof(argv[argi]) * 1000000.0;
		} else if(option == "sc-si") // scale w.r.t. whole source flux
		{
			scaleSource = true;
			scaleWithSI = true;
			++argi; scalePeakA = atof(argv[argi]);
			++argi; scaleFreqA = atof(argv[argi]) * 1000000.0;
			++argi; scalePeakB = atof(argv[argi]);
			++argi; scaleFreqB = atof(argv[argi]) * 1000000.0;
		} else if(option == "set0")
		{
			++argi;
			setPolarization[0] = true;
			setPolFlux[0] = atof(argv[argi]);
		} else if(option == "set1")
		{
			++argi;
			setPolarization[1] = true;
			setPolFlux[1] = atof(argv[argi]);
		} else if(option == "set2")
		{
			++argi;
			setPolarization[2] = true;
			setPolFlux[2] = atof(argv[argi]);
		} else if(option == "set3")
		{
			++argi;
			setPolarization[3] = true;
			setPolFlux[3] = atof(argv[argi]);
		} else if(option == "pl")
		{
			powerlaw = true;
		} else if(option == "t")
		{
			++argi;
			threshold = atof(argv[argi]);
			applyThreshold = true;
		} else if(option == "t-index")
		{
			++argi;
			thresholdIndex = atoi(argv[argi]);
			applyIndexThreshold = true;
		} else if(option == "tc")
		{
			++argi;
			threshold = atof(argv[argi]);
			applyComponentThreshold = true;
		} else if(option == "tcl")
		{
			++argi;
			threshold = atof(argv[argi]);
			applyClusterThreshold = true;
		} else if(option == "r")
		{
			++argi;
			newChannelCount = atoi(argv[argi]);
			resample = true;
		} else if(option == "ravg")
		{
			++argi;
			newChannelCount = atoi(argv[argi]);
			resampleByAveraging = true;
		} else if(option == "unpolarized")
		{
			unpolarized = true;
		} else if(option == "o")
		{
			optimize = true;
		} else if(option == "combine-diff-meas")
		{
			combineStrategy = DifferentFrequencies;
		} else if(option == "combine-avg-meas")
		{
			combineStrategy = AverageFluxes;
		} else if(option == "sort")
		{
			doSort = true;
		} else if(option == "sortbeam")
		{
			++argi;
			sortBeamPrefix = argv[argi];
		} else if(option == "lognlogs")
		{
			++argi;
			logNlogSFrequency = atof(argv[argi]);
			++argi;
			logNlogSBinCount = atoi(argv[argi]);
			logNlogS = true;
		} else if(option == "stats")
		{
			sourceStats = true;
		} else if(option == "delnans")
		{
			delNaNs = true;
		} else if(option == "plotft")
		{
			doFT = true;
		} else if(option == "setfrequency")
		{
			++argi;
			setFrequency = true;
			setFrequencyValue = atof(argv[argi]);
		} else if(option == "min-measurements")
		{
			++argi;
			minMeasurements = atoi(argv[argi]);
		} else if(option == "sagecal")
		{
			++argi;
			sagecalPrefix = argv[argi];
			++argi;
			sagecalFirstClusterChunkSize = atoi(argv[argi]);
		} else if(option == "skymodel" || option == "dppp-model")
		{
			++argi;
			dpppModelFilename = argv[argi];
		} else if(option == "kvis")
		{
			++argi;
			kvisAnnFile = argv[argi];
		} else if(option == "toapp")
		{
			++argi;
			toAppBeamPrefix = argv[argi];
		} else if(option == "select")
		{
			++argi;
			selectSourceName = argv[argi];
		} else if(option == "powersum")
		{
			++argi;
			powerSumLimit = atof(argv[argi]);
		} else if(option == "simuniform")
		{
			simUniformN = atoi(argv[argi+1]);
			simUniformRa = RaDecCoord::ParseRA(argv[argi+2]);
			simUniformDec = RaDecCoord::ParseDec(argv[argi+3]);
			simUniformDist = atof(argv[argi+4]) * (M_PI/180.0);
			simUniformFlux = atof(argv[argi+5]);
			argi+=5;
		} else if(option == "simpopulation")
		{
			simPopulationRa = RaDecCoord::ParseRA(argv[argi+1]);
			simPopulationDec = RaDecCoord::ParseDec(argv[argi+2]);
			simPopulationDist = atof(argv[argi+3]) * (M_PI/180.0);
			simPopulationFluxLo = atof(argv[argi+4]);
			simPopulationFluxHi = atof(argv[argi+5]);
			argi+=5;
		} else if(option == "min-separation")
		{
			++argi;
			minSeparation = Angle::Parse(argv[argi], "min-separation", Angle::Arcseconds);
		}
		else if(option == "replace-si")
		{
			++argi;
			replaceSI = true;
			size_t n = atoi(argv[argi]);
			for(size_t i=0; i!=n; ++i)
			{
				++argi;
				replaceSIvalues.emplace_back(atof(argv[argi]));
			}
		} else if(option == "set-cluster")
		{
			++argi;
			setClusterName = argv[argi];
		} else if(option == "rts")
		{
			++argi;
			rtsFilename = argv[argi];
		} else if(option == "to-powerlaw")
		{
			++argi;
			toPowerlaw = atoi(argv[argi]);
		} else if(option == "shift")
		{
			shiftRA = Angle::Parse(argv[argi+1], "shift RA", Angle::Degrees);
			shiftDec = Angle::Parse(argv[argi+2], "shift Dec",  Angle::Degrees);
			argi += 2;
		} else if(option == "split")
		{
			splits.emplace_back();
			SplitLine& line = splits.back();
			line.ra1 = RaDecCoord::ParseRA(argv[argi+1]);
			line.dec1 = RaDecCoord::ParseDec(argv[argi+2]);
			line.ra2 = RaDecCoord::ParseRA(argv[argi+3]);
			line.dec2 = RaDecCoord::ParseDec(argv[argi+4]);
			argi += 4;
		} else if(option == "split2")
		{
			splits.emplace_back();
			SplitLine& line = splits.back();
			line.ra2 = RaDecCoord::ParseRA(argv[argi+1]);
			line.dec2 = RaDecCoord::ParseDec(argv[argi+2]);
			line.ra1 = RaDecCoord::ParseRA(argv[argi+3]);
			line.dec1 = RaDecCoord::ParseDec(argv[argi+4]);
			argi += 4;
		} else {
			throw std::runtime_error(std::string("Unknown option given: -") + option);
		}
		++argi;
	}
	Model model;
	switch(combineStrategy)
	{
		case AddFluxes:
		case AverageFluxes:
			for(int modelIndex=argi; modelIndex!=argc; ++modelIndex)
			{
				model += Model(argv[modelIndex]);
			}
			break;
		case DifferentFrequencies:
			for(int modelIndex=argi; modelIndex!=argc; ++modelIndex)
			{
				model.CombineMeasurements(Model(argv[modelIndex]));
			}
			break;
	}
	if(combineStrategy == AverageFluxes)
	{
		double fact = 1.0 / (argc - argi);
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				if(!compPtr->HasMeasuredSED())
					throw std::runtime_error("Not a measured SED");
				MeasuredSED &sed = compPtr->MSED();
				for(MeasuredSED::iterator m=sed.begin(); m!=sed.end(); ++m)
				{
					for(size_t p=0; p!=4; ++p)
						m->second.SetFluxDensityFromIndex(p, m->second.FluxDensityFromIndex(p) * fact);
				}
			}
		}
	}
	if(scaleToNTerms > 0)
	{
		Model scaleModel(scaleToModelFilename);
		if(scaleModel.SourceCount() != 1)
			throw std::runtime_error("Scale reference model should have a single source");
		const ModelSource& source = scaleModel.Source(0);
		ao::uvector<double> values(scaleToNTerms), frequencies(scaleToNTerms), weights(scaleToNTerms, 1.0);
		double scaleToBandwidth = scaleToFreqEnd - scaleToFreqStart;
		std::cout << "Scaling model to have the following flux values:\n";
		std::cout << "Frequency\tTo\tOld\n";
		for(size_t term = 0; term!=scaleToNTerms; ++term)
		{
			double frequency = scaleToFreqStart + scaleToBandwidth * double(term) / double(scaleToNTerms-1);
			values[term] = source.TotalFlux(frequency, Polarization::StokesI);
			frequencies[term] = frequency;
			std::cout << round(frequency*1e-6) << " MHz\t" << values[term] << " Jy\t" <<
				model.Source(0).TotalFlux(frequency, Polarization::StokesI) << " Jy\n";
		}
		
		for(ModelSource& source : model)
		{
			ao::uvector<double> scalingFactors(values.size());
			for(size_t term = 0; term!=scaleToNTerms; ++term)
			{
				double frequency = frequencies[term];
				double flux = values[term];
				scalingFactors[term] = flux / source.TotalFlux(frequency, Polarization::StokesI);
			}
			for(ModelComponent& component : source)
			{
				ao::uvector<double> fitValues(scaleToNTerms);
				for(size_t term = 0; term!=scaleToNTerms; ++term)
				{
					double frequency = frequencies[term];
					double oldFlux = component.SED().FluxAtFrequency(frequency, Polarization::StokesI);
					double newFlux = oldFlux * scalingFactors[term];
					fitValues[term] = newFlux;
				}
				SpectralFittingMode mode;
				PowerLawSED& sed = static_cast<PowerLawSED&>(component.SED());
				double refFreq = sed.ReferenceFrequencyHz();
				if(sed.IsLogarithmic())
					mode = LogPolynomialSpectralFitting;
				else
					mode = PolynomialSpectralFitting;
				size_t nTerms = sed.NTerms();
				SpectralFitter fitter(mode, nTerms);
				fitter.SetFrequencies(frequencies.data(), weights.data(), scaleToNTerms);
				ao::uvector<double> newTerms;
				fitter.Fit(newTerms, fitValues.data());
				sed.SetFromStokesIFit(refFreq, newTerms);
			}
		}
	}
	if(shiftRA != 0.0 || shiftDec != 0.0)
	{
		for(ModelSource& source : model)
		{
			for(ModelComponent& component : source)
			{
				component.SetPosRA(component.PosRA() + shiftRA);
				component.SetPosDec(component.PosDec() + shiftDec);
			}
		}
	}
	if(replaceSI)
	{
		for(ModelSource& source : model)
		{
			for(ModelComponent& component : source)
			{
				if(component.HasPowerLawSED())
				{
					PowerLawSED& sed = static_cast<PowerLawSED&>(component.SED());
					long double refFreq = sed.ReferenceFrequencyHz();
					double brightness[4];
					for(size_t p=0; p!=4; ++p)
						brightness[p] = sed.FluxAtFrequencyFromIndex(refFreq, p);
					sed.SetIsLogarithmic(true);
					sed.SetData(refFreq, brightness, replaceSIvalues);
				}
				else throw std::runtime_error("Can only set SI of sources with a power law SED");
			}
		}
	}
	
	if(simUniformN != 0 || simPopulationFluxHi != 0.0)
	{
		size_t simN;
		double simRa, simDec, simDist;
		// According to Franzen et al. (2016) at 154 MHz:
		// dN / dS = 6998 S^(-154) Jy^-1 / Sr^-1
		double popExp = -1.54, popFac = 6998.0;
		bool samplePL;
		if(simUniformN != 0)
		{
			simRa = simUniformRa;
			simDec = simUniformDec;
			simDist = simUniformDist;
			simN = simUniformN;
			samplePL = false;
		}
		else {
			simRa = simPopulationRa;
			simDec = simPopulationDec;
			simDist = simPopulationDist;
			double area = M_PI * simPopulationDist * simPopulationDist;
			// integration of dN / dS over [lo, hi] :
			simN = round(area * popFac / (popExp+1.0) * (pow(simPopulationFluxHi, (popExp+1.0)) - pow(simPopulationFluxLo, (popExp+1.0))));
			std::cout << "N = " << simN << "\n";
			samplePL = true;
		}
		
		std::mt19937 rnd;
		std::uniform_real_distribution<double> dist(0.0, 1.0);
		std::normal_distribution<double> siDist(simSIAverage, simSIStdDev); 
		for(size_t i=0; i!=simN; ++i)
		{
			double ra, dec;
			do {
				ra = dist(rnd) * 2.0*M_PI;
				dec = 0.5*M_PI - acos(2.0*dist(rnd) - 1.0);
			} while(ImageCoordinates::AngularDistance(ra, dec, simRa, simDec) > simDist);
			ModelComponent comp;
			comp.SetPosRA(ra);
			comp.SetPosDec(dec);
			PowerLawSED sed;
			
			double flux;
			if(samplePL)
			{
				double u = dist(rnd);
				double loTerm = pow(simPopulationFluxLo, popExp+1.0);
				flux = pow((pow(simPopulationFluxHi, popExp+1.0) - loTerm) * u + loTerm, 1.0/(popExp+1.0));
				//std::cout << ra << "," << dec << ": " << flux << "Jy\n";
			}
			else {
				flux = simUniformFlux;
				//std::cout << ra << "," << dec << "\n";
			}
			
			double b[] = {flux, 0.0, 0.0, 0.0};
			double spectralIndex = siDist(rnd);
			ao::uvector<double> si(1, spectralIndex);
			sed.SetData(150.0e6, b, si);
			comp.SetSED(sed);
			ModelSource source;
			source.AddComponent(comp);
			std::ostringstream s;
			s << "sim" << (i+1);
			source.SetName(s.str());
			model.AddSource(source);
		}
	}
	
	if(!selectSourceName.empty())
	{
		size_t index = model.FindSourceIndex(selectSourceName);
		ModelSource source = model.Source(index);
		Model newModel;
		newModel.AddSource(source);
		model = newModel;
	}
	
	if(delNaNs)
	{
		for(Model::iterator s=model.begin(); s!=model.end(); ++s)
		{
			for(ModelSource::iterator c=s->begin(); c!=s->end(); ++c)
			{
				if(c->HasMeasuredSED())
					c->MSED().RemoveInvalidMeasurements();
			}
		}
	}
	if(optimize)
		model.Optimize();
	if(doSort)
		model.Sort();
	if(!sortBeamPrefix.empty())
	{
		std::vector<ao::uvector<double>> images(8);
		FitsReader templateReader(sortBeamPrefix+"-xxr.fits");
		readImage(images[0], templateReader, sortBeamPrefix+"-xxr.fits");
		readImage(images[1], templateReader, sortBeamPrefix+"-xxi.fits");
		readImage(images[2], templateReader, sortBeamPrefix+"-xyr.fits"),
		readImage(images[3], templateReader, sortBeamPrefix+"-xyi.fits"),
		readImage(images[4], templateReader, sortBeamPrefix+"-yxr.fits"),
		readImage(images[5], templateReader, sortBeamPrefix+"-yxi.fits"),
		readImage(images[6], templateReader, sortBeamPrefix+"-yyr.fits"),
		readImage(images[7], templateReader, sortBeamPrefix+"-yyi.fits");
		size_t width = templateReader.ImageWidth(), height = templateReader.ImageHeight();
		double refFreq = templateReader.Frequency();
		std::vector<double> appFluxes(model.SourceCount());
		for(size_t i=0; i!=model.SourceCount(); ++i)
		{
			double totalAbsFlux = 0.0;
			for(ModelSource::iterator compPtr = model.Source(i).begin(); compPtr!=model.Source(i).end(); ++compPtr)
			{
				const SpectralEnergyDistribution& sed = compPtr->SED();
				double stokesFlux[4];
				for(size_t p=0; p!=4; ++p)
					stokesFlux[p] = sed.FluxAtFrequencyFromIndex(refFreq, p);
				
				int x, y;
				long double l, m;
				ImageCoordinates::RaDecToLM<long double>(compPtr->PosRA(), compPtr->PosDec(), templateReader.PhaseCentreRA(), templateReader.PhaseCentreDec(), l, m);
				ImageCoordinates::LMToXY<long double>(l, m, templateReader.PixelSizeX(), templateReader.PixelSizeY(), width, height, x, y);
				
				if(x < 0 || y < 0 || x >= int(width) || y >= int(height))
					throw std::runtime_error("Source in model is outside image!");
				std::complex<double> beamValues[4];
				for(size_t p=0; p!=4; ++p)
				{
					beamValues[p] = std::complex<double>(
						images[p*2]  [y*height + x],
						images[p*2+1][y*height + x]);
				}
				std::complex<double> linearFlux[4], tempValues[4];
				Polarization::StokesToLinear(stokesFlux, linearFlux);
				Matrix2x2::ATimesB(tempValues, beamValues, linearFlux);
				Matrix2x2::ATimesHermB(linearFlux, tempValues, beamValues);
				Polarization::LinearToStokes(linearFlux, stokesFlux);
				totalAbsFlux += stokesFlux[0];
			}
			
			ModelSource& s = model.Source(i);
			appFluxes[i] = totalAbsFlux;
			s.SetUserData(&appFluxes[i]);
		}
		
		model.Sort(compareWithBeam);
	}
	
	if(doSearch)
	{
		std::map<double, const ModelSource*> sortedSources;
		for(size_t i=0; i!=model.SourceCount(); ++i)
		{
			const ModelSource& source = model.Source(i);
			double ra = source.MeanRA(), dec = source.MeanDec();
			double dist = ImageCoordinates::AngularDistance(ra, dec, searchRA, searchDec);
			sortedSources.insert(std::pair<double,const ModelSource*>(dist, &source));
		}
		Model newModel;
		size_t count = 0;
		for(std::map<double, const ModelSource*>::const_iterator i=sortedSources.begin(); i!=sortedSources.end() && count < searchCount; ++i, ++count)
		{
			newModel.AddSource(*i->second);
		}
		model = newModel;
	}
	
	if(!toAppBeamPrefix.empty())
	{
		std::vector<ao::uvector<double>> images(8);
		FitsReader templateReader(toAppBeamPrefix+"-xxr.fits");
		readImage(images[0], templateReader, toAppBeamPrefix+"-xxr.fits");
		readImage(images[1], templateReader, toAppBeamPrefix+"-xxi.fits");
		readImage(images[2], templateReader, toAppBeamPrefix+"-xyr.fits"),
		readImage(images[3], templateReader, toAppBeamPrefix+"-xyi.fits"),
		readImage(images[4], templateReader, toAppBeamPrefix+"-yxr.fits"),
		readImage(images[5], templateReader, toAppBeamPrefix+"-yxi.fits"),
		readImage(images[6], templateReader, toAppBeamPrefix+"-yyr.fits"),
		readImage(images[7], templateReader, toAppBeamPrefix+"-yyi.fits");
		size_t width = templateReader.ImageWidth(), height = templateReader.ImageHeight();
		double refFreq = templateReader.Frequency();
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			std::cout << sourcePtr->Name() << ": " << sourcePtr->TotalFlux(refFreq, Polarization::StokesI) << " -> ";
			for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				const SpectralEnergyDistribution& sed = compPtr->SED();
				double stokesFlux[4];
				
				for(size_t p=0; p!=4; ++p)
					stokesFlux[p] = sed.FluxAtFrequencyFromIndex(refFreq, p);
				
				long double l, m;
				ImageCoordinates::RaDecToLM<long double>(compPtr->PosRA(), compPtr->PosDec(), templateReader.PhaseCentreRA(), templateReader.PhaseCentreDec(), l, m);
				
				int x, y;
				ImageCoordinates::LMToXY<long double>(l, m, templateReader.PixelSizeX(), templateReader.PixelSizeY(), width, height, x, y);
				
				if(x < 0 || y < 0 || x >= int(width) || y >= int(height))
					throw std::runtime_error("Source in model is outside image!");
				std::complex<double> beamValues[4];
				for(size_t p=0; p!=4; ++p)
				{
					beamValues[p] = std::complex<double>(
						images[p*2]  [y*height + x],
						images[p*2+1][y*height + x]);
				}
				
				std::complex<double> linearFlux[4], tempValues[4];
				Polarization::StokesToLinear(stokesFlux, linearFlux);
				Matrix2x2::ATimesB(tempValues, beamValues, linearFlux);
				Matrix2x2::ATimesHermB(linearFlux, tempValues, beamValues);
				Polarization::LinearToStokes(linearFlux, stokesFlux);
				PowerLawSED plSED;
				ao::uvector<double> siTerms;
				plSED.SetData(refFreq, stokesFlux, siTerms);
				compPtr->SetSED(plSED);
			}
			std::cout << sourcePtr->TotalFlux(refFreq, Polarization::StokesI) << '\n';
		}
	}
	
	if(applyIndexThreshold)
	{
		size_t index = 0;
		Model temp;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			if(index < thresholdIndex)
			{
				temp.AddSource(*sourcePtr);
			}
			else break;
			++index;
		}
		model = temp;
	}
	
	if(applyThreshold)
	{
		Model temp;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			double totalFlux = sourcePtr->TotalFlux(Polarization::StokesI);
			if(totalFlux >= threshold)
				temp.AddSource(*sourcePtr);
		}
		model = temp;
	}
	
	if(applyComponentThreshold)
	{
		Model temp;
		for(Model::const_iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			ModelSource newSource = *sourcePtr;
			newSource.ClearComponents();
			for(ModelSource::const_iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				if(compPtr->HasSED())
				{
					const SpectralEnergyDistribution& sed = compPtr->SED();
					long double refFreq = sed.ReferenceFrequencyHz();
					long double flux = sed.FluxAtFrequency(refFreq, Polarization::StokesI);
					if(std::isfinite(flux) && flux >= threshold)
						newSource.AddComponent(*compPtr);
				}
			}
			temp.AddSource(newSource);
		}
		model = temp;
	}
	
	if(minSeparation != 0.0)
	{
		Model newModel;
		for(size_t i=0; i!=model.SourceCount(); ++i) {
			const ModelSource &s1 = model.Source(i);
			bool remain = true;
			for(size_t j=0; j!=model.SourceCount(); ++j) {
				if(i != j)
				{
					const ModelSource &s2 = model.Source(j);
					double dist = ImageCoordinates::AngularDistance(s1.MeanRA(), s1.MeanDec(), s2.MeanRA(), s2.MeanDec());
					if(dist < minSeparation)
					{
						if(s1.TotalFlux(Polarization::StokesI) < s2.TotalFlux(Polarization::StokesI))
							remain = false;
					}
				}
			}
			if(remain)
				newModel.AddSource(s1);
		}
		model = newModel;
	}
	
	if(applyClusterThreshold)
	{
		Model temp;
		std::vector<std::string> clusterNames;
		model.GetClusterNames(clusterNames);
		size_t thresholded = 0, total = clusterNames.size();
		for(std::vector<std::string>::const_iterator clName=clusterNames.begin(); clName!=clusterNames.end(); ++clName)
		{
			SourceGroup group;
			model.GetSourcesInCluster(*clName, group);
			
			if(std::fabs(group.TotalFlux(Polarization::StokesI)) >= threshold)
			{
				for(SourceGroup::const_iterator s=group.begin(); s!=group.end(); ++s)
					temp.AddSource(*s);
			}
			else
				++thresholded;
		}
		model = temp;
		std::cout << "Thresholded " << thresholded << " / " << total << " clusters.\n";
	}
	
	if(rndN != 0)
	{
		for(size_t i=0; i!=rndN; ++i)
		{
			double ra, dec;
			do {
				ra = 2.0 * M_PI * (double) rand() / (double) RAND_MAX,
				dec = M_PI * ((double) rand() / (double) RAND_MAX - 0.5);
			} while(ImageCoordinates::AngularDistance(rndRA, rndDec, ra, dec) >= rndDist);
			ModelComponent component;
			component.SetPosRA(ra);
			component.SetPosDec(dec);
			component.SetSED(MeasuredSED(1.0, 150000000.0));
			ModelSource source;
			source.AddComponent(component);
			std::ostringstream name;
			name << "rnd" << i;
			source.SetName(name.str());
			model.AddSource(source);
		}
	}
	if(nearFilter || outsideFilter)
	{
		Model newModel;
		for(size_t i = model.SourceCount(); i>0; --i)
		{
			ModelSource& source = model.Source(i-1);
			bool isNear = false;
			if(source.ComponentCount() > 0)
			{
				double dist = ImageCoordinates::AngularDistance(source.Peak().PosRA(), source.Peak().PosDec(), nearFilterRA, nearFilterDec);
				isNear = (dist <= nearFilterDist);
			}
			if((isNear && nearFilter) || (!isNear && outsideFilter))
			{
				newModel.AddSource(source);
			}
		}
		model = std::move(newModel);
	}
	
	if(uncollect)
	{
		Model newModel;
		for(Model::const_iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			std::string name = sourcePtr->Name();
			ModelSource templateSource(*sourcePtr);
			templateSource.ClearComponents();
			for(size_t i=0; i!=sourcePtr->ComponentCount(); ++i)
			{
				ModelSource newSource(templateSource);
				newSource.AddComponent(sourcePtr->Component(i));
				newSource.SetName(name + std::to_string(i));
				newModel.AddSource(std::move(newSource));
			}
		}
		model = std::move(newModel);
	}
	
	if(doCollect)
	{
		ModelSource newSource;
		newSource.SetName(collectName);
		for(Model::const_iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::const_iterator compPtr = sourcePtr->begin(); compPtr != sourcePtr->end(); ++compPtr)
			{
				newSource.AddComponent(*compPtr);
			}
		}
		model = Model();
		model.AddSource(newSource);
	}
	
	if(!setClusterName.empty())
	{
		for(ModelSource& source : model)
			source.SetClusterName(setClusterName);
	}
	
	if(!splits.empty())
	{
		for(const SplitLine& split : splits)
		{
			Model newModel;
			for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
			{
				double ra = sourcePtr->MeanRA(), dec = sourcePtr->MeanDec();
				if(split.IsPointAboveLine(ra, dec))
				{
					newModel.AddSource(*sourcePtr);
				}
			}
			model = std::move(newModel);
		}
	}
	
	if(powerlaw || resample || resampleByAveraging)
	{
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				const MeasuredSED &sed = compPtr->MSED();
				const long double startFreq = sed.LowestFrequency();
				const long double endFreq = sed.HighestFrequency();
				if(powerlaw)
				{
					long double e, f;
					MeasuredSED newSED;
					Measurement m1, m2;
					m1.SetFrequencyHz(startFreq);
					m2.SetFrequencyHz(endFreq);
					
					for(size_t p=0; p!=4; ++p)
					{
						sed.FitPowerlaw(f, e, Polarization::IndexToStokes(p));
						m1.SetFluxDensityFromIndex(p, f * std::pow(startFreq, e));
						m2.SetFluxDensityFromIndex(p, f * std::pow(endFreq, e));
					}
					newSED.AddMeasurement(m1);
					newSED.AddMeasurement(m2);
					compPtr->SetSED(newSED);
				}
				else if(resample || resampleByAveraging)
				{
					MeasuredSED newSED;
					const long double newBandSize = (endFreq - startFreq) / newChannelCount;
					for(size_t newChIndex=0; newChIndex!=newChannelCount; ++newChIndex)
					{
						long double chStartFreq = startFreq + newBandSize*newChIndex;
						long double chEndFreq = startFreq + newBandSize*(newChIndex+1.0);
						
						Measurement m;
						m.SetFrequencyHz((chStartFreq+chEndFreq)*0.5);
						
						long double flux;
						if(newChannelCount >= sed.MeasurementCount())
							flux = sed.FluxAtFrequency((chStartFreq+chEndFreq)*0.5, Polarization::StokesI);
						else if(resampleByAveraging)
							flux = sed.AverageFlux(chStartFreq, chEndFreq, Polarization::StokesI);
						else 
							flux = sed.IntegratedFlux(chStartFreq, chEndFreq, Polarization::StokesI);
						m.SetFluxDensity(Polarization::StokesI, flux);
						
						PolarizationEnum pols[3] = { Polarization::StokesQ, Polarization::StokesU, Polarization::StokesV };
						for(size_t p=0; p!=3; ++p)
						{
							if(newChannelCount >= sed.MeasurementCount())
								flux = sed.FluxAtFrequency((chStartFreq+chEndFreq)*0.5, pols[p]);
							else
								flux = sed.AverageFlux(chStartFreq, chEndFreq, pols[p]);
							m.SetFluxDensity(pols[p], flux);
						}
						
						newSED.AddMeasurement(m);
					}
					
					compPtr->SetSED(newSED);
				}
			}
		}
	}
	
	for(size_t p=0; p!=4; ++p)
	{
		if(setPolarization[p])
		{
			for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
			{
				for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
				{
					MeasuredSED &sed = compPtr->MSED();
					for(MeasuredSED::iterator m=sed.begin(); m!=sed.end(); ++m)
					{
						m->second.SetFluxDensityFromIndex(p, setPolFlux[p]);
					}
				}
			}
		}
	}
	if(scale != 1.0)
	{
		for(ModelSource& source : model)
		{
			for(ModelComponent& comp : source)
			{
				comp.SED() *= scale;
			}
		}
	}
	if(scaleToFreq != 0.0)
	{
		std::cout << "Source\tScale factor\n";
		for(ModelSource& source : model)
		{
			double flux = source.TotalFlux(scaleToFreq, Polarization::StokesI);
			double factor = scaleToFlux / flux;
			std::cout << source.Name() << '\t' << std::setprecision(15) << factor << '\n';
			for(ModelComponent& comp : source)
			{
				comp.SED() *= factor;
			}
		}
	}
	
	if(scalePeak || scaleSource)
	{
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			long double factorA[4], factorB[4];
			for(size_t p=0; p!=4; ++p)
			{
				if(scalePeak)
				{
					long double
						oldFluxA = sourcePtr->Peak().SED().FluxAtFrequency(scaleFreqA, Polarization::IndexToStokes(p)),
						oldFluxB = sourcePtr->Peak().SED().FluxAtFrequency(scaleFreqB, Polarization::IndexToStokes(p));
					factorA[p] = oldFluxA==0.0 ? 0.0 : scalePeakA / oldFluxA;
					factorB[p] = oldFluxB==0.0 ? 0.0 : scalePeakB / oldFluxB;
				} else {
					long double
						oldFluxA = sourcePtr->TotalFlux(scaleFreqA, Polarization::IndexToStokes(p)),
						oldFluxB = sourcePtr->TotalFlux(scaleFreqB, Polarization::IndexToStokes(p));
					factorA[p] = oldFluxA==0.0 ? 0.0 : scalePeakA / oldFluxA;
					factorB[p] = oldFluxB==0.0 ? 0.0 : scalePeakB / oldFluxB;
				}
			}
			std::cout << "Scale factor for " << sourcePtr->Name() << ": " << factorA[0];
			for(size_t p=1; p!=4; ++p) std::cout << ',' << factorA[p];
			std::cout << " - " << factorB[0];
			for(size_t p=1; p!=4; ++p) std::cout << ',' << factorB[p];
			std::cout << '\n';
			if(scaleWithSI)
			{
				NonLinearPowerLawFitter fitter;
				fitter.AddDataPoint(scaleFreqA, scalePeakA);
				fitter.AddDataPoint(scaleFreqB, scalePeakB);
				double e, f;
				fitter.FastFit(e, f);
				double si = e;
				std::cout << "Spectral index=" << si << '\n';
				ao::uvector<double> sis(1, si);
				double freq = 0.5*(scaleFreqA+scaleFreqB);
				for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
				{
					double brightnessVec[4];
					for(size_t p=0; p!=4; ++p)
					{
						long double oldFluxA = compPtr->SED().FluxAtFrequency(scaleFreqA, Polarization::IndexToStokes(p));
						long double oldFluxB = compPtr->SED().FluxAtFrequency(scaleFreqB, Polarization::IndexToStokes(p));
						NonLinearPowerLawFitter fitter2;
						fitter2.AddDataPoint(scaleFreqA, oldFluxA*factorA[p]);
						fitter2.AddDataPoint(scaleFreqB, oldFluxB*factorB[p]);
						fitter2.FastFit(e, f);
						brightnessVec[p] = NonLinearPowerLawFitter::Evaluate(f, e, freq);
					}
					PowerLawSED sed;
					sed.SetData(freq, brightnessVec, sis);
					compPtr->SetSED(sed);
				}
			}
			else {
				for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
				{
					Measurement measA, measB;
					measA.SetFrequencyHz(scaleFreqA);
					measB.SetFrequencyHz(scaleFreqB);
					for(size_t p=0; p!=4; ++p)
					{
						long double oldFluxA = compPtr->SED().FluxAtFrequency(scaleFreqA, Polarization::IndexToStokes(p));
						long double oldFluxB = compPtr->SED().FluxAtFrequency(scaleFreqB, Polarization::IndexToStokes(p));
						measA.SetFluxDensityFromIndex(p, oldFluxA*factorA[p]);
						measB.SetFluxDensityFromIndex(p, oldFluxB*factorB[p]);
					}
					MeasuredSED sed;
					sed.AddMeasurement(measA);
					sed.AddMeasurement(measB);
					compPtr->SetSED(sed);
				}
			}
		}
	}
		
	if(unpolarized)
	{
		model.SetUnpolarized();
	}
	
	if(toPowerlaw != 0)
	{
		Model newModel;
		for(const ModelSource& source : model)
		{
			ModelSource plSource(source);
			plSource.ClearComponents();
			for(const ModelComponent& comp : source)
			{
				if(comp.HasMeasuredSED())
				{
					const MeasuredSED& sed = comp.MSED();
					double freq = sed.CentreFrequency();
					ao::uvector<double> terms;
					sed.FitLogPolynomial(terms, toPowerlaw, Polarization::StokesI, freq);
					
					PowerLawSED plSed;
					plSed.SetFromStokesIFit(freq, terms);
					
					ModelComponent plComp(comp);
					plComp.SetSED(plSed);
					plSource.AddComponent(std::move(plComp));
				}
				else {
					plSource.AddComponent(comp);
				}
			}
			newModel.AddSource(std::move(plSource));
		}
		model = std::move(newModel);
	}
	
	if(minMeasurements != 0)
	{
		Model newModel;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			bool isValid = true;
			size_t measCount = 0;
			for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				if(!compPtr->HasMeasuredSED())
					throw std::runtime_error("min-measurements requires measured SEDs");
				MeasuredSED& sed = compPtr->MSED();
				if(sed.MeasurementCount() < minMeasurements)
				{
					isValid = false;
					measCount = sed.MeasurementCount();
				}
			}
			if(isValid)
				newModel.AddSource(*sourcePtr);
			else {
				std::cout << "Removing source " << sourcePtr->Name() << ": a component has only " << measCount << " measurements.\n";
			}
		}
		model = std::move(newModel);
	}
	
	if(setFrequency)
	{
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				if(!compPtr->HasMeasuredSED())
					throw std::runtime_error("setfrequency requires measured SEDs");
				MeasuredSED& sed = compPtr->MSED();
				if(sed.MeasurementCount()!=1)
					throw std::runtime_error("setfrequency requires one measurement per component");
				Measurement m = sed.begin()->second;
				m.SetFrequencyHz(setFrequencyValue*1e6);
				MeasuredSED newSED;
				newSED.AddMeasurement(m);
				compPtr->SetSED(newSED);
			}
		}
	}
	
	if(sourceStats)
	{
		const size_t n = model.SourceCount();
		double sum = 0.0, diffSq = 0.0;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			const ModelSource& source = *sourcePtr;
			double stokesI = source.TotalFlux(150000000.0, Polarization::StokesI);
			sum += stokesI;
			diffSq += (stokesI - 1.0) * (stokesI - 1.0);
		}
		double mean = sum / n, stdDevSum = 0.0;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			const ModelSource& source = *sourcePtr;
			double stokesI = source.TotalFlux(150000000.0, Polarization::StokesI);
			stdDevSum += (stokesI - mean) * (stokesI - mean);
		}
		std::cout <<
			"Source statistics:\n"
			" Mean: " << mean << "\n"
			" Standard deviation: " << sqrt(stdDevSum/n) << "\n"
			" Deviation from 1: " << sqrt(diffSq/n) << '\n';
	}
	
	if(logNlogS)
	{
		LogHistogram histogram;
		for(Model::iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			const ModelSource& source = *sourcePtr;
			double v = source.TotalFlux(logNlogSFrequency*1000000.0, Polarization::StokesI);
			histogram.Add(v);
		}
		long double max = histogram.MaxAmplitude(), min = histogram.MinPositiveAmplitude();
		long double step = powl(max/min, 1.0/logNlogSBinCount);
		std::ofstream plotDataStream("lognlogs-data.txt");
		ao::uvector<long double> vals(logNlogSBinCount);
		for(size_t i=0; i!=logNlogSBinCount; ++i)
		{
			long double
				intervalStart = min * powl(step, i),
				mid = min * powl(step, double(i) + 0.5),
				intervalEnd = min * powl(step, i+1),
				val = histogram.NormalizedCount(intervalStart, intervalEnd),
				valNormalized = log10(val * powl(mid, 2.5));
			vals[i] = valNormalized;
			plotDataStream << intervalStart << '\t' << intervalEnd << '\t' << intervalStart * sqrt(step) << '\t' << val << '\t' << valNormalized << '\n';
			intervalStart = intervalEnd;
		}
		size_t avgCount = 0;
		long double avg = 0.0;
		for(size_t i=logNlogSBinCount/5*3; i<logNlogSBinCount; ++i)
		{
			if(std::isfinite(vals[i]))
			{
				avg += vals[i];
				++avgCount;
			}
		}
		avg /= avgCount;
		
		std::ofstream plotStream("lognlogs.plt");
		plotStream <<
			"set terminal postscript enhanced color\n"
			"set logscale x\n"
			"set xrange [" << min << ":" << max << "]\n"
			"set yrange [:1.2]\n"
			"set output \"lognlogs.ps\"\n"
			"set key bottom right\n"
			"set xlabel \"Flux (Jy)\"\n"
			"set ylabel \"log N/S^{-2.5}/avg\"\n"
			"plot \"lognlogs-data.txt\" using 2:(column(5)/" << avg << ") with lines lw 2.0 title \"Measured\", \\\n"
			"1 with lines lw 2.0 lt 3 title \"Euclidean\"\n";
	}
	
	if(powerSumLimit != 0.0)
	{
		double sum = 0.0, powerSum = 0.0;
		Model sorted(model);
		sorted.Sort();
		bool limitReached = false;
		size_t count = 0;
		for(const ModelSource& s : sorted)
		{
			++count;
			const SpectralEnergyDistribution& sed = s.Component(0).SED();
			long double f = sed.FluxAtFrequency(sed.ReferenceFrequencyHz(), Polarization::StokesI);
			sum += f;
			powerSum += f*f;
			if(!limitReached && powerSum >= powerSumLimit)
			{
				limitReached = true;
				std::cout << count << " sources required to reach power sum limit of " << powerSumLimit << " Jy^2.\n";
			}
		}
		std::cout << "Total flux sum = " << sum << " Jy.\n";
		std::cout << "Total power sum = " << powerSum << " Jy^2.\n";
	}
	
	if(!dpppModelFilename.empty())
	{
		NDPPP::SaveSkyModel(dpppModelFilename, model, true);
	}
	
	if(!kvisAnnFile.empty())
	{
		std::ofstream str(kvisAnnFile);
		for(const ModelSource& source : model)
		{
			str << "TEXT "
				<< source.MeanRA()*180.0/M_PI
				<< " " << source.MeanDec()*180.0/M_PI
				<< " " << source.Name() << "\n";
		}
	}
	
	if(!sagecalPrefix.empty())
	{
		std::ofstream sagecalSourceFile(sagecalPrefix + "_model.txt");
		std::map<std::string, std::vector<const ModelSource*>> clusters;
		sagecalSourceFile << "## name h m s d m s I Q U V si0 si1 si2 RM eX(rad) exY(rad) PA(rad) freq0\n";
		for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
		{
			clusters[s->ClusterName()].push_back(&*s);

			for(size_t compIndex=0; compIndex!=s->ComponentCount(); ++compIndex)
			{
				const ModelComponent& c = s->Component(compIndex);
				const PowerLawSED* pl = dynamic_cast<const PowerLawSED*>(&c.SED());
				if(pl == 0)
					throw std::runtime_error("Expecting model with power law SEDs");
				long double ra = c.PosRA(), dec = c.PosDec();
				// I'm not sure this is necessary for Sagecal, but just for
				// sure convert negative ra's to positive ones.
				if(ra < 0.0)
					ra += 2.0 * M_PI;
				int raH, raM, decD, decM;
				double raS, decS;
				RaDecCoord::RAToHMS(ra, raH, raM, raS);
				RaDecCoord::DecToDMS(dec, decD, decM, decS);
				std::string cName = getSagecalSourceName(*s, c, compIndex);
				// Note that the sagecal axes are in radians, and they're not
				// the length of the axes, but half of that.
				// Also, polarization angle is clockwise from North, while normal PA is
				// anti-clockwise from East.
				long double
					axMaj = c.MajorAxis()*0.5,
					axMin = c.MinorAxis()*0.5,
					axPA = -0.5*M_PI + c.PositionAngle();
				double refFreq, iquv[4];
				std::vector<double> si;
				pl->GetData(refFreq, iquv, si);
				if(si.size() > 3)
					throw std::runtime_error("Expecting 1-3 SI terms in model");
				while(si.size() < 3)
					si.push_back(0.0);
				
				sagecalSourceFile
					<< cName << ' '
					<< raH << ' ' << raM << ' ' << raS << ' '
					<< decD << ' ' << decM << ' ' << decS << ' '
					<< iquv[0] << ' ' << iquv[1] << ' '
					<< iquv[2] << ' ' << iquv[3] << ' '
					<< si[0] << ' '
					<< si[1]/M_LN10 << ' '
					<< si[2]/(M_LN10*M_LN10) << ' '
					<< 0.0 << ' ' // rm
					<< axMaj << ' ' << axMin << ' ' << axPA << ' '
					<< refFreq << '\n';
			}
		}
		std::ofstream sagecalClusterFile(sagecalPrefix + "_clustering.txt");
		sagecalClusterFile << "## cID chunk_size source1 source2 ...\n";
		size_t clusterId = 0;
		bool isFirst = true;
		// We need to traverse in original cluster order
		std::set<std::string> clustersAdded;
		for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
		{
			std::string cName = s->ClusterName();
			if(clustersAdded.count(cName) == 0)
			{
				clustersAdded.insert(cName);
				if(cName.empty())
					sagecalClusterFile << clusters.size()+1;
				else
					sagecalClusterFile << clusterId;
				if(isFirst)
				{
					sagecalClusterFile << ' ' << sagecalFirstClusterChunkSize;
					isFirst = false;
				}
				else
					sagecalClusterFile << " 1";
				const std::vector<const ModelSource*>& sources = clusters[cName];
				for(std::vector<const ModelSource*>::const_iterator s=sources.begin(); s!=sources.end(); ++s)
				{
					for(size_t compIndex=0; compIndex!=(*s)->ComponentCount(); ++compIndex)
					{
						const ModelComponent& c = (*s)->Component(compIndex);
						sagecalClusterFile << ' ' << getSagecalSourceName(**s, c, compIndex);
					}
				}
				sagecalClusterFile << '\n';
				++clusterId;
			}
		}
	}
	
	if(outputCsv)
	{
		std::ofstream csvFile(csvFilename);
		csvFile << "freq,i,q,u,v\n";
		if(model.SourceCount() != 0)
		{
			if(model.SourceCount() != 1)
				std::cout << "Warning: multiple sources in model, but will only output a csv file for the first source.\n";
			const ModelSource& source = *model.begin();
			if(source.ComponentCount() != 1)
				std::cout << "Warning: first source has multiple components; only outputting first component.\n";
			if(source.begin()->HasMeasuredSED())
			{
				const MeasuredSED &sed = source.begin()->MSED();

				std::vector<Measurement> measurements;
				sed.GetMeasurements(measurements);
				
				for(MeasuredSED::const_iterator iter=sed.begin(); iter!=sed.end(); ++iter)
				{
					const Measurement& m = iter->second;
					long double
						i = m.FluxDensity(Polarization::StokesI), q = m.FluxDensity(Polarization::StokesQ),
						u = m.FluxDensity(Polarization::StokesU), v = m.FluxDensity(Polarization::StokesV);
					if(std::isfinite(i) && std::isfinite(q) && std::isfinite(u) && std::isfinite(v))
					{
						csvFile
							<< m.FrequencyHz()*1e-6 << ','
							<< i << ',' << q << ',' << u << ',' << v << '\n';
					}
				}
			}
		}
	}
	
	if(outputSICsv)
	{
		std::ofstream csvFile(csvFilename);
		csvFile << "RA,Dec,Brightness,SI\n";
		if(model.SourceCount() != 0)
		{
			for(Model::const_iterator s=model.begin(); s!=model.end(); ++s)
			{
				const ModelSource& source = *s;
				if(source.begin()->HasMeasuredSED())
				{
					const MeasuredSED &sed = source.GetIntegratedMSED();
					double flux = sed.AverageFlux(Polarization::StokesI);
					long double f, si;
					sed.FitPowerlaw(f, si, Polarization::StokesI);

					csvFile
						<< source.MeanRA()*180.0/M_PI << ','
						<< source.MeanDec()*180.0/M_PI << ','
						<< flux << ',' << si << '\n';
				}
			}
		}
	}
	
	if(outputRMPlot)
	{
		std::ofstream plotStream("rmplot.plt");
		plotStream <<
			"set terminal postscript enhanced color\n"
			"#set logscale y\n"
			"#set xrange [0.001:]\n"
			"#set yrange [-8:2]\n"
			"set output \"rmplot.ps\"\n"
			"#set key bottom left\n"
			"set xlabel \"rad/m^2\"\n"
			"set ylabel \"Flux (Jy/beam)\"\n"
			"plot \\\n";
			
		size_t compIndex = 0;
		for(Model::const_iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::const_iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				std::ostringstream dataStreamName;
				dataStreamName << "rm" << compIndex << ".txt";
				std::ofstream dataStream(dataStreamName.str().c_str());
				
				plotStream << "\"" << dataStreamName.str() << "\" using 1:2 with lines lw 1.0 title \"\"\\\n";
				
				if(compIndex != model.ComponentCount()-1)
				{
					plotStream << ",\\";
				}
				plotStream << "\n";
				
				const MeasuredSED &sed = compPtr->MSED();
				RMSynthesis rmSynth(sed);
				rmSynth.Synthesize();
				const ao::uvector<double>& fdf = rmSynth.FDF();
				
				for(size_t i=0; i!=fdf.size(); ++i)
				{
					dataStream
						<< rmSynth.IndexToValue(i) << '\t'
						<< std::abs(fdf[i]) << '\n';
				}
				++compIndex;
			}
		}
	}
	
	if(outputPlot)
	{
		if(!model.Source(0).front().HasMeasuredSED())
			throw std::runtime_error("Plotting works only with models with measured SEDs");
		const MeasuredSED& first = model.Source(0).front().MSED();
		double freqStart = first.LowestFrequency()*1e-6, freqEnd = first.HighestFrequency()*1e-6;
		
		GNUPlot plot("spectrum", "Frequency (MHz)", "Flux (Jy)", true, true);
		plot.SetXRange(freqStart, freqEnd);
		if(!plotTitle.empty())
			plot.SetTitle(plotTitle);

		GNUPlot plotI("spectrum-I", "Frequency (MHz)", "Flux (Jy)", true, true);
		plotI.SetXRange(freqStart, freqEnd);
		plotI.SetKeyBottomLeft();
		if(!plotTitle.empty())
			plotI.SetTitle(plotTitle);
			
		GNUPlot plotIFit("spectrum-I-fit", "Frequency (MHz)", "Flux (Jy)", true, true);
		plotIFit.SetXRange(freqStart, freqEnd);
		plotIFit.SetKeyBottomLeft();
		if(!plotTitle.empty())
			plotIFit.SetTitle(plotTitle);
		
		GNUPlot plot10I("spectrum10-I", "Frequency (MHz)", "Flux (Jy)", true, true);
		plot10I.SetXRange(freqStart, freqEnd);
		plot10I.SetKeyBottomLeft();
		if(!plotTitle.empty())
			plot10I.SetTitle(plotTitle);

		GNUPlot plotFT("spectrumFT", "k\u2225 (h/Mpc)", "Power (Jy^2)", true, true);
		GNUPlot plotFTdelay("spectrumFTdelay", "delay ({/Symbol m}sec)", "Power (Jy^2)", true, true);
		const ModelSource& firstSource = *model.begin();
		const ModelComponent& firstComponent = *firstSource.begin();
		double maxX, maxDelay;
		if(firstComponent.HasMeasuredSED())
		{
			const MeasuredSED& firstSED = firstComponent.MSED();
			SpectrumFT sft(firstSED.MeasurementCount());
			
			maxX = SpectrumFT::GetMaxKParallell(firstSED);
			plotFT.SetXRange(maxX/firstSED.MeasurementCount(), maxX);
			sft.PreparePSPlot(plotFT, firstSED.HighestFrequency(), false);
			
			maxDelay = SpectrumFT::GetMaxDelayInMuSec(firstSED);
			plotFTdelay.SetXRange(maxDelay/firstSED.MeasurementCount(), maxDelay);
			sft.PreparePSPlot(plotFTdelay, firstSED.HighestFrequency(), true);
		} else
		{
			maxX = 1.0;
			maxDelay = 1.0;
		}
		if(!plotTitle.empty())
			plotFT.SetTitle(plotTitle);

		size_t compIndex = 0, sourceIndex = 0;
		for(Model::const_iterator sourcePtr = model.begin(); sourcePtr!=model.end(); ++sourcePtr)
		{
			for(ModelSource::const_iterator compPtr = sourcePtr->begin(); compPtr!=sourcePtr->end(); ++compPtr)
			{
				std::ostringstream dataStreamName;
				dataStreamName << "spectrum" << compIndex;
				
				GNUPlot::Line* lineI = plot.AddLine(dataStreamName.str()+"I.txt", "I");
				GNUPlot::Line* lineQ = plot.AddLine(dataStreamName.str()+"Q.txt", "Q");
				GNUPlot::Line* lineU = plot.AddLine(dataStreamName.str()+"U.txt", "U");
				GNUPlot::Line* lineV = plot.AddLine(dataStreamName.str()+"V.txt", "V");
				
				GNUPlot::Line* lineFT = 0;
				GNUPlot::Line* lineFTdelay = 0;
				if(doFT) lineFT = plotFT.AddLine(dataStreamName.str()+"FT.txt", "");
				if(doFT) lineFTdelay = plotFTdelay.AddLine(dataStreamName.str()+"FT-delay.txt", "");
				
				plotI.AddLineFromExistingFile(dataStreamName.str()+"I.txt", "");
				plotIFit.AddLineFromExistingFile(dataStreamName.str()+"I.txt", "");
				if(sourceIndex < 10 && compPtr == sourcePtr->begin())
				{
					plot10I.AddLineFromExistingFile(dataStreamName.str()+"I.txt", sourcePtr->Name());
				}
				
				if(compPtr->HasMeasuredSED())
				{
					const MeasuredSED &sed = compPtr->MSED();
					
					long double e, f;
					sed.FitPowerlaw(f, e, Polarization::StokesI);
					std::ostringstream funcStr;
					funcStr << f << " * (x*1000000)**" << e;
					plotIFit.AddFunction(funcStr.str());
					ao::uvector<double> terms;
					sed.FitLogPolynomial(terms, 3, Polarization::StokesI, sed.ReferenceFrequencyHz());
					std::ostringstream func2Str;
					func2Str << terms[0] << "* exp(" << terms[1] << "*log(x/" << sed.ReferenceFrequencyHz() << "*1000000) + " << terms[2] << "*log(x/" << sed.ReferenceFrequencyHz() << "*1000000)**2)";
					plotIFit.AddFunction(func2Str.str()); 
					
					if(doFT)
					{
						SpectrumFT sft(sed.MeasurementCount());
						ao::uvector<double> power;
						sft.GetFTPower(power, sed, 3, false, true);
						for(size_t i=0; i!=power.size(); ++i)
						{
							lineFT->AddPoint((i*(maxX)/power.size()), power[i]);
							lineFTdelay->AddPoint((i*(maxDelay)/power.size()), power[i]);
						}
					}
					
					std::vector<Measurement> measurements;
					sed.GetMeasurements(measurements);
					
					for(std::vector<Measurement>::const_iterator i=measurements.begin(); i!=measurements.end(); ++i)
					{
						lineI->AddPoint(i->FrequencyHz()/1000000.0, i->FluxDensity(Polarization::StokesI));
						lineQ->AddPoint(i->FrequencyHz()/1000000.0, i->FluxDensity(Polarization::StokesQ));
						lineU->AddPoint(i->FrequencyHz()/1000000.0, i->FluxDensity(Polarization::StokesU));
						lineV->AddPoint(i->FrequencyHz()/1000000.0, i->FluxDensity(Polarization::StokesV));
					}
					++compIndex;
				}
			}
			++sourceIndex;
		}
	}
	
	if(!saveClusters.empty())
	{
		std::vector<std::string> clusters;
		model.GetClusterNames(clusters);
		Delaunay delaunay;
		for(const std::string& cluster : clusters)
		{
			SourceGroup sourceGroup;
			model.GetSourcesInCluster(cluster, sourceGroup);
			double ra = sourceGroup.MeanRA(), dec = sourceGroup.MeanDec();
			if(!std::isfinite(ra) || !std::isfinite(dec))
				throw std::runtime_error("One of the directions has a non-finite RA or Dec");
			delaunay.AddVertex(ra, dec, nullptr);
		}
		delaunay.Triangulate();
		delaunay.SaveTriangulationAsKvis(saveClusters);
	}
	
	if(outputList) {
		for(Model::const_iterator s=model.begin();
				s!=model.end(); ++s)
		{
			const ModelSource& source = *s;
			std::cout
				<< source.Name() << " "
				<< RaDecCoord::RAToString(source.MeanRA()) << " "
				<< RaDecCoord::DecToString(source.MeanDec()) << " "
				<< source.TotalFlux(Polarization::StokesI) << '\n';
		}
	}
	
	if(!rtsFilename.empty())
	{
		std::ofstream rtsStream(rtsFilename);
		for(const ModelSource& s : model)
		{
			std::string name = s.Name();
			std::replace(name.begin(), name.end(), ' ', '_');
			std::replace(name.begin(), name.end(), '-', 'm');
			ModelSource::const_iterator compIter = s.begin();
			rtsStream << "SOURCE "
				<< std::setprecision(std::numeric_limits<double>::digits10 + 1)
				<< name << " "
				<< compIter->PosRA()/M_PI*12.0 << " "
				<< compIter->PosDec()/M_PI*180.0 << '\n';
			if(compIter->Type() == ModelComponent::GaussianSource)
			{
				double
					rtsMaj = sqrt(M_PI*M_PI / (2.0*M_LN2)) * compIter->MajorAxis() * (180.0*60.0/M_PI),
					rtsMin = sqrt(M_PI*M_PI / (2.0*M_LN2)) * compIter->MinorAxis() * (180.0*60.0/M_PI),
					rtsRot = compIter->PositionAngle()*(180.0/M_PI);
				rtsStream <<
					"GAUSSIAN " << rtsRot << ' ' << rtsMaj << ' ' << rtsMin << '\n';
			}
			rtsStream
				<< "FREQ 150.0e6 "
				<< compIter->SED().FluxAtFrequency(150e6, Polarization::StokesI) << " 0 0 0\n"
				<< "FREQ 200.0e6 "
				<< compIter->SED().FluxAtFrequency(200e6, Polarization::StokesI) << " 0 0 0\n";
			++compIter;
			while(compIter!=s.end())
			{
				rtsStream << "COMPONENT "
					<< compIter->PosRA()/M_PI*12.0 << " "
					<< compIter->PosDec()/M_PI*180.0 << '\n';
				if(compIter->Type() == ModelComponent::GaussianSource)
				{
					double
						rtsMaj = sqrt(M_PI*M_PI / (2.0*M_LN2)) * compIter->MajorAxis() * (180.0*60.0/M_PI),
						rtsMin = sqrt(M_PI*M_PI / (2.0*M_LN2)) * compIter->MinorAxis() * (180.0*60.0/M_PI),
						rtsRot = compIter->PositionAngle()*(180.0/M_PI);
					rtsStream <<
						"GAUSSIAN " << rtsRot << ' ' << rtsMaj << ' ' << rtsMin << '\n';
				}
				rtsStream
					<< "FREQ 150.0e6 "
					<< compIter->SED().FluxAtFrequency(150e6, Polarization::StokesI) << " 0 0 0\n"
					<< "FREQ 200.0e6 "
					<< compIter->SED().FluxAtFrequency(200e6, Polarization::StokesI) << " 0 0 0\n"
					"ENDCOMPONENT\n";
				++compIter;
			}
			rtsStream << "ENDSOURCE\n";
		}
	}
	
	if(evaluateFrequency != 0.0) {
		for(Model::const_iterator s=model.begin();
				s!=model.end(); ++s)
		{
			const ModelSource& source = *s;
			std::cout
				<< source.Name() << " "
				<< RaDecCoord::RAToString(source.MeanRA()) << " "
				<< RaDecCoord::DecToString(source.MeanDec()) << " " << std::setprecision(15)
				<< source.TotalFlux(evaluateFrequency, Polarization::StokesI) << '\n';
		}
	}
	
	if(!outputModel.empty()) {
		std::cout << "Writing model with " << model.SourceCount() << " sources, " << model.ComponentCount() << " components.\n";
		model.Save(outputModel.c_str());
	}
	
	return 0;
}
