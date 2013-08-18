#include "predicter.h"
#include "model.h"
#include "imagecoordinates.h"
#include "beamevaluator.h"
#include "solutionfile.h"

void Predicter::Initialize(ModelSource& source, BeamEvaluator *beamEvaluator)
{
	_beamEvaluator = beamEvaluator;
	
	SourceParameters *parameters = new SourceParameters();
	NumType l, m;
	ImageCoordinates::RaDecToLM<NumType>(source.PosRA(), source.PosDec(), _ra0, _dec0, l, m);
	parameters->l = l;
	parameters->m = m;
	parameters->brightness = new NumType[_channelCount*4];
	parameters->beamValues = new CNumType[_channelCount*4];
	for(size_t ch=0;ch!=_channelCount;++ch)
	{
		for(size_t p=0; p!=4; ++p)
		{
			parameters->brightness[ch*4+p] =
				source.SED().FluxAtChannel(ch, _channelCount, _startFrequency, _endFrequency, p);
		}
		if(beamEvaluator != 0)
		{
			double centreFreq = _startFrequency + (long double) ch * (_endFrequency - _startFrequency) / (long double) (_channelCount-1);
			beamEvaluator->EvaluateAbsToApparentGain(source.PosRA(), source.PosDec(), centreFreq, &parameters->beamValues[ch*4]);
		}
		else {
			parameters->beamValues[ch*4+0] = 1.0; parameters->beamValues[ch*4+1] = 0.0;
			parameters->beamValues[ch*4+2] = 0.0; parameters->beamValues[ch*4+3] = 1.0;
		}
		CNumType temp[4], brightness[4];
		brightness[0] = parameters->brightness[ch*4+0]; brightness[1] = parameters->brightness[ch*4+2];
		brightness[2] = parameters->brightness[ch*4+1]; brightness[3] = parameters->brightness[ch*4+3];
		Matrix2x2::ATimesB(temp, &parameters->beamValues[ch*4], brightness);
		Matrix2x2::PlusATimesHermB(_totalFlux, temp, &parameters->beamValues[ch*4]);
	}
	parameters->lmsqrt = sqrt(1.0 - l*l - m*m);
	
	source.SetUserData(parameters);
}

void Predicter::Initialize(Model& model, const std::string& solutionFile, BeamEvaluator *beamEvaluator)
{
	_beamEvaluator = beamEvaluator;
	
	for(Model::iterator i=model.begin(); i!=model.end(); ++i)
		Initialize(*i, beamEvaluator);
	if(!solutionFile.empty())
		readSolutions(solutionFile);
}

void Predicter::ReportSources(Model& model)
{
	std::cout << "Model predicter initialized with " << model.SourceCount() << " sources of apparent brightness [";
	std::cout << (_totalFlux[0].real() / _channelCount);
	for(size_t p=1; p!=4; ++p)
		std::cout << ',' << (_totalFlux[p].real() / _channelCount);
	std::cout << "]\n";
	
	std::cout << "(absolute brightness: [";
	std::cout << model.TotalFlux(_startFrequency, _endFrequency, 0);
	for(size_t p=1; p!=4; ++p)
		std::cout << ',' << model.TotalFlux(_startFrequency, _endFrequency, p);
	std::cout << "], app at low freq: " << model.TotalFlux(_startFrequency, 0) << ", app at high freq: " << model.TotalFlux(_endFrequency, 0) << ")\n";
}

void Predicter::predict4(CNumType *dest, const ModelSource& source, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2)
{
	switch(source.Type())
	{
		case ModelSource::PointSource:
		{
			SourceParameters *parameters = reinterpret_cast<SourceParameters *>(source.UserData());
			NumType l = parameters->l, m = parameters->m;
			NumType lmsqrt = parameters->lmsqrt;
			NumType angle = 2.0*M_PI*(u*l + v*m + w*(lmsqrt-1.0));
			double sinangle, cosangle;
			sincos(angle, &sinangle, &cosangle);
			for(size_t p=0; p!=4; ++p)
			{
				NumType fact = parameters->brightness[channelIndex*4+p] / lmsqrt;
				dest[p] = CNumType(fact * cosangle, fact * sinangle);
			}
		}
	}
	
	CNumType temp[4];
	if(!_rhsSolutions.empty())
	{
		std::complex<double>
			*antenna1Sol = &_rhsSolutions[a1*_channelCount*4],
			*antenna2Sol = &_rhsSolutions[a2*_channelCount*4];
		Matrix2x2::ATimesB(temp, antenna1Sol, dest);
		Matrix2x2::ATimesHermB(dest, temp, antenna2Sol);
	}
	if(_beamEvaluator != 0)
	{
		SourceParameters *parameters = reinterpret_cast<SourceParameters *>(source.UserData());
		Matrix2x2::ATimesB(temp, &parameters->beamValues[channelIndex*4], dest);
		Matrix2x2::ATimesHermB(dest, temp, &parameters->beamValues[channelIndex*4]);
	}
}

void Predicter::Predict4(Predicter::CNumType* dest, const ModelSource& source, Predicter::NumType u, Predicter::NumType v, Predicter::NumType w, size_t channelIndex, size_t a1, size_t a2)
{
	predict4(dest, source, u, v, w, channelIndex, a1, a2);
}

void Predicter::Predict4(CNumType *dest, const Model& model, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2)
{
	for(size_t p=0; p!=4; ++p)
		dest[p] = 0.0;
	for(Model::const_iterator i=model.begin(); i!=model.end(); ++i)
	{
		CNumType temp[4];
		predict4(temp, *i, u, v, w, channelIndex, a1, a2);
		for(size_t p=0; p!=4; ++p)
			dest[p] += temp[p];
	}
}

void Predicter::readSolutions(const std::string& solutionFile)
{
	SolutionFile file;
	file.OpenForReading(solutionFile.c_str());
	if(file.PolarizationCount() != 4) throw std::runtime_error("Polarization counts in solution file do not match");
	if(_channelCount != file.ChannelCount()) throw std::runtime_error("Channel counts in solution file do not match");
	size_t antennaCount = file.AntennaCount();
	
	_rhsSolutions.resize(antennaCount*_channelCount*4);
	for(size_t a = 0; a!=antennaCount; ++a) {
		std::complex<double> *antennaSol = &_rhsSolutions[a*_channelCount*4];
		for(size_t ch = 0; ch!=_channelCount; ++ch) {
			for(size_t p = 0; p!=4; ++p) {
				antennaSol[ch*4+p] = file.ReadNextSolution();
			}
		}
	}
}
