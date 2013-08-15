#include "predicter.h"
#include "model.h"
#include "imagecoordinates.h"
#include "beamevaluator.h"
#include "solutionfile.h"

void Predicter::applyGain(double *dataVal, const std::complex<double> *gain)
{
  double gainSq[4];
  gainSq[0] = std::fabs(gain[0] * gain[0] + gain[1] * gain[1]);
  gainSq[1] = std::fabs(gain[0] * gain[2] + gain[1] * gain[3]);
  gainSq[2] = std::fabs(gain[2] * gain[0] + gain[3] * gain[1]);
  gainSq[3] = std::fabs(gain[2] * gain[2] + gain[3] * gain[3]);

  dataVal[0] = dataVal[0] * gainSq[0] + dataVal[1] * gainSq[2];
  dataVal[1] = dataVal[0] * gainSq[1] + dataVal[1] * gainSq[3];
  dataVal[2] = dataVal[2] * gainSq[0] + dataVal[3] * gainSq[2];
  dataVal[3] = dataVal[2] * gainSq[1] + dataVal[3] * gainSq[3];
}

void Predicter::Initialize(ModelSource& source, BeamEvaluator *beamEvaluator)
{
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
			beamEvaluator->EvaluateGain(source.PosRA(), source.PosDec(), centreFreq, &parameters->beamValues[ch*4]);
		}
		else {
			for(size_t p=0; p!=4; ++p)
				parameters->beamValues[ch*4+p] = 0.0;
		}
		for(size_t p=0; p!=4; ++p)
			_totalFlux[p] += parameters->brightness[ch*4+p];
	}
	parameters->lmsqrt = sqrt(1.0 - l*l - m*m);
	
	source.SetUserData(parameters);
}

void Predicter::Initialize(Model& model, const std::string& solutionFile, BeamEvaluator *beamEvaluator)
{
	for(Model::iterator i=model.begin(); i!=model.end(); ++i)
		Initialize(*i, beamEvaluator);
	if(!solutionFile.empty())
		readSolutions(solutionFile);
}

void Predicter::ReportSources(Model& model)
{
	std::cout << "Model predicter initialized with " << model.SourceCount() << " sources of apparent brightness [";
	std::cout << (_totalFlux[0] / _channelCount);
	for(size_t p=1; p!=4; ++p)
		std::cout << ',' << (_totalFlux[p] / _channelCount);
	std::cout << "]\n";
	
	std::cout << "(absolute brightness: [";
	std::cout << model.TotalFlux(_startFrequency, _endFrequency, 0);
	for(size_t p=1; p!=4; ++p)
		std::cout << ',' << model.TotalFlux(_startFrequency, _endFrequency, p);
	std::cout << "], app at low freq: " << model.TotalFlux(_startFrequency, 0) << ", app at high freq: " << model.TotalFlux(_endFrequency, 0) << ")\n";
}

Predicter::CNumType Predicter::Predict(const ModelSource& source, NumType u, NumType v, NumType w, size_t channelIndex, size_t polarizationIndex)
{
	switch(source.Type())
	{
		case ModelSource::PointSource:
		{
			SourceParameters *parameters = reinterpret_cast<SourceParameters *>(source.UserData());
			NumType l = parameters->l, m = parameters->m;
			NumType lmsqrt = parameters->lmsqrt;
			NumType angle = 2.0*M_PI*(u*l + v*m + w*(lmsqrt-1.0));
			NumType fact = parameters->brightness[channelIndex*4+polarizationIndex] / lmsqrt;
			double sinangle, cosangle;
			sincos(angle, &sinangle, &cosangle);
			return CNumType(fact * cosangle, fact * sinangle);
		}
	}
	return 0.0;
}

Predicter::CNumType Predicter::Predict(const Model& model, NumType u, NumType v, NumType w, size_t channelIndex, size_t polarizationIndex)
{
	CNumType sum(0.0, 0.0);
	for(Model::const_iterator i=model.begin(); i!=model.end(); ++i)
		sum += Predict(*i, u, v, w, channelIndex, polarizationIndex);
	return sum;
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
		Matrix2x2::ATimesB(temp, &parameters->beamValues[_channelCount*4], dest);
		Matrix2x2::ATimesHermB(dest, temp, &parameters->beamValues[_channelCount*4]);
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
