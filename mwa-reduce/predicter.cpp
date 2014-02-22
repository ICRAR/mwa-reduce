#include "predicter.h"
#include "model.h"
#include "imagecoordinates.h"
#include "beamevaluator.h"
#include "solutionfile.h"

void Predicter::initialize(ModelComponent& component)
{
	SourceParameters *parameters = new SourceParameters();
	component.SetUserData(parameters);
	NumType l, m;
	ImageCoordinates::RaDecToLM<NumType>(component.PosRA(), component.PosDec(), _ra0, _dec0, l, m);
	parameters->l = l;
	parameters->m = m;
	parameters->brightness = new CNumType[_channelCount*4];
	parameters->beamValues = new CNumType[_channelCount*4];
	parameters->appBrightness = new CNumType[_channelCount*4];
	for(size_t ch=0;ch!=_channelCount;++ch)
	{
		NumType stokesValues[4];
		
		for(size_t p=0; p!=4; ++p)
		{
			stokesValues[p] = component.SED().FluxAtChannel(ch, _channelCount, _startFrequency, _endFrequency, Polarization::IndexToStokes(p));
		}
		Polarization::StokesToLinear(stokesValues, &parameters->brightness[ch*4]);
	}

	updateBeam(component);
	
	for(size_t ch=0;ch!=_channelCount;++ch)
	{
		CNumType temp[4];
		Matrix2x2::ATimesB(temp, &parameters->beamValues[ch*4], &parameters->brightness[ch*4]);
		Matrix2x2::PlusATimesHermB(_totalFlux, temp, &parameters->beamValues[ch*4]);
	}
	parameters->lmsqrt = sqrt(1.0 - l*l - m*m);
}

void Predicter::Initialize(ModelSource& source, BeamEvaluator *beamEvaluator)
{
	_beamEvaluator = beamEvaluator;
	for(ModelSource::iterator i=source.begin(); i!=source.end(); ++i)
		initialize(*i);
}

void Predicter::Initialize(Model& model, const std::string& solutionFile, BeamEvaluator *beamEvaluator)
{
	_beamEvaluator = beamEvaluator;
	
	for(Model::iterator src=model.begin(); src!=model.end(); ++src)
	{
		for(ModelSource::iterator i=src->begin(); i!=src->end(); ++i)
			initialize(*i);
	}
	if(!solutionFile.empty())
		readSolutions(solutionFile);
}

void Predicter::updateBeam(ModelComponent& component)
{
	SourceParameters *parameters = reinterpret_cast<SourceParameters *>(component.UserData());
	BeamEvaluator::PrecalcPosInfo posInfo;
	if(_beamEvaluator != 0) {
		_beamEvaluator->PrecalculatePositionInfo(posInfo, component.PosRA(), component.PosDec());
	}
	for(size_t ch=0;ch!=_channelCount;++ch)
	{
		if(_beamEvaluator != 0)
		{
			double chCentreFreq = _startFrequency + (long double) ch * (_endFrequency - _startFrequency) / (long double) (_channelCount-1);
			_beamEvaluator->EvaluateAbsToApparentGain(posInfo, chCentreFreq, &parameters->beamValues[ch*4]);
			CNumType temp[4];
			Matrix2x2::ATimesB(temp, &parameters->beamValues[ch*4], &parameters->brightness[ch*4]);
			Matrix2x2::ATimesHermB(&parameters->appBrightness[ch*4], temp, &parameters->beamValues[ch*4]);
		}
		else {
			parameters->beamValues[ch*4+0] = 1.0; parameters->beamValues[ch*4+1] = 0.0;
			parameters->beamValues[ch*4+2] = 0.0; parameters->beamValues[ch*4+3] = 1.0;
			parameters->appBrightness[ch*4+0] = 1.0; parameters->beamValues[ch*4+1] = 0.0;
			parameters->appBrightness[ch*4+2] = 0.0; parameters->beamValues[ch*4+3] = 1.0;
		}
	}
}

void Predicter::UpdateBeam(Model& model)
{
	for(Model::iterator src=model.begin(); src!=model.end(); ++src)
	{
		for(ModelSource::iterator i=src->begin(); i!=src->end(); ++i)
			updateBeam(*i);
	}
}

void Predicter::ReportSources(Model& model)
{
	std::cout << "Model predicter initialized with " << model.SourceCount() << " sources of apparent brightness [";
	std::cout << (_totalFlux[0].real() / _channelCount);
	for(size_t p=1; p!=4; ++p)
		std::cout << ',' << (_totalFlux[p].real() / _channelCount);
	std::cout << "]\n";
	
	std::cout << "(absolute brightness: [";
	std::cout << model.TotalFlux(_startFrequency, _endFrequency, Polarization::StokesI);
	std::cout << ',' << model.TotalFlux(_startFrequency, _endFrequency, Polarization::StokesQ);
	std::cout << ',' << model.TotalFlux(_startFrequency, _endFrequency, Polarization::StokesU);
	std::cout << ',' << model.TotalFlux(_startFrequency, _endFrequency, Polarization::StokesV);
	std::cout << "], abs at low freq: " << model.TotalFlux(_startFrequency, Polarization::StokesI) << ", abs at high freq: " << model.TotalFlux(_endFrequency, Polarization::StokesI) << ")\n";
}

void Predicter::predict4(CNumType *dest, const ModelComponent& component, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2)
{
	switch(component.Type())
	{
		case ModelComponent::PointSource:
		{
			SourceParameters *parameters = reinterpret_cast<SourceParameters *>(component.UserData());
			NumType l = parameters->l, m = parameters->m;
			NumType lmsqrt = parameters->lmsqrt;
			NumType angle = 2.0*M_PI*(u*l + v*m + w*(lmsqrt-1.0));
			double sinangleOverLMS, cosangleOverLMS;
			sincos(angle, &sinangleOverLMS, &cosangleOverLMS);
			sinangleOverLMS /= lmsqrt;
			cosangleOverLMS /= lmsqrt;
			if(_beamEvaluator != 0)
			{
				for(size_t p=0; p!=4; ++p)
				{
					CNumType fact = parameters->appBrightness[channelIndex*4+p];
					dest[p] = CNumType(
						fact.real() * cosangleOverLMS - fact.imag() * sinangleOverLMS,
						fact.real() * sinangleOverLMS + fact.imag() * cosangleOverLMS);
				}
			}
			else {
				for(size_t p=0; p!=4; ++p)
				{
					CNumType fact = parameters->brightness[channelIndex*4+p];
					dest[p] = CNumType(fact.real() * cosangleOverLMS - fact.imag() * sinangleOverLMS,
														 fact.real() * sinangleOverLMS + fact.imag() * cosangleOverLMS);
				}
			}
		}
	}
	
	if(!_rhsSolutions.empty())
	{
		CNumType temp[4];
		std::complex<double>
			*antenna1Sol = &_rhsSolutions[a1*_channelCount*4],
			*antenna2Sol = &_rhsSolutions[a2*_channelCount*4];
		Matrix2x2::ATimesB(temp, antenna1Sol, dest);
		Matrix2x2::ATimesHermB(dest, temp, antenna2Sol);
	}
	/*if(_beamEvaluator != 0)
	{
		SourceParameters *parameters = reinterpret_cast<SourceParameters *>(component.UserData());
		Matrix2x2::ATimesB(temp, &parameters->beamValues[channelIndex*4], dest);
		Matrix2x2::ATimesHermB(dest, temp, &parameters->beamValues[channelIndex*4]);
	}*/
}

void Predicter::Predict4(Predicter::CNumType* dest, const ModelSource& source, Predicter::NumType u, Predicter::NumType v, Predicter::NumType w, size_t channelIndex, size_t a1, size_t a2)
{
	for(size_t p=0; p!=4; ++p)
		dest[p] = 0.0;
	for(ModelSource::const_iterator i=source.begin(); i!=source.end(); ++i)
	{
		CNumType temp[4];
		predict4(dest, *i, u, v, w, channelIndex, a1, a2);
		for(size_t p=0; p!=4; ++p)
			dest[p] += temp[p];
	}
}

void Predicter::Predict4(CNumType *dest, const Model& model, NumType u, NumType v, NumType w, size_t channelIndex, size_t a1, size_t a2)
{
	for(size_t p=0; p!=4; ++p)
		dest[p] = 0.0;
	for(Model::const_iterator i=model.begin(); i!=model.end(); ++i)
	{
		for(ModelSource::const_iterator j=i->begin(); j!=i->end(); ++j)
		{
			CNumType temp[4];
			predict4(temp, *j, u, v, w, channelIndex, a1, a2);
			for(size_t p=0; p!=4; ++p)
				dest[p] += temp[p];
		}
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
