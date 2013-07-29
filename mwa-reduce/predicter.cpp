#include "predicter.h"
#include "model.h"
#include "imagecoordinates.h"
#include "beamevaluator.h"

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
			beamEvaluator->AbsToApparent(source.PosRA(), source.PosDec(), centreFreq, &parameters->brightness[ch*4]);
		}
		for(size_t p=0; p!=4; ++p)
			_totalFlux[p] += parameters->brightness[ch*4+p];
	}
	parameters->lmsqrt = sqrt(1.0 - l*l - m*m);
	
	source.SetUserData(parameters);
}

void Predicter::Initialize(Model& model, BeamEvaluator *beamEvaluator)
{
	for(Model::iterator i=model.begin(); i!=model.end(); ++i)
		Initialize(*i, beamEvaluator);
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
	std::cout << "])\n";
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
