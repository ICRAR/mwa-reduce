#include "predicter.h"
#include "model.h"
#include "imagecoordinates.h"

void Predicter::Initialize(ModelSource& source)
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
				source.SED().FluxAtFrequency(ch, _channelCount, _startFrequency, _endFrequency, p);
		}
	}
	parameters->lmsqrt = sqrt(1.0 - l*l - m*m);
	
	source.SetUserData(parameters);
}

void Predicter::Initialize(Model& model)
{
	for(Model::iterator i=model.begin(); i!=model.end(); ++i)
		Initialize(*i);
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
