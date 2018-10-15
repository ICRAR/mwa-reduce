#include "predicter.h"
#include "beamevaluator.h"
#include "solutionfile.h"

#include "units/imagecoordinates.h"

#include "model/model.h"

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
	
	if(component.Type() == ModelComponent::GaussianSource)
	{
		// Using the FWHM formula for a Gaussian:
		double sigmaMaj = component.MajorAxis() / (2.0L * sqrtl(2.0L * logl(2.0L)));
		double sigmaMin = component.MinorAxis() / (2.0L * sqrtl(2.0L * logl(2.0L)));
		// Position angle is angle from North:
		double paSin, paCos;
		sincos(component.PositionAngle()+0.5*M_PI, &paSin, &paCos);
		// Make rotation matrix
		long double transf[4];
		transf[0] = paCos;
		transf[1] = -paSin;
		transf[2] = paSin;
		transf[3] = paCos;
		// Multiply with scaling matrix to make variance 1.
		// sigmamaj/min are multiplications and include pi^2 factor, because the sigma
		// of the Fourier transform of a Gaus is 1/sigma of the normal Gaus and has a sqrt(2 pi^2) factor.
		parameters->gausTransf[0] = transf[0] * sigmaMaj * M_PI * sqrt(2.0);
		parameters->gausTransf[1] = transf[1] * sigmaMaj * M_PI * sqrt(2.0);
		parameters->gausTransf[2] = transf[2] * sigmaMin * M_PI * sqrt(2.0);
		parameters->gausTransf[3] = transf[3] * sigmaMin * M_PI * sqrt(2.0);
	}
	
	for(size_t ch=0;ch!=_channelCount;++ch)
	{
		NumType stokesValues[4];
		
		for(size_t p=0; p!=4; ++p)
		{
			stokesValues[p] = component.SED().FluxAtChannel(ch, _channelCount, _startFrequency, _endFrequency, Polarization::IndexToStokes(p));
		}
		Polarization::StokesToLinear(stokesValues, &parameters->brightness[ch*4]);
	}

	updateBeam(component, 0, _channelCount);
	
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

void Predicter::Initialize(Model& model, BeamEvaluator *beamEvaluator)
{
	_beamEvaluator = beamEvaluator;
	
	for(Model::iterator src=model.begin(); src!=model.end(); ++src)
	{
		for(ModelSource::iterator i=src->begin(); i!=src->end(); ++i)
			initialize(*i);
	}
}

void Predicter::updateBeam(ModelComponent& component, size_t startChannel, size_t endChannel)
{
	SourceParameters *parameters = reinterpret_cast<SourceParameters *>(component.UserData());
	BeamEvaluator::PrecalcPosInfo posInfo;
	if(_beamEvaluator != nullptr) {
		_beamEvaluator->PrecalculatePositionInfo(posInfo, component.PosRA(), component.PosDec());
	}
	
	_threads.For(startChannel, endChannel, [&](size_t ch, size_t /*thread*/)
	{
		if(_beamEvaluator != nullptr)
		{
			double chCentreFreq = (_channelCount>1) ? (_startFrequency + (long double) ch * (_endFrequency - _startFrequency) / (long double) (_channelCount-1)) : _startFrequency;
			_beamEvaluator->EvaluateAbsToApparentGain(posInfo, chCentreFreq, &parameters->beamValues[ch*4]);
			CNumType temp[4];
			Matrix2x2::ATimesB(temp, &parameters->beamValues[ch*4], &parameters->brightness[ch*4]);
			Matrix2x2::ATimesHermB(&parameters->appBrightness[ch*4], temp, &parameters->beamValues[ch*4]);
		}
		else {
			parameters->beamValues[ch*4+0] = 1.0; parameters->beamValues[ch*4+1] = 0.0;
			parameters->beamValues[ch*4+2] = 0.0; parameters->beamValues[ch*4+3] = 1.0;
			parameters->appBrightness[ch*4+0] = 1.0; parameters->appBrightness[ch*4+1] = 0.0;
			parameters->appBrightness[ch*4+2] = 0.0; parameters->appBrightness[ch*4+3] = 1.0;
		}
	});
}

void Predicter::UpdateBeam(Model& model, size_t startChannel, size_t endChannel)
{
	for(Model::iterator src=model.begin(); src!=model.end(); ++src)
	{
		for(ModelSource::iterator i=src->begin(); i!=src->end(); ++i)
			updateBeam(*i, startChannel, endChannel);
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
		case ModelComponent::GaussianSource:
		{
			SourceParameters *parameters = reinterpret_cast<SourceParameters *>(component.UserData());
			NumType l = parameters->l, m = parameters->m;
			NumType lmsqrt = parameters->lmsqrt;
			NumType angle = 2.0*M_PI*(u*l + v*m + w*(lmsqrt-1.0));
			double sinangleOverLMS, cosangleOverLMS;
			sincos(angle, &sinangleOverLMS, &cosangleOverLMS);
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
	if(component.Type() == ModelComponent::GaussianSource)
	{
		SourceParameters& parameters = *reinterpret_cast<SourceParameters*>(component.UserData());
		
		NumType
			uTransf = u*parameters.gausTransf[0] + v*parameters.gausTransf[1],
			vTransf = u*parameters.gausTransf[2] + v*parameters.gausTransf[3];
		//std::cout << uTransf << ',' << vTransf << '\n';
		NumType g = exp(-uTransf*uTransf - vTransf*vTransf);
		for(size_t p=0; p!=4; ++p)
			dest[p] *= g;
	}
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
