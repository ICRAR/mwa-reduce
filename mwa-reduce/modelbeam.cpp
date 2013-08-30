#include "model.h"
#include "tilebeam.h"
#include "banddata.h"
#include "beamevaluator.h"

#include <ms/MeasurementSets/MeasurementSet.h>

#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include <measures/TableMeasures/ScalarMeasColumn.h>
#include <measures/Measures/MPosition.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MDirection.h>

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		std::cout << "Usage: modelbeam <ms> <model-in> <model-out>\n";
		return -1;
	}
	
	casa::MeasurementSet ms(argv[1]);
	BeamEvaluator beamEval(ms);
		
	Model model(argv[2]);
	std::cout << "Calculating beam for " << model.SourceCount() << " sources...\n";
	for(Model::iterator sourceIter=model.begin(); sourceIter!=model.end(); ++sourceIter)
	{
		for(ModelSource::iterator compIter=sourceIter->begin(); compIter!=sourceIter->end(); ++compIter)
		{
			double ra = compIter->PosRA();
			double dec = compIter->PosDec();
			SpectralEnergyDistribution &sed = compIter->SED();
			
			std::complex<double> gains[4];
			beamEval.EvaluateAbsToApparentGain(ra, dec, gains);
			sed.SetConstantBeam(std::abs(gains[0]), std::abs(gains[1]), std::abs(gains[2]), std::abs(gains[3]));
		}
	}
	
	model.Save(argv[3]);
}
