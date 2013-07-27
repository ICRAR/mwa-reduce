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
		double ra = sourceIter->PosRA();
		double dec = sourceIter->PosDec();
		SpectralEnergyDistribution &sed = sourceIter->SED();
		
		double gains[4];
		beamEval.EvaluateGain(ra, dec, gains);
		sed.SetConstantBeam(gains[0], gains[1], gains[2], gains[3]);
	}
	
	model.Save(argv[3]);
}
