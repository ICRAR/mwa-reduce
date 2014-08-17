#include "mrccatalogue.h"

struct MrcEntry
{
	// 1 - 10
	char name[9];
	char undefined1;
	
	// 11 - 21
	char raHoursB1950[2];
	char undefined2;
	char raMinutesB1950[2];
	char undefined3;
	char raSecsB1950[4];
	char undefined4;
	
	// 22 - 31
	char declinationSignB1950;
	char decDegB1950[2];
	char undefined5;
	char decMinB1950[2];
	char undefined6;
	char decSecB1950[2];
	char undefined7;
	
	// 32 - 42
	char raHours[2];
	char undefined8;
	char raMins[2];
	char undefined9;
	char raSec[4];
	char undefined10;
	
	// 43 - 51
	char decSign;
	char decDeg[2];
	char undefined11;
	char decMin[2];
	char undefined12;
	char decSec[2];
	
	// 52 - 58
	char stdErrs[7];
	
	// 59 - 65
	char fluxDensity[7];
	
	// 66 - 70
	char fluxDensityStdErr[4];
	char scatteringFlag;
	
	// 71 - 76
	char morphologyFlags[3];
	char crossReferenceFlags[3];
	
	// 77 - 87
	char galacticLongitude[5];
	char galacticLattitude[6];
	char endLineSymbol;
};

MRCCatalogue::MRCCatalogue(const std::string& filename) :
	_filename(filename),
	_file(filename.c_str())
{
}

bool MRCCatalogue::ReadNext(ModelSource& source)
{
	source = ModelSource();
	if(!_file.good()) return false;
	MrcEntry entry;
	_file.read(reinterpret_cast<char*>(&entry), sizeof(entry));
	
	ModelComponent component;
	int raHrs = (entry.raHours[0]-'0')*10 + (entry.raHours[1]-'0');
	int raMins = (entry.raMins[0]-'0')*10 + (entry.raMins[1]-'0');
	long double raSec = (entry.raSec[0]-'0')*10 + (entry.raSec[1]-'0') +
		(entry.raSec[3]-'0')*0.1;
	component.SetPosRA(raHrs*(M_PI/12.0) + raMins*(M_PI/(12.0*60.0)) + raSec*(M_PI/(12.0*60.0*60.0)));
	
	long double decSign = (entry.decSign == '-') ? -1.0 : 1.0;
	int decHrs = (entry.decDeg[0]-'0')*10 + (entry.decDeg[1]-'0');
	int decMins = (entry.decMin[0]-'0')*10 + (entry.decMin[1]-'0');
	int decSec = (entry.decSec[0]-'0')*10 + (entry.decSec[1]-'0');
	component.SetPosDec((decHrs*(M_PI/180.0) + decMins*(M_PI/(180.0*60.0)) + decSec*(M_PI/(180.0*60.0*60.0))) * decSign);
	
	long double fluxDensity =
		digToVal(entry.fluxDensity[0])*1000.0 +
		digToVal(entry.fluxDensity[1])*100.0 +
		digToVal(entry.fluxDensity[2])*10.0 +
		digToVal(entry.fluxDensity[3])*1.0 +
		digToVal(entry.fluxDensity[5])*0.1 +
		digToVal(entry.fluxDensity[6])*0.01;
	component.SetSED(SpectralEnergyDistribution(fluxDensity, 408000000.0));
	
	char nameStr[11];
	memcpy(nameStr, entry.name, 10);
	nameStr[9] = 0;
	if(nameStr[8] == 32) nameStr[8] = 0;
	
	source.SetName(nameStr);
	source.AddComponent(component);
	
	return true;
}
