#include <casacore/ms/MeasurementSets/MeasurementSet.h>

using namespace casacore;

int main(int argc, char *argv[])
{
	// Create a MeasurementSet and cast it to a Table pointer
	// Table is the rootclass of MeasurementSet
	Table *ms = new MeasurementSet(argv[1]);
	
	// Delete the measurement set. This will only call ~Table(),
	// and not ~MeasurementSet(), because ~Table() is non-virtual.
	delete ms;
	
	// This leaks a lot of stuff!
	return 0;
}
