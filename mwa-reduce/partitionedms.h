#ifndef PARTITIONED_MS
#define PARTITIONED_MS

#include <string>
#include <vector>

#include <ms/MeasurementSets/MeasurementSet.h>

class PartitionedMS
{
public:
	PartitionedMS(const std::string& filename);
	
	static void Partition(std::vector<std::string>& result, const std::string& source, const std::string& dataColumn, bool includeModel);
	
	void OpenMS(casa::MeasurementSet* dest);
	
	void ReadData();
	
	void ReadModel();
	
	
};

#endif
