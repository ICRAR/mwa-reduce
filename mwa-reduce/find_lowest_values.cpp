#include <iostream>
#include <set>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <tables/Tables/ArrayColumn.h>
#include <tables/Tables/ScalarColumn.h>

#include "progressbar.h"

struct Visibility
{
	bool operator<(const Visibility& other) const { return amplitudeSq < other.amplitudeSq; }
	
	double amplitudeSq;
	casa::Complex value;
	size_t timestep;
	size_t a1, a2;
	size_t channel;
	size_t polarization;
};

int main(int argc, char* argv[])
{
	casa::MeasurementSet ms(argv[1]);
	casa::ROArrayColumn<casa::Complex> dataCol(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::DATA));
	casa::ROScalarColumn<int>
		ant1Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA1)),
		ant2Col(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::ANTENNA2));
	casa::ROScalarColumn<double>
		timeCol(ms, casa::MeasurementSet::columnName(casa::MSMainEnums::TIME));
	
	casa::IPosition shape = dataCol.shape(0);
	casa::Array<casa::Complex> dataArray(shape);
	size_t polarizationCount = shape(0);
	size_t channelCount = shape(1);
	ProgressBar progress("Searching");
	
	size_t count = atoi(argv[2]);
	
	std::set<Visibility> lowest;
	double curTime = timeCol(0);
	size_t timestep = 0;
	for(size_t row=0; row!=ms.nrow(); ++row)
	{
		progress.SetProgress(row, ms.nrow());
		if(timeCol(row) != curTime)
		{
			++timestep;
			curTime = timeCol(row);
		}
		dataCol.get(row, dataArray);
		Visibility vis;
		vis.timestep = timestep;
		vis.a1 = ant1Col(row);
		vis.a2 = ant2Col(row);
		casa::Array<casa::Complex>::const_contiter iter = dataArray.cbegin();
		for(size_t ch=0; ch!=channelCount; ++ch)
		{
			for(size_t p=0; p!=polarizationCount; ++p)
			{
				double r = iter->real(), im = iter->imag();
				double amplitudeSq = r*r + im*im;
				if(amplitudeSq != 0.0)
				{
					if(lowest.size() < count || lowest.rbegin()->amplitudeSq > amplitudeSq)
					{
						vis.amplitudeSq = amplitudeSq;
						vis.value = *iter;
						vis.channel = ch;
						vis.polarization = p;
						lowest.insert(vis);
					}
					if(lowest.size() > count)
					{
						std::set<Visibility>::iterator last = lowest.end();
						--last;
						lowest.erase(last);
					}
				}
				++iter;
			}
		}
	}
	progress.SetProgress(ms.nrow(), ms.nrow());
	
	size_t index = 0;
	for(std::set<Visibility>::const_iterator i=lowest.begin(); i!=lowest.end(); ++i)
	{
		std::cout <<
			"Index: " << index << '\n' <<
			"Value: " << sqrt(i->amplitudeSq) << " (" << i->value.real() << " + " << i->value.imag() << "i)\n" <<
			"Timestep: " << i->timestep << '\n' <<
			"Antenna1: " << i->a1 << '\n' <<
			"Antenna2: " << i->a2 << '\n' <<
			"Channel: " << i->channel << '\n' <<
			"Polarization: " << i->polarization << "\n\n";
		++index;
	}
	
	return 0;
}
