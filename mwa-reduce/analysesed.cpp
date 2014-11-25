#include <iostream>

#include "model.h"
#include "sedanalyser.h"

int main(int argc, char* argv[])
{
	bool withRM = false;
	
	SEDAnalyser analyser;
	
	std::cout << "Reading model... " << std::flush;
	const char* filename = argv[1];
	analyser.Read(filename);
	std::cout << "DONE\n";
	
	std::cout << "Determining centre frequency... " << std::flush;
	long double centreFrequency = analyser.CentralFrequency();
	std::cout << centreFrequency*1e-6 << " MHz.\n";
	
	std::cout << "Fitting power laws... " << std::flush;
	analyser.Process();
	std::cout << "DONE\n";
	
	if(withRM)
		analyser.ProcessRM();
	
	analyser.SaveResults();
}
