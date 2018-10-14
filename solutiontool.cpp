#include <iostream>

#include "solutionfile.h"

int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		std::cout << "Syntax: solutiontool <solutions.bin> <antenna>\n";
	}
	else {
		const char* filename(argv[1]);
		size_t selectedAnt = atoi(argv[2]);
		SolutionFile file;
		file.OpenForReading(filename);
		size_t
			nAnt = file.AntennaCount(),
		//	nInt = file.IntervalCount(),
			nChanBlock = file.ChannelCount();
		if(nAnt < selectedAnt)
			throw std::runtime_error("Antenna not found");
		
		for(size_t a = 0; a!=selectedAnt; ++a) {
			for(size_t cb = 0; cb!=nChanBlock; ++cb) {
				for(size_t p = 0; p!=4; ++p) {
					file.ReadNextSolution();
				}
			}
		}
		
		for(size_t cb = 0; cb!=nChanBlock; ++cb) {
			std::cout << cb;
			for(size_t p = 0; p!=4; ++p) {
				std::cout << '\t' << file.ReadNextSolution();
			}
			std::cout << '\n';
		}
	}
}

